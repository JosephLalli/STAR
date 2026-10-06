/* Thin resident adapter to pinned upstream bcftools consensus.c (MIT).
 * Allele selection is performed by PersonalizedInput; upstream applies the
 * selected allele and constructs its exact forward chain. No CLI is invoked. */
#include "ConsensusSymbols.h"
#define error ce_upstream_fatal
#define main_consensus ce_unused_cli
#include "upstream/consensus.c"
#undef main_consensus
#undef error
#include "ConsensusEmbedded.h"
#include <limits.h>

void ce_upstream_fatal(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    fputs("STAR embedded bcftools consensus: ", stderr);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(3);
}

struct ConsensusEmbedded {
    args_t args;
    bcf1_t *record;
    int reference_length, finished, nblocks;
    int *starts, *ends, *query_starts;
    char last_base;  /* last base of the previous edited fragment */
};

ConsensusEmbedded *ce_create(const char *contig, int reference_length)
{
    if (!contig || !*contig || reference_length <= 0) return NULL;
    ConsensusEmbedded *c = calloc(1, sizeof(*c));
    if (!c) return NULL;
    c->reference_length = reference_length;
    c->args.hdr = bcf_hdr_init("w");
    c->record = bcf_init();
    c->args.chr = strdup(contig);
    if (!c->args.hdr || !c->record || !c->args.chr) { ce_destroy(c); return NULL; }
    kstring_t line = {0,0,NULL};
    ksprintf(&line,"##contig=<ID=%s,length=%d>",contig,reference_length);
    if (bcf_hdr_append(c->args.hdr,line.s) || bcf_hdr_sync(c->args.hdr)) {
        free(line.s); ce_destroy(c); return NULL;
    }
    free(line.s);
    c->args.rid = bcf_hdr_name2id(c->args.hdr,contig);
    c->args.isample = -1;
    c->args.fa_frz_pos = -1;
    c->args.prev_base_pos = -1;
    c->args.fa_length = reference_length;
    c->args.chain = init_chain(NULL,0);
    return c;
}

int ce_apply(ConsensusEmbedded *c, int pos, const char *ref, const char *alt,
             const char **result, int *length)
{
    if (!c || c->finished || !ref || !alt || !result || !length) return -1;
    const size_t nr = strlen(ref), na = strlen(alt);
    args_t *a = &c->args;
    const int64_t query_pos = (int64_t)pos + a->fa_mod_off;
    const int64_t query_length = (int64_t)c->reference_length + a->fa_mod_off + (int64_t)na - nr;
    if (!nr || !na || nr > INT_MAX || na > INT_MAX || pos < 0 ||
        (int64_t)pos + nr > c->reference_length || pos < a->fa_frz_pos ||
        query_pos < 0 || query_pos > INT_MAX || query_length <= 0 || query_length > INT_MAX)
        return -1;
    if (strspn(ref,"ACGTN") != nr || strspn(alt,"ACGTN") != na) return -1;
    bcf_clear(c->record);
    c->record->rid = a->rid;
    c->record->pos = pos;
    const char *alleles[2] = {ref,alt};
    if (bcf_update_alleles(a->hdr,c->record,alleles,2) < 0) return -1;
    bcf_unpack(c->record,BCF_UN_STR);
    /* A call on the frozen position shares its first base with the previous fragment, which the caller
     * has emitted: upstream sees that edited base, as in the CLI's running buffer, and decides whether the
     * call is an anchored indel it applies; the returned fragment starts after the shared base. */
    const int shared = pos == a->fa_frz_pos;
    a->fa_ori_pos = (int)query_pos;
    a->fa_buf.l = 0;
    if (shared ? kputc(c->last_base,&a->fa_buf) < 0 || kputsn(ref+1,nr-1,&a->fa_buf) < 0
               : kputsn(ref,nr,&a->fa_buf) < 0) return -1;
    const int before = a->napplied;
    apply_variant(a,c->record);
    if (a->napplied != before+1 || a->fa_buf.l != na) return -1;
    c->last_base = a->fa_buf.s[na-1];
    *result = a->fa_buf.s + shared;
    *length = (int)a->fa_buf.l - shared;
    return 0;
}

void ce_finish(ConsensusEmbedded *c)
{
    if (!c || c->finished) return;
    const chain_t *chain = c->args.chain;
    const int n = chain->num + 1;
    c->starts = malloc((size_t)n*sizeof(int));
    c->ends = malloc((size_t)n*sizeof(int));
    c->query_starts = malloc((size_t)n*sizeof(int));
    if (!c->starts || !c->ends || !c->query_starts) ce_upstream_fatal("out of memory retaining chain\n");
    int r = 0, q = 0;
    for (int i=0; i<n; ++i) {
        int length = i<chain->num ? chain->block_lengths[i] : c->reference_length-chain->ref_last_block_ori;
        if (length>0) {
            int j=c->nblocks++;
            c->starts[j]=r; c->ends[j]=r+length; c->query_starts[j]=q;
        }
        r+=length; q+=length;
        if (i<chain->num) { r+=chain->ref_gaps[i]; q+=chain->alt_gaps[i]; }
    }
    if (r != c->reference_length || (int64_t)q != (int64_t)c->reference_length+c->args.fa_mod_off)
        ce_upstream_fatal("inconsistent consensus chain\n");
    c->finished=1;
}

int ce_block_count(const ConsensusEmbedded *c) { return c && c->finished ? c->nblocks : -1; }
void ce_block(const ConsensusEmbedded *c, int i, int *r, int *e, int *q)
{
    if (!c || !c->finished || i<0 || i>=c->nblocks) ce_upstream_fatal("invalid chain block\n");
    *r=c->starts[i]; *e=c->ends[i]; *q=c->query_starts[i];
}
void ce_destroy(ConsensusEmbedded *c)
{
    if (!c) return;
    free(c->args.fa_buf.s); free(c->args.chr); destroy_chain(c->args.chain);
    if (c->args.hdr) bcf_hdr_destroy(c->args.hdr);
    if (c->record) bcf_destroy(c->record);
    free(c->starts); free(c->ends); free(c->query_starts); free(c);
}
