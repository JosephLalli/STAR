#ifndef STAR_CONSENSUS_EMBEDDED_H
#define STAR_CONSENSUS_EMBEDDED_H
#ifdef __cplusplus
extern "C" {
#endif
typedef struct ConsensusEmbedded ConsensusEmbedded;
ConsensusEmbedded *ce_create(const char *contig, int reference_length);
/* Sorted concrete biallelic calls; 0 means success. A call may start on the previous call's last REF
 * base only as an indel upstream applies there; its result then omits that shared base.
 * result is borrowed until the next apply/destroy; it is length-delimited. */
int ce_apply(ConsensusEmbedded *, int pos0, const char *ref, const char *alt,
             const char **result, int *length);
void ce_finish(ConsensusEmbedded *);
int ce_block_count(const ConsensusEmbedded *);
void ce_block(const ConsensusEmbedded *, int i, int *ref_start, int *ref_end, int *query_start);
void ce_destroy(ConsensusEmbedded *);
#ifdef __cplusplus
}
#endif
#endif
