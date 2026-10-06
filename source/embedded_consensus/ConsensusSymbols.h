/* bcftools carries a separate region-index implementation from HTSlib.
 * Keep its symbols private to the embedded consensus translation units. */
#define regidx_set ce_regidx_set
#define regidx_parse_bed ce_regidx_parse_bed
#define regidx_parse_tab ce_regidx_parse_tab
#define regidx_parse_reg ce_regidx_parse_reg
#define regidx_parse_vcf ce_regidx_parse_vcf
#define regidx_init ce_regidx_init
#define regidx_init_string ce_regidx_init_string
#define regidx_destroy ce_regidx_destroy
#define regidx_overlap ce_regidx_overlap
#define regidx_insert ce_regidx_insert
#define regidx_insert_list ce_regidx_insert_list
#define regidx_push ce_regidx_push
#define regidx_seq_names ce_regidx_seq_names
#define regidx_seq_nregs ce_regidx_seq_nregs
#define regidx_nregs ce_regidx_nregs
#define regitr_init ce_regitr_init
#define regitr_destroy ce_regitr_destroy
#define regitr_reset ce_regitr_reset
#define regitr_overlap ce_regitr_overlap
#define regitr_loop ce_regitr_loop
#define regitr_copy ce_regitr_copy
