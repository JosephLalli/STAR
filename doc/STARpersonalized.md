# Exact fast construction and optional resident personalization

The normal STAR index-and-align commands remain available. To select exact libsais
construction, add `--genomeGenerateMethod libsais` to `--runMode genomeGenerate`.
The resulting index uses STAR's ordinary files and loader. STAR's native constructor
remains the default. Set `--limitGenomeGenerateRAM` for the memory available to the
build; it remains binding for the constructor allocation estimate.

For one-off mapping from an available FASTA, GTF and reads, the resident mode avoids
writing and loading that index:

```sh
STAR --runMode genomeGenerateAndAlign --genomeGenerateMethod libsais \
     --runThreadN 32 --limitGenomeGenerateRAM 400000000000 \
     --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 149 \
     --readFilesIn reads_1.fastq reads_2.fastq --twopassMode Basic \
     --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts \
     --outFileNamePrefix sample.
```

The memory value is an example allocation budget, not a portable memory requirement
or an operating-system hard cap. Size it for the machine and concurrent work. The
resident mode requires a Full genome, `NoSharedMemory`, no output coordinate transform,
and `sjdbInsertSave` other than `All`. It does not serialize a reusable index. To retain
one, use the normal generation command and then the normal alignment command.
Fast construction requires dense suffixes (`genomeSAsparseD=1`), unrestricted suffix
length, a 64-bit little-endian platform and at most 125 contigs. Other configurations
can use the native constructor. Supplied personalized FASTA/GTF files need no regeneration.

## Optional reference-plus-variants route

Add both options to the resident command:

```sh
--personalizationVcf calls.vcf.gz --personalizationSample donor
```

Use coordinate-sorted, indexed BGZF VCF or BCF. Variant contig names and REF alleles
must match the supplied reference FASTA. No external consensus or GTF liftover step is
required: allele application and forward coordinate chains remain resident, and GTF
rows are projected as they are consumed. This route duplicates contigs into `_L` and
`_R` and applies the same suffix to gene, transcript, exon and protein identifiers.

Allele policy:

- All FILTER states are considered; this is not a PASS-only selection.
- A diploid genotype's written GT slots select its two alleles. Homozygous calls apply
  to both; a haploid genotype supplies the same allele to both unless its contig is
  explicitly emitted once.
- Unphased heterozygotes keep written GT order and are counted in the log. This does
  not infer phase or establish biological haplotype identity.
- Missing, reference and symbolic/star alleles leave the relevant sequence unchanged.
  Incompatible overlap and past-end calls are counted and skipped. The embedded
  bcftools application retains its shared-anchor indel rules.
- Concrete REF mismatch, malformed genotypes, unsupported ploidy and unsorted input
  are fatal rather than silently changing the reference.

Contig policies are optional explicit lists: `--personalizationHaploidContigs` emits
those contigs once; `--personalizationSkipAnnotationContigs` retains their sequence
while omitting annotation; `--personalizationExcludeContigs` omits their sequence and
annotation. By default, no contig receives one of these policies. They do not apply
an automatic mappability mask.

## Optional exports

`--personalizationOutputPrefix prefix` writes
`<prefix><donor>.personalized.fa.gz` and `.gtf.gz`, plus FASTA `.fai` and `.gzi` indexes.
`--personalizationTranscriptFasta Yes` additionally writes
`<prefix><donor>.personalized.transcripts.fa.gz`, with records in transcriptome-BAM
header order, plus `.fai` and `.gzi`. Transcript export requires a GTF and output prefix.
These exports are off by default; requested background writers finish before success
is reported. The resulting BAM and transcript-reference identifiers must be kept
consistent when configuring downstream quantification.

Linux builds and the recorded small-fixture comparisons were validated. Other
platforms and optional HTSlib configurations need their corresponding validation.
