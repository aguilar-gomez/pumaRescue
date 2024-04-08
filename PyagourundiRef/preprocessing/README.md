1. **CallAnnotate**: individual gvcfs
  2. joint gvcfs
  3. call genotypes
  4. left-align, trim
  5. annotate with snpEff
6. **RepeatMask**:
- Mask repeats that were soft masked by NCBI, using WindowMasker
- Exclude Sex chromosomes and mitochondria
- RepeatMasker
7. **customFiltering**: Jacqueline Robinson's script
  - Sites where REF is in [A,C,G,T] and ALT is in [A,C,G,T,.] go on to genotype filtering if AD and DP present in FORMAT
  - Filtered out genotypes are changed to './.', all others reported
  - Genotypes failing filters are set to missing (./.)
  - Applies individual min and max depth filters
  - Filters heterozygotes if the allele balance (REF/DP) is <20% or >80%
  - Filters homozygotes if more than 10% of the alleles are different type
  - Fail sites with poor variant metrics:
    - 'QD' (Quality by Depth)
    - 'FS' (Fisher Strand bias)
    - 'MQ' (Mapping Quality)
    - 'SOR' (Strand Odds Ratio)
8. **SimplifyExtract**: remove sites that fail filters and merge all scaffolds
