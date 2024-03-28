1. **CallAnnotate**: call genotypes, left-align, trim and annotate with snpEff
2. **customFiltering**: Jacqueline Robinson's script
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
