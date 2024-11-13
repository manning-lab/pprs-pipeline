# 1000G files no IDs, so we fill the IDs with chr:pos:ref:alt.
Rscript ../pprs.R "
  --geno_files http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/*.vcf.gz$
  --fill_vcf_ids_with chr:pos:ref:alt

  --score_file ../score_files/alzheimers.csv

  --score_file_chr_col Chr
  --score_file_pos_col Pos
  --score_file_id_col chr_cpra_id
  --score_file_ref_col Ref
  --score_file_alt_col Alt
  --score_file_ea_col Effect_Allele

  --output_fnm alzheimers-1000G_high_coverage.txt
"
# --threads 4
# --ldlink_token <YOUR TOKEN>
# --ldlink_pop GBR
