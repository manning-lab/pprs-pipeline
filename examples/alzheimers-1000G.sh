# 1000 Genomes files have no IDs, so we fill the IDs with chr:pos:ref:alt.
Rscript ../pprs.R "
  --geno_files http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/*.vcf.gz
  --fill_vcf_ids_with chr:pos:ref:alt

  --score_file ../score_files/alzheimers.csv

  --score_file_chr_col Chr
  --score_file_pos_col Pos
  --score_file_id_col chr_cpra_id
  --score_file_ref_col Ref
  --score_file_alt_col Alt
  --score_file_ea_col Effect_Allele
  --score_file_weight_cols 10 Total_GRS 12 13 14 15 11

  --output_fnm alzheimers-1000G.txt
"
# --threads 4 
# --ldlink_token <YOUR TOKEN>
# --ldlink_pop GBR
