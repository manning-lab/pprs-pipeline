Rscript ../pprs.R "
  --geno_files http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/*.vcf.gz

  --score_file score_files/alzheimers.csv

  --score_file_chr_col Chr
  --score_file_pos_col Pos
  --score_file_ref_col Ref
  --score_file_alt_col Alt
  --score_file_ea_col Effect_Allele
  --score_file_weight_cols Apolipoprotein_B_Pos CRP_Neg Diffusivity_inferior_longitudinal_fasciculus_Neg ICVF_Cingulum_CingulateR Proinsulin_Neg T2stat_PutamenL_Neg Total_GRS

  --output_fnm alzheimers-1000G.txt
"

# Other args you could add:
# --threads 4 
# --ldlink_token <YOUR TOKEN>
# --ldlink_pop GBR
