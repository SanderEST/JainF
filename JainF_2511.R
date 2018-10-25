# Preprocessing
# 
#
# bcftools to split multiallelic and leftalign indels (recommended for Annovar)
#
# > prefix=JainF_2511
# > ref=/Users/sander/data_references/Homo_sapiens_assembly19.fasta
# > bcftools norm -m-both -o ${prefix}.step1.vcf ${prefix}.raw.vcf.gz
# > bcftools norm -f ${ref} -o ${prefix}.normalized.vcf ${prefix}.step1.vcf
# 
# Annovar to annotate
# > perl ~/data_tools/annovar/table_annovar.pl ${wd}${prefix}.normalized.vcf ~/data_tools/annovar/humandb/ \
# > -buildver hg19 -out ${wd}${family} -remove -protocol refGene,gnomad_genome,gnomad_exome,exac03,dbnsfp33a,clinvar_20180603 \
# > -operation g,f,f,f,f,f -arg '-hgvs,,,,,,' -nastring . -otherinfo -vcfinput


# Load libraries
library(tidyverse)
library(vcfR)

#Read input

vcf <- read.vcfR("JainF_2511.hg19_multianno.vcf", verbose = FALSE )
genes <- read_tsv('genes.tsv', col_names = c('gene', 'disorder', 'inheritance'))

##Specify info_columns

types <- c(AC= 'i',
  AF= 'n',
  AN= 'i',
  BaseQRankSum= 'n',
  ClippingRankSum= 'n',
  DB= 'n',
  DP= 'n',
  DS= 'n',
  END= 'i',
  ExcessHet= 'n',
  FS= 'n',
  HaplotypeScore= 'n',
  InbreedingCoeff= 'n',
  MLEAC= 'n',
  MLEAF= 'n',
  MQ= 'n',
  MQRankSum= 'n',
  QD= 'n',
  RAW_MQ= 'n',
  ReadPosRankSum= 'n',
  SOR= 'n',
  gnomAD_genome_ALL= 'n',
  gnomAD_genome_AFR= 'n',
  gnomAD_genome_AMR= 'n',
  gnomAD_genome_ASJ= 'n',
  gnomAD_genome_EAS= 'n',
  gnomAD_genome_FIN= 'n',
  gnomAD_genome_NFE= 'n',
  gnomAD_genome_OTH= 'n',
  gnomAD_exome_ALL= 'n',
  gnomAD_exome_AFR= 'n',
  gnomAD_exome_AMR= 'n',
  gnomAD_exome_ASJ= 'n',
  gnomAD_exome_EAS= 'n',
  gnomAD_exome_FIN= 'n',
  gnomAD_exome_NFE= 'n',
  gnomAD_exome_OTH= 'n',
  gnomAD_exome_SAS= 'n',
  ExAC_ALL= 'n',
  ExAC_AFR= 'n',
  ExAC_AMR= 'n',
  ExAC_EAS= 'n',
  ExAC_FIN= 'n',
  ExAC_NFE= 'n',
  ExAC_OTH= 'n',
  ExAC_SAS= 'n',
  SIFT_score= 'n',
  SIFT_converted_rankscore= 'n',
  Polyphen2_HDIV_score= 'n',
  Polyphen2_HDIV_rankscore= 'n',
  Polyphen2_HVAR_score= 'n',
  Polyphen2_HVAR_rankscore= 'n',
  LRT_score= 'n',
  LRT_converted_rankscore= 'n',
  MutationTaster_score= 'n',
  MutationTaster_converted_rankscore= 'n',
  MutationAssessor_score= 'n',
  MutationAssessor_score_rankscore= 'n',
  FATHMM_score= 'n',
  FATHMM_converted_rankscore= 'n',
  PROVEAN_score= 'n',
  PROVEAN_converted_rankscore= 'n',
  VEST3_score= 'n',
  VEST3_rankscore= 'n',
  MetaSVM_score= 'n',
  MetaSVM_rankscore= 'n',
  MetaLR_score= 'n',
  MetaLR_rankscore= 'n',
  'M-CAP_score'= 'n',
  'M-CAP_rankscore'= 'n',
  CADD_raw= 'n',
  CADD_raw_rankscore= 'n',
  CADD_phred= 'n',
  DANN_score= 'n',
  DANN_rankscore= 'n',
  'fathmm-MKL_coding_score'= 'n',
  'fathmm-MKL_coding_rankscore'= 'n',
  Eigen_coding_or_noncoding= 'n',
  'Eigen-raw'= 'n',
  'Eigen-PC-raw'= 'n',
  GenoCanyon_score= 'n',
  GenoCanyon_score_rankscore= 'n',
  integrated_fitCons_score= 'n',
  integrated_fitCons_score_rankscore= 'n',
  integrated_confidence_value= 'n',
  'GERP++_RS'= 'n',
  'GERP++_RS_rankscore'= 'n',
  phyloP100way_vertebrate= 'n',
  phyloP100way_vertebrate_rankscore= 'n',
  phyloP20way_mammalian= 'n',
  phyloP20way_mammalian_rankscore= 'n',
  phastCons100way_vertebrate= 'n',
  phastCons100way_vertebrate_rankscore= 'n',
  phastCons20way_mammalian= 'n',
  phastCons20way_mammalian_rankscore= 'n',
  SiPhy_29way_logOdds= 'n',
  SiPhy_29way_logOdds_rankscore= 'n')

# VcfR format to tidy dataframes

df <- vcfR2tidy(vcf, dot_is_NA = T, info_types = types)

vars<- df[[1]]
gts <- df[[2]]
meta <- df[[3]]

#Clean up the workspace
rm(vcf)
rm(df)

#Define max allowed AF in population databases to be considered as pathogenic.
popAF_lim = 0.001

#Filter for possible pathogenic mutations (filter in exonic, conventional splice sites, filter out common and synonymous) 

core_vars <- vars %>% filter(gnomAD_genome_ALL < popAF_lim | is.na(gnomAD_genome_ALL),
                gnomAD_exome_ALL < popAF_lim | is.na(gnomAD_exome_ALL),
                ExAC_ALL < popAF_lim | is.na(ExAC_ALL),
                Func.refGene == 'exonic' | Func.refGene == 'splicing' | Func.refGene == 'exonic\\x3bsplicing',
                ExonicFunc.refGene != 'synonymous_SNV')


# Annotate variants with pathogenicity estimations
lof = c("frameshift_deletion", "frameshift_insertion", "stopgain")
splicing = c('splicing','exonic\\x3bsplicing')

core_vars_class <- core_vars %>% mutate(Classification = case_when( 
  ExonicFunc.refGene %in% lof | Func.refGene %in% splicing | CLNSIG == "Pathogenic" ~ 'Pathogenic',
  CLNSIG %in% c("Likely_pathogenic", "Pathogenic/Likely_pathogenic")  ~ 'Likely pathogenic',
  CLNSIG %in% c("Likely_benign") ~ 'Likely benign',
  CLNSIG %in% c("Benign", "Benign/Likely_benign", "Benign/Likely_benign,_other") ~ 'Benign',
  TRUE ~ 'VUS')
  ) 

# Select cols you want to keep
tojoin <- core_vars_class %>% select(ChromKey, POS, Gene.refGene, Classification)


# Annotate genes with info
tojoin <- full_join(tojoin, genes, by = c('Gene.refGene' = 'gene'))


# Join genotypes table with variants

gts_add <- inner_join(gts, tojoin)

# Transform genotype to variant count - 2 for homozygous/hemizygous and 1 for heterozygous
gts_add <- gts_add %>% mutate(vars = case_when(
  gt_GT == '1/1' ~ 2L,
  gt_GT == '0/1' | gt_GT == '1/0' ~ 1L,
  gt_GT == '0/0' ~ 0L,
  TRUE ~ 0L
))

#select AR/XL genes and AD genes, and remove benign/likely benign
### AR/AD genes and DG genes are left out at this point!

ar_pat <- gts_add %>% filter(inheritance %in% c('AR', 'XL'), Classification %in% c('Likely pathogenic', 'Pathogenic'))
ad_pat <- gts_add %>% filter(inheritance %in% c('AD'), Classification %in% c('Likely pathogenic', 'Pathogenic'))

#Summarize to sum detected pathogenic and likely pathogenic variants. Select for possibly diagnostic findings (AR >= 2, AD >= 1)

ar_summary <- ar_pat %>% group_by(Indiv, Gene.refGene, Classification) %>%
  summarise(var_count = sum(vars)) %>% 
  group_by(Indiv, Gene.refGene) %>% 
  summarise(pathogenic_vars = sum(var_count)) %>%
  filter(pathogenic_vars >= 2)

ad_summary  <- ad_pat %>% group_by(Indiv, Gene.refGene, Classification) %>%
  summarise(var_count = sum(vars)) %>% 
  group_by(Indiv, Gene.refGene) %>% 
  summarise(pathogenic_vars = sum(var_count)) %>%
  filter(pathogenic_vars >= 1)

# Output per sample var_counts for AR and AD

var_counts <- gts_add %>% group_by(Indiv, Gene.refGene, Classification) %>%
  summarise(var_count = sum(vars)) %>% filter(var_count != 0)

write_tsv(var_counts, 'JainF_2511_var_counts.tsv')

# make data frame containing all samples and all genes 2511 x 34

samples <- unique(gts$Indiv)
gene_list <- unlist(genes['gene'])

gene <- rep(gene_list, length(samples))
sample_id <- rep(samples, each = length(gene_list))

all_samples_genes <- data_frame(sample_id,gene)

# Prepare path_var and vus_var counts
path_vars <- gts_add %>% group_by(Indiv, Gene.refGene, Classification) %>% 
  summarise(var_count = sum(vars)) %>% 
  filter(Classification %in% c('Likely pathogenic', 'Pathogenic')) %>% 
  group_by(Indiv, Gene.refGene) %>% 
  summarise(pathogenic_vars = sum(var_count))

vus_vars <- gts_add %>% group_by(Indiv, Gene.refGene, Classification) %>% 
  summarise(vus_vars = sum(vars)) %>% 
  filter(Classification %in% c('VUS')) %>%
  select(-Classification)

# Join the counts

all_samples_genes <- full_join(all_samples_genes, path_vars, by= c("sample_id" = "Indiv", "gene" = "Gene.refGene")) %>% full_join(vus_vars, by= c("sample_id" = "Indiv", "gene" = "Gene.refGene"))

# Add the legend

all_samples_genes <- all_samples_genes %>% mutate(Legend = case_when(
  (pathogenic_vars == 0 | is.na(pathogenic_vars)) & (vus_vars == 0 | is.na(vus_vars)) ~ 'NONE',
  (pathogenic_vars == 0 | is.na(pathogenic_vars)) & (vus_vars == 1 ) ~ 'VUS',
  (pathogenic_vars == 0 | is.na(pathogenic_vars)) & (vus_vars >= 2 ) ~ 'BI_VUS',
  (pathogenic_vars == 1) & (vus_vars == 0 | is.na(pathogenic_vars)) ~ 'PATH',
  (pathogenic_vars == 1) & (vus_vars >= 1 ) ~ 'VUS_PATH',
  (pathogenic_vars >= 2) ~ 'BI_PATH'
))

# Plot

p <- ggplot(all_samples_genes, aes(gene, sample_id)) + geom_tile(aes(fill = Legend,width=0.9,height=0.9),colour = "grey95")
p <- p + scale_x_discrete("") + scale_y_discrete("Jain Foundation Cohort, n = 2511")
p <- p + scale_fill_manual(values=c("BI_VUS" = rgb(0,0,1,alpha=1.0),"VUS" = rgb(0,0,1,alpha=0.2), "NONE" = "white","BI_PATH" = "red","VUS_PATH" = "purple","PATH" = rgb(1,0,0,alpha=0.2)))
p <- p + theme(axis.text.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
p
