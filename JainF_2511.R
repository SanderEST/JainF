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
# > -buildver hg19 -out ${wd}${prefix} -remove -protocol refGene,gnomad_genome,gnomad_exome,exac03,dbnsfp33a,clinvar_20180603 \
# > -operation g,f,f,f,f,f -arg '-hgvs,,,,,,' -nastring . -otherinfo -vcfinput
#
# Compress 
# bgzip JainF_2511.hg19_multianno.vcf
# tabix JainF_2511.hg19_multianno.vcf.gz


# Load libraries
library(tidyverse)
library(vcfR)
library(scales)

#Read input

vcf <- read.vcfR("JainF_2511.hg19_multianno.vcf.gz", verbose = FALSE )
genes <- read_tsv('genes.tsv', col_names = c('gene', 'disorder', 'inheritance'))

#Extract indels and snps

indels <- extract.indels(vcf, return.indels = T)

snps <- extract.indels(vcf, return.indels = F)


##Specify info_columns

types <- c(AC= 'i', AF= 'n', AN= 'i', BaseQRankSum= 'n', ClippingRankSum= 'n', DB= 'n', DP= 'n', 
           DS= 'n', END= 'i', ExcessHet= 'n', FS= 'n', HaplotypeScore= 'n', InbreedingCoeff= 'n', 
           MLEAC= 'n', MLEAF= 'n', MQ= 'n', MQRankSum= 'n', QD= 'n', RAW_MQ= 'n', ReadPosRankSum= 'n', 
           SOR= 'n', gnomAD_genome_ALL= 'n', gnomAD_genome_AFR= 'n', gnomAD_genome_AMR= 'n', 
           gnomAD_genome_ASJ= 'n', gnomAD_genome_EAS= 'n', gnomAD_genome_FIN= 'n', gnomAD_genome_NFE= 'n', 
           gnomAD_genome_OTH= 'n', gnomAD_exome_ALL= 'n', gnomAD_exome_AFR= 'n', gnomAD_exome_AMR= 'n', 
           gnomAD_exome_ASJ= 'n', gnomAD_exome_EAS= 'n', gnomAD_exome_FIN= 'n', gnomAD_exome_NFE= 'n', 
           gnomAD_exome_OTH= 'n', gnomAD_exome_SAS= 'n', ExAC_ALL= 'n', ExAC_AFR= 'n', ExAC_AMR= 'n', 
           ExAC_EAS= 'n', ExAC_FIN= 'n', ExAC_NFE= 'n', ExAC_OTH= 'n', ExAC_SAS= 'n', SIFT_score= 'n', 
           SIFT_converted_rankscore= 'n', Polyphen2_HDIV_score= 'n', Polyphen2_HDIV_rankscore= 'n', 
           Polyphen2_HVAR_score= 'n', Polyphen2_HVAR_rankscore= 'n', LRT_score= 'n', LRT_converted_rankscore= 'n', 
           MutationTaster_score= 'n', MutationTaster_converted_rankscore= 'n', MutationAssessor_score= 'n', 
           MutationAssessor_score_rankscore= 'n', FATHMM_score= 'n', FATHMM_converted_rankscore= 'n', 
           PROVEAN_score= 'n', PROVEAN_converted_rankscore= 'n', VEST3_score= 'n', VEST3_rankscore= 'n', 
           MetaSVM_score= 'n', MetaSVM_rankscore= 'n', MetaLR_score= 'n', MetaLR_rankscore= 'n', 
           'M-CAP_score'= 'n', 'M-CAP_rankscore'= 'n', CADD_raw= 'n', CADD_raw_rankscore= 'n', CADD_phred= 'n', 
           DANN_score= 'n', DANN_rankscore= 'n', 'fathmm-MKL_coding_score'= 'n', 'fathmm-MKL_coding_rankscore'= 'n', 
           Eigen_coding_or_noncoding= 'n', 'Eigen-raw'= 'n', 'Eigen-PC-raw'= 'n', GenoCanyon_score= 'n', 
           GenoCanyon_score_rankscore= 'n', integrated_fitCons_score= 'n', integrated_fitCons_score_rankscore= 'n', 
           integrated_confidence_value= 'n', 'GERP++_RS'= 'n', 'GERP++_RS_rankscore'= 'n', phyloP100way_vertebrate= 'n', 
           phyloP100way_vertebrate_rankscore= 'n', phyloP20way_mammalian= 'n', phyloP20way_mammalian_rankscore= 'n', 
           phastCons100way_vertebrate= 'n', phastCons100way_vertebrate_rankscore= 'n', phastCons20way_mammalian= 'n', 
           phastCons20way_mammalian_rankscore= 'n', SiPhy_29way_logOdds= 'n', SiPhy_29way_logOdds_rankscore= 'n')

# VcfR format to tidy dataframes

df <- vcfR2tidy(vcf, dot_is_NA = T, info_types = types)

vars<- df[[1]]
gts <- df[[2]]
meta <- df[[3]]

#index varints
vars$index <- 1:dim(vars)[1]
gts$index <- rep(1:dim(vars)[1], 2511)


#Clean up the workspace
rm(vcf)
rm(df)

#Define max allowed AF in population databases to be considered as pathogenic.
popAF_lim = 0.005

#Filter for possible pathogenic mutations (filter in exonic, conventional splice sites, filter out common and synonymous) 

core_vars <- vars %>% filter(gnomAD_genome_ALL < popAF_lim | is.na(gnomAD_genome_ALL),
                gnomAD_exome_ALL < popAF_lim | is.na(gnomAD_exome_ALL),
                ExAC_ALL < popAF_lim | is.na(ExAC_ALL))


# Annotate variants with pathogenicity estimations
lof = c("frameshift_deletion", "frameshift_insertion", "stopgain")
splicing = c('splicing','exonic\\x3bsplicing')

core_vars_class <- core_vars %>% mutate(Classification = case_when( 
  ExonicFunc.refGene %in% lof | Func.refGene %in% splicing | CLNSIG == "Pathogenic" ~ 'Pathogenic',
  CLNSIG %in% c("Likely_pathogenic", "Pathogenic/Likely_pathogenic")  ~ 'Likely pathogenic',
  CLNSIG %in% c("Likely_benign") ~ 'Likely benign',
  CLNSIG %in% c("Benign", "Benign/Likely_benign", "Benign/Likely_benign,_other") ~ 'Benign',
  ExonicFunc.refGene == 'synonymous_SNV' ~ 'Likely benign',
  !(Func.refGene %in% c('exonic', 'splicing','exonic\\x3bsplicing')) ~ 'Likely benign',
  TRUE ~ 'VUS')
  ) 

core_vars_class$Classification[which(core_vars_class$POS == 47259533)] <- 'Pathogenic'

# Select cols you want to keep
tojoin <- core_vars_class %>% select(ChromKey, POS, Gene.refGene, Classification, index)


# Annotate genes with info
tojoin <- full_join(tojoin, genes, by = c('Gene.refGene' = 'gene'))


# Join genotypes table with variants

gts_add <- inner_join(gts, tojoin)

#filter for quality

gts_add <- gts_add %>% separate(gt_AD, c("AD_REF", 'AD_ALT'), sep = ',', remove = F) %>% 
  mutate(AD_ALT = as.integer(AD_ALT), AD_REF = as.integer(AD_REF)) %>% 
  mutate(altvar_pct = case_when(AD_ALT > 0 ~ AD_ALT / (AD_ALT + AD_REF) *100, TRUE ~ 0))

gts_add <- gts_add %>% filter(gt_DP >= 5, altvar_pct > 30)

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

#write_tsv(var_counts, 'JainF_2511_var_counts.tsv')

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

all_samples_genes <- full_join(all_samples_genes, path_vars, by= c("sample_id" = "Indiv", "gene" = "Gene.refGene")) %>% 
                      full_join(vus_vars, by= c("sample_id" = "Indiv", "gene" = "Gene.refGene"))

#Replace NAs with 0.
all_samples_genes <- all_samples_genes %>% mutate_all(funs(replace(., is.na(.), 0)))

# Add the legend and gene data

all_samples_genes <- all_samples_genes %>% mutate(Legend = case_when(
  (pathogenic_vars == 0 ) & (vus_vars == 0 ) ~ 'NONE',
  (pathogenic_vars == 0 ) & (vus_vars == 1 ) ~ 'VUS',
  (pathogenic_vars == 0 ) & (vus_vars >= 2 ) ~ 'BI_VUS',
  (pathogenic_vars == 1) & (vus_vars == 0 ) ~ 'PATH',
  (pathogenic_vars == 1) & (vus_vars >= 1 ) ~ 'VUS_PATH',
  (pathogenic_vars >= 2) ~ 'BI_PATH'
))

all_samples_genes <- all_samples_genes %>% full_join(genes)

# Separate recessive genes

all_samples_ARgenes <- all_samples_genes %>% filter(inheritance == 'AR')

all_samples_ARgenes$gene <- factor(all_samples_ARgenes$gene, levels= rev(c("CAPN3", "DYSF", "SGCG",  "SGCA",  "SGCB",  "SGCD", "TCAP", "TRIM32", 
                                                                           "FKRP", "TTN", "POMT1", "ANO5",  "FKTN",  "POMT2", "POMGNT1", "DAG1", "PLEC", "ISPD", "GAA", "GNE"))) 


path_vus_summary <- all_samples_ARgenes %>% group_by(gene) %>% summarise(path_vus_sum = sum(pathogenic_vars, na.rm=TRUE) + sum(vus_vars, na.rm=TRUE), 
                                                                         path_sum = sum(pathogenic_vars, na.rm=TRUE), 
                                                                         vus_sum = sum(vus_vars, na.rm=TRUE) ) %>% arrange(desc(path_vus_sum))

all_samples_ARgenes$gene2 <- factor(all_samples_ARgenes$gene, levels= rev(as.character(path_vus_summary$gene))) 

# Add diagnosis_status 

all_samples_genes <- all_samples_genes %>% mutate(diagnosis = case_when(
  inheritance == 'AD' & pathogenic_vars >= 1 ~ 'AD_diagnosis',
  inheritance == 'AD' & vus_vars >= 1 ~ 'AD_vus',
  inheritance %in% c('AR', 'XL') & pathogenic_vars >= 2 ~ 'AR_diagnosis',
  inheritance %in% c('AR', 'XL') & (pathogenic_vars == 1 & vus_vars >= 1) ~ 'AR_pathogenic_and_vus',
  inheritance %in% c('AR', 'XL') & (pathogenic_vars == 1 & vus_vars == 0) ~ 'AR_single_pathogenic',
  inheritance %in% c('AR', 'XL') & (pathogenic_vars == 0 & vus_vars >= 2) ~ 'AR_bi_vus',
  inheritance %in% c('AR', 'XL') & (pathogenic_vars == 0 & vus_vars == 1) ~ 'AR_single_vus',
  inheritance == 'AD/AR' & pathogenic_vars >= 2 ~ 'AD_AR_diagnosis',
  inheritance == 'AD/AR' & (pathogenic_vars == 1 & vus_vars >= 1) ~ 'AD_AR_pathogenic_and_vus',
  inheritance == 'AD/AR' & (pathogenic_vars == 1 & vus_vars == 0) ~ 'AD_AR_single_pathogenic',
  inheritance == 'AD/AR' & (pathogenic_vars == 0 & vus_vars >= 2) ~ 'AD_AR_bi_vus',
  inheritance == 'AD/AR' & (pathogenic_vars == 0 & vus_vars == 1) ~ 'AD_AR_single_vus',
  TRUE ~ 'no_significant_findings'
))

diagnosis_df <- all_samples_genes %>% select(sample_id, gene, diagnosis) %>% 
  spread(key = gene, value = diagnosis)

#count the per gene conclusions and spread the data to get a per sample table
dx_count <- all_samples_genes %>% count(sample_id, diagnosis) %>% spread(diagnosis, n)

dx_count <- dx_count %>% mutate(conclusion = case_when(
  AD_diagnosis >= 2 | AR_diagnosis >= 2 | AD_AR_diagnosis >= 2 | 
    (AD_diagnosis == 1 & AR_diagnosis == 1) |
    (AD_AR_diagnosis == 1 & AR_diagnosis == 1) |
    (AD_diagnosis == 1 & AD_AR_diagnosis == 1) |
    (AD_diagnosis == 1 & AR_diagnosis == 1 & AD_AR_diagnosis == 1)~ 'Double_dx',
  AD_diagnosis == 1 ~ 'AD_dx',
  AR_diagnosis  == 1 | AD_AR_diagnosis >= 1 ~ 'AR_dx',
  AR_pathogenic_and_vus >= 1 | AD_AR_pathogenic_and_vus >= 1 ~ 'AR_pathogenic_and_vus',
  AR_single_pathogenic >= 1 | AD_AR_single_pathogenic >= 1 ~ 'AR_single_pathogenic',
  AR_bi_vus >= 1 | AD_vus >= 1 | AD_AR_bi_vus >= 1 ~ 'AR_bi_vus/AD_vus',
  AR_single_vus >= 1 | AD_AR_single_vus >= 1 ~ 'AR_single_vus',
  TRUE ~ 'No_findings'
))

dx_c <- dx_count %>% count(conclusion)

dx_c$conclusion <- factor(dx_c$conclusion, levels = c("Double_dx", "AD_dx", "AR_dx", "AR_pathogenic_and_vus", "AR_single_pathogenic", "AR_bi_vus/AD_vus", "AR_single_vus", "No_findings"))

#AR one pathogenic hit

ar1hit <- dx_count %>% filter(conclusion == "AR_single_pathogenic") %>% select(sample_id, AR_single_pathogenic, AD_AR_single_pathogenic, conclusion)

ardf <- all_samples_genes %>% select(sample_id, gene, diagnosis) %>% inner_join(ar1hit) %>% filter(diagnosis %in% c('AR_single_pathogenic', 'AD_AR_single_pathogenic'))

ar1hitgenes <- ardf %>% count(gene)

ar1hitgts <- gts_add %>% right_join(ar1hit, by = c('Indiv' = 'sample_id')) %>% filter(vars > 0, Classification %in% c('Pathogenic', 'Likely pathogenic')) 

ar1hitgts_i <- ar1hitgts %>% count(index) %>% arrange(desc(n))

ar1hitvars <- core_vars_class %>% inner_join(ar1hitgts_i)


write_excel_csv(ar1hitvars, 'AR_1hit_vars_Jainf2511.xls')

##
arpathvus <- dx_count %>% filter(conclusion == "AR_pathogenic_and_vus") %>% select(sample_id, AR_pathogenic_and_vus, AD_AR_pathogenic_and_vus, conclusion)

arpathvusdf <- all_samples_genes %>% select(sample_id, gene, diagnosis) %>% inner_join(arpathvus) %>% filter(diagnosis %in% c('AD_AR_pathogenic_and_vus', 'AR_pathogenic_and_vus'))

arpathvusgenes <- arpathvusdf %>% count(gene)





# Sanity check

neg <- read_tsv('Negative_samples.tsv')

neg <- rename(neg, sample_id = 'Sample')


conclusions <- dx_count %>% select(sample_id, conclusion) %>% separate(sample_id, c('sample_id', 'index'), sep = '_')

neg_conclusions <- neg %>% left_join(conclusions)

#neg_conclusions %>% count(conclusion)

fnegs <- conclusions %>% anti_join(neg) %>% filter(conclusion == 'No_findings')


# Plots

p <- ggplot(all_samples_genes, aes(gene, sample_id)) + geom_tile(aes(fill = Legend,width=0.9,height=0.9),colour = "grey95")
p <- p + scale_x_discrete("") + scale_y_discrete("Jain Foundation Cohort, n = 2511")
p <- p + scale_fill_manual(values=c("BI_VUS" = rgb(0,0,1,alpha=1.0),"VUS" = rgb(0,0,1,alpha=0.2), "NONE" = "white","BI_PATH" = "red","VUS_PATH" = "purple","PATH" = rgb(1,0,0,alpha=0.2)))
p <- p + theme(axis.text.y=element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
p

r <- ggplot(all_samples_ARgenes, aes(sample_id, gene2)) + geom_tile(aes(fill = Legend,width=0.9,height=0.9),colour = "grey100", size = 0.001)
r <- r + scale_y_discrete("") + scale_x_discrete("Jain Foundation Cohort, n = 2511")
r <- r + scale_fill_manual(values=c("BI_VUS" = rgb(0,0,1,alpha=1.0),"VUS" = rgb(0,0,1,alpha=0.2), "NONE" = "white","BI_PATH" = "red","VUS_PATH" = "purple","PATH" = "orange"))
r <- r + theme(axis.text.x=element_blank())
r

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pie <- ggplot(dx_c, aes(x="", y=n, fill=conclusion)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start=0, direction = -1) 
pie + blank_theme +
  theme(axis.text.x=element_blank()) + scale_fill_brewer(palette = 'YlOrRd', direction = -1)



