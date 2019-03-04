library(readxl)
results_tab_D <- read_excel("G:/My Drive/research/shogun_metagenome/ARG/resfinder/results_tab_D.xlsx")
D_abricate_tab <- read_excel("G:/My Drive/research/shogun_metagenome/ARG/resfinder/D_abricate.tab.xlsx")
  
#results_tab_D <- read_excel("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/ARG/resfinder/results_tab_D.xlsx")
#D_abricate_tab <- read_excel("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/ARG/resfinder/D_abricate.tab.xlsx")
library(dplyr)
library(tidyr)
D <- merge(results_tab_D, D_abricate_tab, by="Contig")
Df <- full_join(results_tab_D, D_abricate_tab, by="Contig")

#setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/MEGAHIT_D.contigs.FASTA/abricate")
setwd("G:/My Drive/research/shogun_metagenome/MEGAHIT_D.contigs.FASTA/abricate")

megahit_D_argannot <- read_excel("megahit_D.argannot.xlsx")
colnames(megahit_D_argannot) <- c("#FILE_argannot","SEQUENCE","START_argannot","END_argannot","GENE_argannot","COVERAGE_argannot","COVERAGE_MAP_argannot","GAPS_argannot","%COVERAGE_argannot","%IDENTITY_argannot","DATABASE_argannot","ACCESSION_argannot","PRODUCT_argannot")
megahit_D_card <- read_excel("megahit_D.card.xlsx")
colnames(megahit_D_card) <- c("#FILE_card","SEQUENCE","START_card","END_card","GENE_card","COVERAGE_card","COVERAGE_MAP_card","GAPS_card","%COVERAGE_card","%IDENTITY_card","DATABASE_card","ACCESSION_card","PRODUCT_card")
megahit_D_ncbi <- read_excel("megahit_D.ncbi.xlsx")
colnames(megahit_D_ncbi) <- c("#FILE_ncbi","SEQUENCE","START_ncbi","END_ncbi","GENE_ncbi","COVERAGE_ncbi","COVERAGE_MAP_ncbi","GAPS_ncbi","%COVERAGE_ncbi","%IDENTITY_ncbi","DATABASE_ncbi","ACCESSION_ncbi","PRODUCT_ncbi")
megahit_D_resfinder <- read_excel("megahit_D.resfinder.xlsx")
colnames(megahit_D_resfinder) <- c("#FILE_resfinder","SEQUENCE","START_resfinder","END_resfinder","GENE_resfinder","COVERAGE_resfinder","COVERAGE_MAP_resfinder","GAPS_resfinder","%COVERAGE_resfinder","%IDENTITY_resfinder","DATABASE_resfinder","ACCESSION_resfinder","PRODUCT_resfinder")

a1 <- full_join(megahit_D_argannot,megahit_D_card,by="SEQUENCE")
a2 <- full_join(a1,megahit_D_ncbi,by="SEQUENCE")
a3 <- full_join(a2,megahit_D_resfinder,by="SEQUENCE")
#write.csv(a3,"combine_D.csv",row.names = F)

#setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C")
setwd("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_C")
megahit_C_gene_argannot <- read_excel("megahit_C_gene.argannot.xlsx")
colnames(megahit_C_gene_argannot) <- c("#FILE_argannot","SEQUENCE","START_argannot","END_argannot","GENE_argannot","COVERAGE_argannot","COVERAGE_MAP_argannot","GAPS_argannot","%COVERAGE_argannot","%IDENTITY_argannot","DATABASE_argannot","ACCESSION_argannot","PRODUCT_argannot")
megahit_C_gene_card <- read_excel("megahit_C_gene.card.xlsx")
colnames(megahit_C_gene_card) <- c("#FILE_card","SEQUENCE","START_card","END_card","GENE_card","COVERAGE_card","COVERAGE_MAP_card","GAPS_card","%COVERAGE_card","%IDENTITY_card","DATABASE_card","DNA Accession","PRODUCT_card")
megahit_C_gene_ncbi <- read_excel("megahit_C_gene.ncbi.xlsx")
colnames(megahit_C_gene_ncbi) <- c("#FILE_ncbi","SEQUENCE","START_ncbi","END_ncbi","GENE_ncbi","COVERAGE_ncbi","COVERAGE_MAP_ncbi","GAPS_ncbi","%COVERAGE_ncbi","%IDENTITY_ncbi","DATABASE_ncbi","ACCESSION_ncbi","PRODUCT_ncbi")
megahit_C_gene_resfinder <- read_excel("megahit_C_gene.resfinder.xlsx")
colnames(megahit_C_gene_resfinder) <- c("#FILE_resfinder","SEQUENCE","START_resfinder","END_resfinder","GENE_resfinder","COVERAGE_resfinder","COVERAGE_MAP_resfinder","GAPS_resfinder","%COVERAGE_resfinder","%IDENTITY_resfinder","DATABASE_resfinder","ACCESSION_resfinder","PRODUCT_resfinder")
C1 <- full_join(megahit_C_gene_argannot,megahit_C_gene_card,by="SEQUENCE")
C2 <- full_join(C1,megahit_C_gene_ncbi,by="SEQUENCE")
C3 <- full_join(C2,megahit_C_gene_resfinder,by="SEQUENCE")
#write.csv(C3,"combine_C.csv",row.names = F)

#setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_D")
setwd("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_D")
getwd()
megahit_D_gene_argannot <- read_excel("megahit_D_gene.argannot.xlsx")
colnames(megahit_D_gene_argannot) <- c("#FILE_argannot","SEQUENCE","START_argannot","END_argannot","GENE_argannot","COVERAGE_argannot","COVERAGE_MAP_argannot","GAPS_argannot","%COVERAGE_argannot","%IDENTITY_argannot","DATABASE_argannot","ACCESSION_argannot","PRODUCT_argannot")
megahit_D_gene_card <- read_excel("megahit_D_gene.card.xlsx")
colnames(megahit_D_gene_card) <- c("#FILE_card","SEQUENCE","START_card","END_card","GENE_card","COVERAGE_card","COVERAGE_MAP_card","GAPS_card","%COVERAGE_card","%IDENTITY_card","DATABASE_card","DNA Accession","PRODUCT_card")
megahit_D_gene_ncbi <- read_excel("megahit_D_gene.ncbi.xlsx")
colnames(megahit_D_gene_ncbi) <- c("#FILE_ncbi","SEQUENCE","START_ncbi","END_ncbi","GENE_ncbi","COVERAGE_ncbi","COVERAGE_MAP_ncbi","GAPS_ncbi","%COVERAGE_ncbi","%IDENTITY_ncbi","DATABASE_ncbi","ACCESSION_ncbi","PRODUCT_ncbi")
megahit_D_gene_resfinder <- read_excel("megahit_D_gene.resfinder.xlsx")
colnames(megahit_D_gene_resfinder) <- c("#FILE_resfinder","SEQUENCE","START_resfinder","END_resfinder","GENE_resfinder","COVERAGE_resfinder","COVERAGE_MAP_resfinder","GAPS_resfinder","%COVERAGE_resfinder","%IDENTITY_resfinder","DATABASE_resfinder","ACCESSION_resfinder","PRODUCT_resfinder")
D1 <- full_join(megahit_D_gene_argannot,megahit_D_gene_card,by="SEQUENCE")
D2 <- full_join(D1,megahit_D_gene_ncbi,by="SEQUENCE")
D3 <- full_join(D2,megahit_D_gene_resfinder,by="SEQUENCE")
#write.csv(D3,"combine_D.csv",row.names = F)


#setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C/kaiju_arg_gene_C")
setwd("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_C/kaiju_arg_gene_C")
kaiju_arg_genes_C_names_out <- read_excel("kaiju_arg_genes_C.names.out.xlsx")
arg_bac_C <- full_join(megahit_C_gene_card[,c(2,5)],kaiju_arg_genes_C_names_out ,by = "SEQUENCE") 

for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,10] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,10] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,4] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,4] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,5] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,5] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,6] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,6] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,7] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,7] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,8] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,8] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,9] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,9] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,10] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,10] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,11] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,11] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,12] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,12] <- 'unknown'}}
for (i in 1:nrow(arg_bac_C)){
  if (arg_bac_C[i,13] %in% c(NA, 'NA', '')) {
    arg_bac_C[i,13] <- 'unknown'}}


for (i in 1:nrow(arg_bac_C)) {
  if (!(grepl('group$',arg_bac_C[i,5]) | grepl('^un',arg_bac_C[i,5]))){
    for (k in 15:6){
      arg_bac_C[i,k] <- arg_bac_C[i,k-1]
    }
    arg_bac_C[i,5] <- 'unknown'
    }}
#table(arg_bac_C$group)
# clean phylum
for (i in 1:nrow(arg_bac_C)) {
  if (grepl('phyla$',arg_bac_C[i,6]) | grepl('^unclassified Bacteria',arg_bac_C[i,6])){
    for (k in 15:10){
    arg_bac_C[i,k] <- arg_bac_C[i,k-3] #15=12,10=7
  }
  arg_bac_C[i,6] <- 'unknown'}}

for (i in 1:nrow(arg_bac_C)) {
  if (grepl('group$',arg_bac_C[i,6])){
    for (k in 6:14){
      arg_bac_C[i,k] <- arg_bac_C[i,k+1] #6=7, 7=8
    }}}
#table(arg_bac_C$Phylum)

table(arg_bac_C$Class)

for (i in 1:nrow(arg_bac_C)) {
  if (grepl('delta/epsilon subdivisions',arg_bac_C[i,7])){
    for (k in 7:14){
      arg_bac_C[i,k] <- arg_bac_C[i,k+1] #6=7, 7=8
    }}}


# clean order
#table(arg_bac_C$Order)
#x <- subset(arg_bac_C,!grepl('ales$',arg_bac_C$Order))

for (i in 1:nrow(arg_bac_C)) {
  if (grepl('Sphaerobacteridae',arg_bac_C[i,8])){
    for (k in 8:14){
      arg_bac_C[i,k] <- arg_bac_C[i,k+1] #8=9
    }}}

pattern='Verrucomicrobiae bacterium DG1235|Thermorudis|'
for (i in 1:nrow(arg_bac_C)) {
  if (grepl(pattern,arg_bac_C[i,9])){
      arg_bac_C[i,10] <- arg_bac_C[i,9] #8=9
    }}

# clean family
#table(arg_bac_C$Family)
#x <- subset(arg_bac_C,!grepl('ceae$',arg_bac_C$Family)&!grepl('^un',arg_bac_C$Family))

for (i in 1:nrow(arg_bac_C)) {
  if (arg_bac_C[i,9] ==arg_bac_C[i,10]){
    arg_bac_C[i,10] <- arg_bac_C[i,11]
    arg_bac_C[i,11] <- 'unknown'}
  }

for (i in 1:nrow(arg_bac_C)) {
  if (grepl('^un',arg_bac_C[i,10]) & !grepl('^un',arg_bac_C[i,12])){
    arg_bac_C[i,10] <- arg_bac_C[i,12]}
}

#x <- subset(arg_bac_C,!grepl('^un',arg_bac_C$Family)&!grepl('ceae$',arg_bac_C$Family)&!grepl('neae$',arg_bac_C$Family)&grepl('^un',arg_bac_C$Genus))

for (i in 1:nrow(arg_bac_C)) {
  if (!grepl('^un',arg_bac_C[i,9])&!grepl('ceae$',arg_bac_C[i,9])&!grepl('neae$',arg_bac_C[i,9])&grepl('^un',arg_bac_C[i,10])){
    arg_bac_C[i,10] <- arg_bac_C[i,9]}
}



#kaiju_arg_genes_D_names_out <- read_excel("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_D/kaiju_arg_genes_D/kaiju_arg_genes_D.names.out.xlsx")
kaiju_arg_genes_D_names_out <- read_excel("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_D/kaiju_arg_genes_D/kaiju_arg_genes_D.names.out.xlsx")
arg_bac_D <- full_join(megahit_D_gene_card[,c(2,5)],kaiju_arg_genes_D_names_out,by = "SEQUENCE") 
#colnames(arg_bac_D) <- c('sequence','gene_card','cellular organisms','domain','group','Phylum','Class','Order','Family','Genus','Species','Strain','Strain1')
table(arg_bac_C$`cellular organisms`)
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,10] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,10] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,4] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,4] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,5] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,5] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,6] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,6] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,7] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,7] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,8] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,8] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,9] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,9] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,10] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,10] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,11] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,11] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,12] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,12] <- 'unknown'}}
for (i in 1:nrow(arg_bac_D)){
  if (arg_bac_D[i,13] %in% c(NA, 'NA', '')) {
    arg_bac_D[i,13] <- 'unknown'}}


for (i in 1:nrow(arg_bac_D)) {
  if (!(grepl('group$',arg_bac_D[i,5]) | grepl('^un',arg_bac_D[i,5]))){
    for (k in 15:6){
      arg_bac_D[i,k] <- arg_bac_D[i,k-1]
    }
    arg_bac_D[i,5] <- 'unknown'
  }}
#table(arg_bac_D$group)
# clean phylum
for (i in 1:nrow(arg_bac_D)) {
  if (grepl('phyla$',arg_bac_D[i,6]) | grepl('^unclassified Bacteria',arg_bac_D[i,6])){
    for (k in 15:10){
      arg_bac_D[i,k] <- arg_bac_D[i,k-3] #15=12,10=7
    }
    arg_bac_D[i,6] <- 'unknown'}}

for (i in 1:nrow(arg_bac_D)) {
  if (grepl('group$',arg_bac_D[i,6])){
    for (k in 6:14){
      arg_bac_D[i,k] <- arg_bac_D[i,k+1] #6=7, 7=8
    }}}
#table(arg_bac_D$Phylum)

table(arg_bac_D$Class)

for (i in 1:nrow(arg_bac_D)) {
  if (grepl('delta/epsilon subdivisions',arg_bac_D[i,7])){
    for (k in 7:14){
      arg_bac_D[i,k] <- arg_bac_D[i,k+1] #6=7, 7=8
    }}}


# clean order
#table(arg_bac_D$Order)
#x <- subset(arg_bac_D,!grepl('ales$',arg_bac_D$Order))

for (i in 1:nrow(arg_bac_D)) {
  if (grepl('Sphaerobacteridae',arg_bac_D[i,8])){
    for (k in 8:14){
      arg_bac_D[i,k] <- arg_bac_D[i,k+1] #8=9
    }}}

pattern='Verrucomicrobiae bacterium DG1235|Thermorudis|'
for (i in 1:nrow(arg_bac_D)) {
  if (grepl(pattern,arg_bac_D[i,9])){
    arg_bac_D[i,10] <- arg_bac_D[i,9] #8=9
  }}

# clean family
#table(arg_bac_D$Family)
#x <- subset(arg_bac_D,!grepl('ceae$',arg_bac_D$Family)&!grepl('^un',arg_bac_D$Family))

for (i in 1:nrow(arg_bac_D)) {
  if (arg_bac_D[i,9] ==arg_bac_D[i,10]){
    arg_bac_D[i,10] <- arg_bac_D[i,11]
    arg_bac_D[i,11] <- 'unknown'}
}

for (i in 1:nrow(arg_bac_D)) {
  if (grepl('^un',arg_bac_D[i,10]) & !grepl('^un',arg_bac_D[i,12])){
    arg_bac_D[i,10] <- arg_bac_D[i,12]}
}

#x <- subset(arg_bac_D,!grepl('^un',arg_bac_D$Family)&!grepl('ceae$',arg_bac_D$Family)&!grepl('neae$',arg_bac_D$Family)&grepl('^un',arg_bac_D$Genus))

for (i in 1:nrow(arg_bac_D)) {
  if (!grepl('^un',arg_bac_D[i,9])&!grepl('ceae$',arg_bac_D[i,9])&!grepl('neae$',arg_bac_D[i,9])&grepl('^un',arg_bac_D[i,10])){
    arg_bac_D[i,10] <- arg_bac_D[i,9]}
}

y <- as.data.frame(table(arg_bac_D$Genus))
z <- as.data.frame(table(arg_bac_D$GENE_card))
z1 <- as.data.frame(table(arg_bac_C$GENE_card))
z0 <- full_join(z,z1, by="Var1")

#setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/card-data/")
setwd("G:/My Drive/research/shogun_metagenome/card-data/")
aro_categories_index <- read_excel("aro_categories_index.xlsx")
colnames(aro_categories_index) <- c("Protein Accession","ACCESSION","AMR Gene Family","Drug Class","Resistance Mechanism" )

#setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C")
setwd("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_C")
megahit_C_gene_card <- read_excel("megahit_C_gene.card.xlsx")
megahit_C_gene_card_1 <- megahit_C_gene_card %>% dplyr::select(SEQUENCE,GENE,ACCESSION)
colnames(megahit_C_gene_card) <- c("SEQUENCE",'GENE','ACCESSION')
head(aro_categories_index)
library(splitstackshape)
megahit_C_gene_card_2 <- cSplit(megahit_C_gene_card_1,"ACCESSION",sep=":")
megahit_C_gene_card_2 <- megahit_C_gene_card_2[,-4]
head(megahit_C_gene_card_2)
colnames(megahit_C_gene_card_2) <- c("SEQUENCE",'GENE','ACCESSION')
accession_card_C <-  inner_join(megahit_C_gene_card_2,aro_categories_index,by="ACCESSION")
dim(megahit_C_gene_card_2)
dim(aro_categories_index)
aro_index <- read_excel("~/Downloads/card-data/aro_index.xlsx")


# DeepARG result
library(readxl)
library(dplyr)
deeparg_gene_C_ARG <- read_excel("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_C/deeparg/deeparg_gene_C.ARG.xlsx")
deeparg_gene_D_ARG <- read_excel("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_D/deeparg_out/deeparg_gene_D.ARG.xlsx")

#deeparg_gene_C_ARG <- read_excel("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C/deeparg/deeparg_gene_C.ARG.xlsx")
#deeparg_gene_D_ARG <- read_excel("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_D/deeparg_out/deeparg_gene_D.ARG.xlsx")
C_gene <- as.data.frame(deeparg_gene_C_ARG[,1]) 
D_gene <- as.data.frame(deeparg_gene_D_ARG[,1]) 
t_C_gene <- as.data.frame(table(C_gene))
colnames(t_C_gene) <- c("gene","C_count")
t_D_gene <- as.data.frame(table(D_gene))
head(t_D_gene)
head(t_C_gene)
colnames(t_D_gene) <- c("gene","D_count")
m_C_D_gene <- full_join(t_C_gene,t_D_gene,by='gene')
m_C_D_gene[is.na(m_C_D_gene)] <- 0
heat_C_D_gene2 <- m_C_D_gene %>% filter(C_count==0 | D_count==0)

C_class <- as.data.frame(deeparg_gene_C_ARG[,5]) 
D_class <- as.data.frame(deeparg_gene_D_ARG[,5]) 
t_C_class <- as.data.frame(table(C_class))
colnames(t_C_class) <- c("class","C_count")
t_D_class <- as.data.frame(table(D_class))
colnames(t_D_class) <- c("class","D_count")
m_C_D_class <- full_join(t_C_class,t_D_class,by='class')
m_C_D_class[is.na(m_C_D_class)] <- 0
m_C_D_gene <- m_C_D_gene %>% mutate(normC=C_count/133,normD=D_count/155) # 133 16S rRNA gene of C, 155 of D analyzed from barrnap
m_C_D_class <- m_C_D_class %>% mutate(normC=C_count/133,normD=D_count/155)
sum(as.numeric(m_C_D_class$C_count), na.rm = TRUE)# 738 genes were found from 694 contigs, 738 ORFs 
sum(as.numeric(m_C_D_class$D_count), na.rm = TRUE) # 946 genes were found from 890 contigs, 946 ORFs

heat_C_D_gene <- m_C_D_gene[,c(1,4:5)]
heat_C_D_class <- m_C_D_class[,c(1,4:5)]
colnames(heat_C_D_gene) <- c("gene","C","D")
colnames(heat_C_D_class) <- c("ARGtype","C","D")
heat_C_D_gene$gene[heat_C_D_gene$gene =='vanU'] <- 'VanUG'


heat_gene_diff <- heat_C_D_gene %>% mutate(diff= D-C,diff_p= (D-C)/C) # find out the gene difference before and after treatment
heat_gene_diff_fiter <- heat_gene_diff %>% filter(diff > 0.01 | diff< -0.01)
heat_gene_plot <- heat_gene_diff_fiter %>%  select(gene,diff) %>% filter(diff >=0.05) 
heat_gene_plot$gene <- factor(heat_gene_plot$gene, levels = heat_gene_plot$gene[order(heat_gene_plot$diff)])
ggplot(data=heat_gene_plot, aes(x=gene, y=diff)) + 
  geom_bar(stat="identity",width=0.6) +
  ylab("Normalized abundance(ARGs/16S rRNA gene copies)") +
  xlab("ARGs") +
  theme_bw() +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size=10),
        axis.title=element_text(size=10,face="bold",family='Times'),
        axis.text.x = element_text(hjust = 0,size=10,color = "black",family='Times'),
        axis.text.y = element_text(size=10,color = "black",family='Times'))


library(reshape)
library(ggplot2)
setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C/deeparg/")
setwd("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_C/deeparg/")
#write.csv(heat_C_D_gene,"heat_C_D_gene.csv",row.names = F)
heat_C_D_gene1 <- heat_C_D_gene[apply(heat_C_D_gene[,-1], MARGIN = 1, function(x) any(x >0.05)), ]
#heat_C_D_gene1 <- heat_C_D_gene[apply(heat_C_D_gene[,-1], MARGIN = 1, function(x) all(x <= 0.2)), ]

heat_gene <- melt(heat_C_D_gene1)
heat_class <- melt(heat_C_D_class)
# VanU and VanUg are the same thing, just from different database. Merge them

heat_gene$gene[39] = "VanUG"

heat_gene$gene[heat_gene$gene =='vanU'] <- 'VanUG'
heat_gene <- heat_gene %>% group_by(gene,variable) %>% summarise_each(funs(sum))

ggplot(heat_gene, aes(x=gene,y=value, fill=variable)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_discrete()+
  ylab("Normalized abundance") +
  xlab("ARG subtypes") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        axis.text.x = element_text(hjust = 1,size=10,color = "black"),
        axis.text.y = element_text(size=10,face='italic',color = "black")) +
  labs(fill = "Digestion")+coord_flip()

ggplot(heat_gene, aes(variable, gene)) +
  geom_tile(aes(fill = value), color = "black") +
  scale_fill_gradient(low = "#FFFFCC", high = "#00703C") +
  ylab("ARGs") +
  xlab("sample") +
  theme(legend.title = element_text(size = 12,family = "Times New Roman"),
        legend.text = element_text(size = 12,family = "Times"),
        plot.title = element_text(size=12,family = "Times"),
        axis.title=element_text(size=12,family = "Times"),
        axis.text.x = element_text(angle = 90, hjust = 1,size=10,face="italic",color = "black",family = "Times New Roman"),
        text=element_text(family = "Times",colour = "black")) +
  labs(fill = "Normalized Abundance")+coord_flip()


heat_class[heat_class=="macrolide-lincosamide-streptogramin"] <- "MLS" #MLS is short for macrolide-lincosamide-streptogramin

ggplot(heat_class, aes(variable, ARGtype )) +
  geom_tile(aes(fill = value), color = "black") +
  scale_fill_gradient(low = "#FDDB6D", high = "#00703C") +
  ylab("ARGTypes ") +
  xlab("Digester") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1,size=12,face="italic"),
        text=element_text(family = "Times")) +
  labs(fill = "Normalized Abundance")+coord_flip()


#======================================= Venn Diagram ==========================================
library(systemPipeR)
library(VennDiagram)
heat_C_D_gene <- m_C_D_gene[,c(1,4:5)]
for (i in 1:nrow(heat_C_D_gene)){
  for (j in 2:3){
    if (heat_C_D_gene[i,j] == 0){
      heat_C_D_gene[i,j] = paste('FALSE')}
    else (heat_C_D_gene[i,j] <-  paste(heat_C_D_gene[i,1]))}
}

head(heat_C_D_gene)
colnames(heat_C_D_gene) <- c("gene","C","D")

dat <-heat_C_D_gene[,-1] # delete the the first row, since your first row is row names
da <- as.list(dat)
for (i in (1:length(da))){
  x <- da[[i]] 
  da[[i]]=x[x!="FALSE"]
} # I delete all the 'FLASE', which means 0 in the list

ven <- overLapper(da[1:2], type="vennsets") # here da[1:3] means the first 3 columns in your raw dataframe, you can change based on your own settings
vennPlot(ven)


#============================= Mapping Nanopore =======================================

mapping_C_nanopore <- read_excel("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_C/mapping_C_nanopore.xlsx")
head(mapping_C_nanopore)
C_MGE <- mapping_C_nanopore %>% filter(Group=="MGEs")
C_MGE1 <- C_MGE %>% filter(`E-value`< 1E-10)
C_t_MGE <- as.data.frame(table(C_MGE1$Category))  
C_MGE2 <- C_MGE %>% filter(`E-value`< 1E-10 & identity >25 & coverage>50)
C_t_MGE1 <- as.data.frame(table(C_MGE2$Category))
head(C_t_MGE1)

mapping_nanopore_ARG_D <- read_excel("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_D/mapping_nanopore_ARG_D.xlsx")
head(mapping_nanopore_ARG_D)
D_MGE <- mapping_nanopore_ARG_D %>% filter(Group=="MGEs")
head(D_MGE)
D_MGE1 <- D_MGE %>% filter(`E-value`< 1E-10)
D_t_MGE <- as.data.frame(table(D_MGE1$Category))  
D_MGE2 <- D_MGE %>% filter(`E-value`< 1E-10 & identity >25 & coverage>50)
D_t_MGE1 <- as.data.frame(table(D_MGE2$Category))
head(D_t_MGE1)
summary(D_MGE2)

D_C_MGE <- full_join(D_t_MGE1,C_t_MGE1,by="Var1")
D_C_MGE[is.na(D_C_MGE)] <- 0
colnames(D_C_MGE) <- c("MGE","D","C")
m_D_C_MGE <- melt(D_C_MGE,id="MGE")
m_D_C_MGE <- m_D_C_MGE %>% filter(value>2)
ggplot(m_D_C_MGE, aes(x=MGE,y=log(value))) +
  geom_bar(stat="identity",color="#696969",position=position_dodge(width = 0.8),aes(fill = variable),width = 0.6) +
  scale_fill_brewer(palette="Paired")+
  ylab("log(#genes)") +
  xlab("") +
  theme_bw()+
  theme(legend.title = element_text(size = 10,family = "Times",colour = "black"),
        legend.text = element_text(size = 10,family = "Times",colour = "black"),
        plot.title = element_text(size=16),
        axis.title=element_text(size=10,family = "Times",colour = "black"),
        axis.text.x = element_text(size=10,family = "Times",colour = "black"),
        axis.text.y = element_text(size=10,family = "Times",colour = "black")) +
  labs(fill = "Digester")

#============================== Circos plot ===================================
library(circlize)
m_C_D_class <- full_join(t_C_class,t_D_class,by='class')
m_C_D_class[is.na(m_C_D_class)] <- 0
m_C_D_class[m_C_D_class=="macrolide-lincosamide-streptogramin"] <- "MLS"

#temp = m_C_D_class
#idx = temp$C_count<10 & temp$D_count<10
#temp[idx, 'class'] <- 'others'
#x = aggregate(.~class, temp, sum)
#m_C_D_class<- x

rownames(m_C_D_class) <- m_C_D_class$class
m_C_D_class <- m_C_D_class[,-1]
colnames(m_C_D_class) <- c("C","D")
cis <- data.matrix(m_C_D_class)
grid.col = c(C = "grey", D = "grey",
             aminocoumarin = "red", aminoglycoside = "blue", bacitracin = "green", beta_lactam = "purple", 
             bleomycin ='pink', chloramphenicol= "yellow", fosfomycin = "#B0BF1A",fosmidomycin='#7CB9E8',
             glycopeptide="#FFBF00",kasugamycin="#804040",MLS="#993300",multidrug="#34B334",
             mupirocin="#FF8B00",peptide="#FF9899",polymyxin="#551B8C",quinolone="#CD9575",
             rifampin="#D0FF14",sulfonamide="#E9D66B",tetracycline="#00DEA4",trimethoprim="#21ABCD",
             acriflavin="#FF3988",streptothricin="#BF4F51",tunicamycin="#F4BBFF")

#grid.col = c(C = "grey", D = "grey",
 #            aminocoumarin = "red", aminoglycoside = "blue", bacitracin = "green",others="yellow",
  #           fosmidomycin='#7CB9E8',glycopeptide="#FFBF00",MLS="#993300",multidrug="#34B334",
   #          mupirocin="#FF8B00",polymyxin="#551B8C",rifampin="#D0FF14",tetracycline="#00DEA4",
    #         trimethoprim="#21ABCD")
circos.par(gap.after = c(rep(3, nrow(cis)-1), 10, rep(2, ncol(cis)-1), 15))
chordDiagram(cis, grid.col = grid.col, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = uh(6, "mm")))

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 2)
}

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  xplot = get.cell.meta.data("xplot")
  
  circos.text(mean(xlim+1),ylim[1] +1, cex = 0.8, sector.name, niceFacing = T, adj = c(0, 0.1),facing = "clockwise")
}, bg.border = NA)


circos.bar(xlim, c(mean(ylim), mean(ylim)), lty = 3) # dotted line
by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.2, 0.5)
for(p in seq(by, 1, by = by)) {
  circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.2, 
              paste0(p*100, "%"), cex = 0.3, adj = c(0.5, 0), niceFacing = TRUE,facing = "clockwise")}



#===================== DeepARG run bacteria=======================
arg_bac_C
library(splitstackshape)

#kaiju_deeparg_C_names_out <- read_excel("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C/kaiju_deeparg_C.names.out.xlsx")

kaiju_deeparg_C_names_out <- read_excel("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_C/kaiju_deeparg_C.names.out.xlsx")
kaiju_deeparg_C_names_out <- cSplit(kaiju_deeparg_C_names_out,"tax",";")
colnames(kaiju_deeparg_C_names_out) <- c('contig','cellular organisms','domain','group','Phylum','Class','Order','Family','Genus','Species','Strain','Strain1')

deeparg_gene_C_ARG1 <- deeparg_gene_C_ARG [,c(1,4:5)]
colnames(deeparg_gene_C_ARG1) <- c("ARG",'contig','ARGtype')
deep <- full_join(kaiju_deeparg_C_names_out,deeparg_gene_C_ARG1,by="contig")
deep <- data.frame(lapply(deep, as.character), stringsAsFactors=FALSE)
str(deep)

for (i in 1:nrow(deep)){
  if (deep[i,4] %in% c(NA, 'NA', '')) {
    deep[i,4] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,5] %in% c(NA, 'NA', '')) {
    deep[i,5] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,6] %in% c(NA, 'NA', '')) {
    deep[i,6] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,7] %in% c(NA, 'NA', '')) {
    deep[i,7] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,8] %in% c(NA, 'NA', '')) {
    deep[i,8] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,9] %in% c(NA, 'NA', '')) {
    deep[i,9] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,10] %in% c(NA, 'NA', '')) {
    deep[i,10] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,11] %in% c(NA, 'NA', '')) {
    deep[i,11] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,12] %in% c(NA, 'NA', '')) {
    deep[i,12] <- 'unknown'}}


for (i in 1:nrow(deep)) {
  if (!(grepl('group$',deep[i,4]))){
    for (k in 12:5){
      deep[i,k] <- deep[i,k-1]
    }
    deep[i,4] <- 'unknown'
  }}
#table(deep$Phylum)
# clean phylum
for (i in 1:nrow(deep)) {
  if (grepl('group$',deep[i,5])){
    for (k in 5:11){
      deep[i,k] <- deep[i,k+1] #6=7, 7=8
    }}}

x <- subset(deep,(grepl('unclassified Bacteria',deep$Phylum)))


for (i in 1:nrow(deep)) {
  if (grepl('phyla$',deep[i,5]) | grepl('^unclassified Bacteria',deep[i,5])){
    for (k in 12:9){
      deep[i,k] <- deep[i,k-3] #15=12,10=7
    }
    deep[i,5] <- 'unknown'}}

#table(deep$Phylum)

table(deep$Class)

for (i in 1:nrow(deep)) {
  if (grepl('delta/epsilon subdivisions',deep[i,6])){
    for (k in 6:11){
      deep[i,k] <- deep[i,k+1] #6=7, 7=8
    }}}

#x <- subset(deep,grepl('Marinimicrobia bacterium SCGC AAA298-D23', deep$Class))

for (i in 1:nrow(deep)){
  if (grepl('Marinimicrobia bacterium SCGC AAA298-D23', deep[i,6])){
    deep[i,10] <- deep[i,6]
    deep[i,6] <- 'unknown'
  }
}

for (i in 1:nrow(deep)){
  if (grepl('ales$', deep[i,6])){
    deep[i,7] <- deep[i,6]
    deep[i,6] <- 'unknown'
  }
}

for (i in 1:nrow(deep)){
  if (grepl('Candidatus Omnitrophus',deep[i,6])){
    for (k in 9:11){ #9=6,10=7,11=8,12=9
      deep[i,k] <- deep[i,k-3]}}}
# clean order
#table(deep$Order)
x <- subset(deep,!grepl('ales$',deep$Order)&!grepl('^un',deep$Order) & !grepl('^Ca',deep$Order))

for (i in 1:nrow(deep)) {
  if (!grepl('ales$',deep[i,7])&!grepl('^un',deep[i,7]) & !grepl('^Ca',deep[i,7])){
    for (k in 11:9){ #9=7,10=8,11=9
      deep[i,k] <- deep[i,k-2]
    }
    deep[i,7] = 'unknown'
    deep[i,8] = 'unknown'
    }}

#x <- subset(deep,grepl('^un',deep$Genus))

for (i in 1:nrow(deep)){
  if (grepl('neae',deep[i,8]) & grepl('ceae',deep[i,9])){
    deep[i,9] <- deep[i,8]
    deep[i,8] <- 'unknown'
  }
}

for (i in 1:nrow(deep)){
  if (grepl('^un',deep[i,9]) & !grepl('^un',deep[i,10])){
    deep[i,9] <- deep[i,10]
  }
}

for (i in 1:nrow(deep)){
  if (grepl('^un',deep[i,9]) & !grepl('^un',deep[i,11])){
    deep[i,9] <- deep[i,11]
  }
}

for (i in 1:nrow(deep)){
  if (grepl('^un',deep[i,9]) & !grepl('^un',deep[i,8])& !grepl('^ceae',deep[i,8])){
    deep[i,9] <- deep[i,8]
  }
}

deep_C <- deep
#-----------------END C -------------------
# ---------------start D ----------------

#kaiju_deeparg_D_names_out <- read_excel("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_D/kaiju_deeparg_D.names.out.xlsx")
kaiju_deeparg_D_names_out <- read_excel("G:/My Drive/research/shogun_metagenome/annotation_results_megahit_D/kaiju_deeparg_D.names.out.xlsx")
kaiju_deeparg_D_names_out <- cSplit(kaiju_deeparg_D_names_out,"tax",";")
colnames(kaiju_deeparg_D_names_out) <- c('contig','cellular organisms','domain','group','Phylum','Class','Order','Family','Genus','Species','Strain','Strain1')

deeparg_gene_D_ARG1 <- deeparg_gene_D_ARG [,c(1,4:5)]
colnames(deeparg_gene_D_ARG1) <- c("ARG",'contig','ARGtype')
deep <- full_join(kaiju_deeparg_D_names_out,deeparg_gene_D_ARG1,by="contig")
deep <- data.frame(lapply(deep, as.character), stringsAsFactors=FALSE)
str(deep)

for (i in 1:nrow(deep)){
  if (deep[i,4] %in% c(NA, 'NA', '')) {
    deep[i,4] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,5] %in% c(NA, 'NA', '')) {
    deep[i,5] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,6] %in% c(NA, 'NA', '')) {
    deep[i,6] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,7] %in% c(NA, 'NA', '')) {
    deep[i,7] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,8] %in% c(NA, 'NA', '')) {
    deep[i,8] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,9] %in% c(NA, 'NA', '')) {
    deep[i,9] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,10] %in% c(NA, 'NA', '')) {
    deep[i,10] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,11] %in% c(NA, 'NA', '')) {
    deep[i,11] <- 'unknown'}}
for (i in 1:nrow(deep)){
  if (deep[i,12] %in% c(NA, 'NA', '')) {
    deep[i,12] <- 'unknown'}}


for (i in 1:nrow(deep)) {
  if (!(grepl('group$',deep[i,4]))){
    for (k in 12:5){
      deep[i,k] <- deep[i,k-1]
    }
    deep[i,4] <- 'unknown'
  }}
#table(deep$Phylum)
# clean phylum
for (i in 1:nrow(deep)) {
  if (grepl('group$',deep[i,5])){
    for (k in 5:11){
      deep[i,k] <- deep[i,k+1] #6=7, 7=8
    }}}

for (i in 1:nrow(deep)) {
  if (grepl('phyla$',deep[i,5]) | grepl('^unclassified Bacteria',deep[i,5])){
    for (k in 12:9){
      deep[i,k] <- deep[i,k-3] #15=12,10=7
    }
    deep[i,5] <- 'unknown'}}

#table(deep$Phylum)

table(deep$Class)

for (i in 1:nrow(deep)) {
  if (grepl('delta/epsilon subdivisions',deep[i,6])){
    for (k in 6:11){
      deep[i,k] <- deep[i,k+1] #6=7, 7=8
    }}}

#x <- subset(deep,grepl('Marinimicrobia bacterium SCGC AAA298-D23', deep$Class))
table(deep$Class)
for (i in 1:nrow(deep)){
  if (grepl('Marinimicrobia bacterium SCGC AAA298-D23', deep[i,6])){
    deep[i,10] <- deep[i,6]
    deep[i,6] <- 'unknown'
  }
}

for (i in 1:nrow(deep)){
  if (grepl('ales$', deep[i,6])){
    deep[i,7] <- deep[i,6]
    deep[i,6] <- 'unknown'
  }
}


for (i in 1:nrow(deep)){
  if (grepl('Candidatus Omnitrophus',deep[i,6])){
    for (k in 9:11){ #9=6,10=7,11=8,12=9
      deep[i,k] <- deep[i,k-3]}}}
# clean order
#table(deep$Order)
x <- subset(deep,!grepl('ales$',deep$Order)&!grepl('^un',deep$Order) & !grepl('^Ca',deep$Order))

for (i in 1:nrow(deep)) {
  if (!grepl('ales$',deep[i,7])&!grepl('^un',deep[i,7]) & !grepl('^Ca',deep[i,7])){
    for (k in 11:9){ #9=7,10=8,11=9
      deep[i,k] <- deep[i,k-2]
    }
    deep[i,7] = 'unknown'
    deep[i,8] = 'unknown'
  }}

#x <- subset(deep,grepl('^un',deep$Genus))

for (i in 1:nrow(deep)){
  if (grepl('neae',deep[i,8]) & grepl('ceae',deep[i,9])){
    deep[i,9] <- deep[i,8]
    deep[i,8] <- 'unknown'
  }
}

for (i in 1:nrow(deep)){
  if (grepl('^un',deep[i,9]) & !grepl('^un',deep[i,10])){
    deep[i,9] <- deep[i,10]
  }
}

for (i in 1:nrow(deep)){
  if (grepl('^un',deep[i,9]) & !grepl('^un',deep[i,11])){
    deep[i,9] <- deep[i,11]
  }
}

for (i in 1:nrow(deep)){
  if (grepl('^un',deep[i,9]) & !grepl('^un',deep[i,8])& !grepl('^ceae',deep[i,8])){
    deep[i,9] <- deep[i,8]
  }
}

x <- subset(deep,(grepl('unclassified Bacteria',deep$Phylum)))
table(deep$Phylum)
table(deep$Class)

deep_D <- deep

# merge C and D in phylum level for making the plot
phylum_deep_D <- deep_D %>% select(Phylum)
phylum_deep_D1 <- as.data.frame(table(phylum_deep_D))
colnames(phylum_deep_D1) <- c("phylum","D")
phylum_deep_C <- deep_C %>% select(Phylum)
phylum_deep_C1 <- as.data.frame(table(phylum_deep_C))
colnames(phylum_deep_C1) <- c("phylum","C")
phylum_merge <- inner_join(phylum_deep_C1,phylum_deep_D1,by='phylum')
#write.csv(phylum_merge,"phylum_merge.csv",row.names = F)
#write.csv(deep_D,"deep_D.csv",row.names = F)

# merge C and D in family level for making the plot
family_deep_D <- deep_D %>% select(Family)
family_deep_D1 <- as.data.frame(table(family_deep_D))
colnames(family_deep_D1) <- c("family","D")
family_deep_C <- deep_C %>% select(Family)
family_deep_C1 <- as.data.frame(table(family_deep_C))
colnames(family_deep_C1) <- c("family","C")
family_merge <- inner_join(family_deep_C1,family_deep_D1,by='family')
#write.csv(family_merge,"family_merge.csv",row.names = F)
#--------------- D END----------------------------------

#============================network bacteria and ARG======
library(igraph)
nC <- deep_C[,c(9,13:14)]
nC$ARG[nC$ARG == 'vanU'] <- 'VanUG'
#nC <- nC %>% group_by(Genus,ARG,ARGtype) %>% summarise_each(funs(sum))
head(nC);dim(nC)
nC$count <- 1
nC2 <- nC[,c(1,4)]
head(nC2)
nC2$ARGtype <- "genus"
head(nC2)
nC2 <- nC2 %>% select(Genus,ARGtype,count)
head(nC2)

nC3 <- nC[,c(2:4)]
head(nC3)
colnames(nC2) <- c('genus','type','readsmapped')
colnames(nC3) <- c('genus','type','readsmapped')
nC23 <- rbind(nC2,nC3)
head(nC23)
nC23_nodes<- aggregate(readsmapped~ genus+type,nC23,sum) # nodes
head(nC23_nodes);dim(nC23_nodes)

link_C <- deep_C[,c(9,13)] 

link_C$ARG[link_C$ARG == 'vanU'] <- 'VanUG'
link_C$count <- 1
head(link_C)
link_C1 <- aggregate(count~.,link_C,sum)
head(link_C1)
colnames(link_C1) <- c("from","to","weight")

unique(nC23_nodes$type)
nC23_nodes$col <- as.numeric(as.factor(nC23_nodes$type))
nC23_nodes <- nC23_nodes %>% mutate(logreads=(readsmapped/133)^0.5)
colnames(nC23_nodes) <- c("genus","phylumfamily","readsmapped","type","logreads")

#setwd("/Users/luyang/Downloads/annotation_results_megahit_C")
#write.csv(nC23_nodes,"nC23_nodes.csv",row.names = F)
#write.csv(link_C1,"link_C1.csv",row.names = F)
colrs <- c("#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F",
           "#FF3F34","#76D7EA","#8359A3","grey","#FBE870",
           "#C5E17A","#03BB85","#FFDF00","#0A6B0D","#8FD8D8",
           "#A36F40","#F653A6","#CA3435","#FFCBA4","#FF99CC",
           "#FA9D5A")

nrow(nC23_nodes); length(unique(nC23_nodes$genus))
nrow(link_C1); nrow(unique(link_C1[,c("from", "to")]))
link_C1 <- aggregate(link_C1[,3], link_C1[,-3], sum)
link_C1 <- link_C1[order(link_C1$from, link_C1$to),]
#colnames(link_C1)[4] <- "weight"
rownames(link_C1) <- NULL
head(link_C1);dim(link_C1)

net <- graph_from_data_frame(d=link_C1, vertices=nC23_nodes, directed=F) 
#net <- delete_edges(net, E(net)[x<2])
#plot(net, edge.arrow.size=.4,vertex.label=NA)
# nodes color = grey(virus) and colorful(bacteria), colorful=c(different color is different phylum)
unique(nC23_nodes$phylumfamily)
#
E(net)
V(net)$color <- colrs[V(net)$type]
V(net)$size <- V(net)$logreads*15
V(net)$label.color <- "grey20"
V(net)$label <- NA
E(net)$width <- E(net)$x
E(net)$arrow.size <- .3
E(net)$edge.color <- "white"
E(net)$width <- 1+E(net)$x/7
plot(net, vertex.label.font=1,vertex.label.cex=.7)

#vertex.label=nC23_nodes$genus,
legend(x=0.95, y=0, c("aminocoumarin","aminoglycoside","bacitracin","beta_lactam","bleomycin","chloramphenicol",
                        "fosfomycin","fosmidomycin","genus","glycopeptide","kasugamycin","macrolide-lincosamide-streptogramin",
                        "multidrug","mupirocin","peptide","polymyxin","quinolone","rifampin",
                        "sulfonamide","tetracycline","trimethoprim"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=1.5, cex=0.9, bty="n", ncol=1,y.intersp=0.3)
par(family="Times")
tkid <- tkplot(net) 


#----- code for filtered network-------
link_Ca <- link_C1[which(link_C1$x>1),]
colnames(link_Ca) <- c("from","to","weight")
from <- as.data.frame(link_Ca[,1]) 
colnames(from) <- "genus"
#from <- from %>% filter(genus != "unknown")
to <- as.data.frame(link_Ca[,2]) 
colnames(to) <- "genus"
genus <- rbind(from,to)
genus <- unique(genus)



nC23_nodes$name <- NA
nC23_nodes[nC23_nodes$phylumfamily!="genus","name"] <- nC23_nodes[nC23_nodes$phylumfamily!="genus","genus"]
nC23_nodes[nC23_nodes$phylumfamily=="genus","name"] <- ""
#nC23_nodes[nC23_nodes$phylumfamily=="genus","name"] <- as.numeric(as.factor(nC23_nodes[nC23_nodes$phylumfamily=="genus","genus"]))

nC23_nodes$shape <- NA
nC23_nodes[nC23_nodes$phylumfamily!="genus","shape"] <- "circle"
nC23_nodes[nC23_nodes$phylumfamily=="genus","shape"] <- "square"

n=1
for (i in 1:nrow(nC23_nodes)){
  if (nC23_nodes[i,2]=='genus'){
    nC23_nodes[i,6]=n
    n=n+1
  }
}
#write.csv(nC23_nodes,"nC23_nodes_label_genus_name.csv",row.names = F)

nodes_C <- nC23_nodes
nodes_Ca <- inner_join(nodes_C,genus)
head(nodes_Ca)
dim(nodes_Ca)
nodes_Ca <- nodes_Ca[,-4]
nodes_Ca$type <- as.numeric(as.factor(nodes_Ca$phylumfamily))

unique(nodes_Ca$phylumfamily)
colrs <- c("#76D7EA","#8359A3","white","#AF593E","#C5E17A",
           "#03BB85","#FFFDD0","#0A6B0D","#FA9D5A","#FFFF33","#F653A6","#CA3435","#FFCBA4")
nrow(nodes_Ca); length(unique(nodes_Ca$genus))
nrow(link_Ca); nrow(unique(link_Ca[,c("from", "to")]))
link_Ca <- aggregate(link_Ca[,3], link_Ca[,-3], sum)
link_Ca <- link_Ca[order(link_Ca$from, link_Ca$to),]
#colnames(link_Ca)[4] <- "weight"
rownames(link_Ca) <- NULL
net <- graph_from_data_frame(d=link_Ca, vertices=nodes_Ca, directed=F) 
#plot(net, edge.arrow.size=.4,vertex.label=NA)
# nodes color = grey(virus) and colorful(bacteria), colorful=c(different color is different phylum)

#
E(net)
V(net)$color <- colrs[V(net)$type]
V(net)$size <- V(net)$logreads*20
V(net)$label.color <- "black"
V(net)$label <- NA
E(net)$width <- E(net)$x
E(net)$arrow.size <- .3
E(net)$edge.color <- "red"
E(net)$width <- 1+E(net)$x/6
#V(net)$shape <- ifelse(V(net)$shape, "circle", "square")

#l <- layout_with_fr(net)
#l <- norm_coords(l, ymin=-3, ymax=3, xmin=-3, xmax=3)
#plot(net, vertex.label.font=1,vertex.label.cex=.7)
#plot(net, vertex.label.font=1,edge.color="304",rescale=F,vertex.label.cex=.7,vertex.frame.color="gray",vertex.label.dist=0.6,vertex.label=nodes_Ca$name,layout=l*0.45)
#vertex.label.dist=1,
plot(net, vertex.label.font=1,vertex.frame.color="gray",vertex.label.cex=.8,vertex.label=nodes_Ca$name)
#plot(net, vertex.label.font=1,vertex.label.cex=.7,vertex.label=nodes_Ca$name,layout=layout_on_sphere)
#plot(net, vertex.label.font=1,vertex.label.cex=.7,vertex.label=nodes_Ca$name,layout=layout_with_kk)

legend(x=1.1, y=0.3, c("aminoglycoside","fosfomycin","genus","glycopeptide","kasugamycin","macrolide-lincosamide-streptogramin",
                        "multidrug","mupirocin","peptide","polymyxin","rifampin","tetracycline","trimethoprim"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=1.2, cex=1, bty="n", ncol=1,y.intersp=0.45)
par(family="Times")
tkid <- tkplot(net, vertex.label.font=1,vertex.label.cex=.6,vertex.label=nodes_Ca$name,layout=layout_with_fr) 
setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_D/")
#write.csv(nodes_Ca,'name_nodes.csv',row.names = F)
#-------------D digester----------------
nD <- deep_D[,c(9,13:14)]
head(nD)
nD$ARG[nD$ARG == 'vanU'] <- 'VanUG'
nD <- nD %>% group_by(Genus,ARG,ARGtype) %>% summarise_each(funs(sum))
nD$count <- 1
nD2 <- nD[,c(1,4)]
nD2$ARGtype <- "genus"
head(nD2)
nD2 <- nD2 %>% select(Genus,ARGtype,count)
nD3 <- nD[,c(2:4)]
head(nD2)
head(nD3)
colnames(nD2) <- c('genus','type','readsmapped')
colnames(nD3) <- c('genus','type','readsmapped')
nD23 <- rbind(nD2,nD3)
head(nD23)
nD23_nodes<- aggregate(readsmapped~ genus+type,nD23,sum) # nodes


link_D <- deep_D[,c(9,13)] 
link_D$ARG[link_D$ARG == 'vanU'] <- 'VanUG'
link_D$count <- 1
head(link_D)
link_D1 <- aggregate(count~.,link_D,sum)


head(link_D1)
colnames(link_D1) <- c("from","to","weight")

unique(nD23_nodes$type)
nD23_nodes$col <- as.numeric(as.factor(nD23_nodes$type))
nD23_nodes <- nD23_nodes %>% mutate(logreads=(readsmapped/155)^0.5)
colnames(nD23_nodes) <- c("genus","phylumfamily","readsmapped","type","logreads")

setwd("/Users/luyang/Downloads/annotation_results_megahit_D")
#write.csv(nD23_nodes,"nD23_nodes.csv",row.names = F)
#write.csv(link_D1,"link_D1.csv",row.names = F)
unique(nodes_Da$phylumfamily)
colrs <- c("#0066FF","#AF593E","#01A368","#FF861F","#ED0A3F",
           "#FF3F34","#76D7EA","#8359A3","grey","#FBE870",
           "#C5E17A","#03BB85","#FFDF00","#0A6B0D","#8FD8D8",
           "#A36F40","#F653A6","#CA3435","#FFCBA4","#FF99CC",
           "#FA9D5A")

nrow(nD23_nodes); length(unique(nD23_nodes$genus))
nrow(link_D1); nrow(unique(link_D1[,c("from", "to")]))
link_D1 <- aggregate(link_D1[,3], link_D1[,-3], sum)
link_D1 <- link_D1[order(link_D1$from, link_D1$to),]
#colnames(link_D1)[4] <- "weight"
rownames(link_D1) <- NULL
net <- graph_from_data_frame(d=link_D1, vertices=nD23_nodes, directed=F) 
#net <- delete_edges(net, E(net)[x<2])
#plot(net, edge.arrow.size=.4,vertex.label=NA)
# nodes color = grey(virus) and colorful(bacteria), colorful=c(different color is different phylum)
unique(nD23_nodes$phylumfamily)
#
E(net)
V(net)$color <- colrs[V(net)$type]
V(net)$size <- V(net)$logreads*15
V(net)$label.color <- "black"
V(net)$label <- NA
E(net)$width <- E(net)$x
E(net)$arrow.size <- .3
E(net)$edge.color <- "white"
E(net)$width <- 1+E(net)$x/7
plot(net, vertex.label.font=1,vertex.label.cex=.7)


#vertex.label=nC23_nodes$genus,
legend(x=0.95, y=0, c("acriflavin","aminocoumarin","aminoglycoside",
                      "bacitracin","beta_lactam","bleomycin",
                      "chloramphenicol","fosfomycin","fosmidomycin",
                      "genus","glycopeptide","kasugamycin",
                      "macrolide-lincosamide-streptogramin","multidrug","mupirocin",
                      "peptide","polymyxin","quinolone",
                      "rifampin","streptothricin","sulfonamide",
                      "tetracycline","trimethoprim","tunicamycin"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=1.5, cex=0.9, bty="n", ncol=1,y.intersp=0.3)
par(family="Times")
#----- code for filtered network-------
link_Da <- link_D1[which(link_D1$x>1),] # pay attention to link_D1$
colnames(link_Da) <- c("from","to","weight")
head(link_Da)
from <- as.data.frame(link_Da[,1]) 
colnames(from) <- "genus"
#from <- from %>% filter(genus != "unknown")
to <- as.data.frame(link_Da[,2]) 
colnames(to) <- "genus"
genus <- rbind(from,to)
genus <- unique(genus)
head(genus)

nD23_nodes$name <- NA
nD23_nodes[nD23_nodes$phylumfamily!="genus","name"] <- nD23_nodes[nD23_nodes$phylumfamily!="genus","genus"]
nD23_nodes[nD23_nodes$phylumfamily=="genus","name"] <- ""
nD23_nodes[nD23_nodes$phylumfamily=="genus","name"] <- nD23_nodes[nD23_nodes$phylumfamily!="genus","name"]
#nD23_nodes[nD23_nodes$phylumfamily=="genus","name"] <- as.numeric(as.factor(nD23_nodes[nD23_nodes$phylumfamily=="genus","genus"]))
head(nD23_nodes)

nD23_nodes$shape <- NA
nD23_nodes[nD23_nodes$phylumfamily!="genus","shape"] <- "circle"
nD23_nodes[nD23_nodes$phylumfamily=="genus","shape"] <- "square"
n=1
for (i in 1:nrow(nD23_nodes)){
  if (nD23_nodes[i,2]=='genus'){
    nD23_nodes[i,6]=n
    n=n+1
  }
}
#write.csv(nD23_nodes,"nD23_nodes_label_genus_name.csv",row.names = F)

nodes_D <- nD23_nodes
nodes_Da <- inner_join(nodes_D,genus)
head(nodes_Da)
dim(nodes_Da)
nodes_Da <- nodes_Da[,-4]
nodes_Da$type <- as.numeric(as.factor(nodes_Da$phylumfamily))

unique(nodes_Da$phylumfamily)
unique(nodes_Ca$phylumfamily)

colrs <- c("grey","#76D7EA","pink","#8359A3","green",
           "white","#AF593E","#C5E17A","#03BB85","#FFFDD0",
           "#0A6B0D","#FA9D5A",	"#FFFF33","tomato","#F653A6"
           ,"#CA3435","#FFCBA4")

#,"#FFCBA4","#FF99CC","#FA9D5A","#CA3435"
#c:::colrs <- c("#76D7EA","#8359A3","white","#AF593E","#C5E17A",
 #          "#03BB85","#FFFDD0","#0A6B0D","#FA9D5A",	"#FFFF33",
  #         "#F653A6","#CA3435","#FFCBA4")

nrow(nodes_Da); length(unique(nodes_Da$genus))
name_genus_D <- nodes_Da %>% filter(nodes_Da$phylumfamily=='genus') %>% select(genus,name)
head(nodes_Da)
colnames(name_genus_D) <- c('genus','name_D')
name_genus_C <- nodes_Ca %>% filter(nodes_Ca$phylumfamily=='genus') %>% select(genus,name)
colnames(name_genus_C) <- c('genus','name_C')
name_D_C <- full_join(name_genus_D,name_genus_C,by='genus')
head(name_D_C)
getwd()
setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C/kaiju_arg_gene_C")
#write.csv(name_D_C,'name_D_C.csv',row.names = F) # This csv will cover C and D genus name label in the networks
#write.csv(nodes_Da, 'nodes_Da.csv',row.names = F) # this csv will cover the the number lable explanation of genus level name
nrow(link_Da); nrow(unique(link_Da[,c("from", "to")]))
link_Da <- aggregate(link_Da[,3], link_Da[,-3], sum)
link_Da <- link_Da[order(link_Da$from, link_Da$to),]
#colnames(link_Da)[4] <- "weight"
rownames(link_Da) <- NULL
net <- graph_from_data_frame(d=link_Da, vertices=nodes_Da, directed=F) 
#plot(net, edge.arrow.size=.4,vertex.label=NA)
# nodes color = grey(virus) and colorful(bacteria), colorful=c(different color is different phylum)

#
E(net)
V(net)$color <- colrs[V(net)$type]
V(net)$size <- V(net)$logreads*20
V(net)$label.color <- "black"
V(net)$label <- NA
E(net)$width <- E(net)$x
E(net)$arrow.size <- .3
E(net)$edge.color <- "red"
E(net)$width <- 1+E(net)$x/6
#V(net)$shape <- ifelse(V(net)$shape, "circle", "square")

#l <- layout_with_fr(net)
#l <- norm_coords(l, ymin=-3, ymax=3, xmin=-3, xmax=3)
plot(net, vertex.label.font=1,vertex.label.cex=.7,vertex.frame.color="gray",vertex.label=nodes_Da$name) # PLS take care of the C and D lables

#vertex.label.dist=1,
#plot(net, vertex.label.font=1,vertex.label.cex=.7,vertex.frame.color="gray",vertex.label=nodes_Da$name,layout=layout_with_fr)
#plot(net, vertex.label.font=1,vertex.label.cex=.7,vertex.label=nodes_Da$name,layout=layout_on_sphere)
#plot(net, vertex.label.font=1,vertex.label.cex=.7,vertex.label=nodes_Da$name,layout=layout_with_kk)

# the color below is the same order and the color as the C net
legend(x=1, y=0.35, c("aminocoumarin","aminoglycoside",
                         "bacitracin","fosfomycin",
                         "fosmidomycin","genus",
                         "glycopeptide","kasugamycin",
                         "macrolide-lincosamide-streptogramin","multidrug",
                         "mupirocin","peptide",
                         "polymyxin","quinolone",
                         "rifampin","tetracycline","trimethoprim"), pch=21,col="#777777",
       pt.bg=colrs, pt.cex=1.2, cex=1, bty="n", ncol=1,y.intersp=0.4)
par(family="Times")
tkid <- tkplot(net) 
#tiff('net_D.tiff', units="in", width=10, height=8, res=300)
dev.off()
# ----------This part will create the nodes name taht existed in the network of both D and C, is the appendix in the manussript
#setwd("/Volumes/GoogleDrive/My Drive/research/shogun_metagenome/annotation_results_megahit_C/")
nodes_name_C <- nodes_Ca %>% filter(phylumfamily=='genus') %>% mutate(name_C=name) %>% select(genus,name_C) 
nodes_name_D <- nodes_Da %>% filter(phylumfamily=='genus') %>% mutate(name_D=name) %>% select(genus,name_D) 
node_C_D <- full_join(nodes_name_C,nodes_name_D,by='genus')
head(node_C_D)
#write.csv(node_C_D,'name_nodes_D_C.csv',row.names = F)

# --------------------------------------join D table and C table to see how many bacteria and ARG are connected------------------
head(link_Ca)
head(link_Da)
colnames(link_Ca) <- c("from",'to','x-c')
colnames(link_Da) <- c("from",'to','x-d')
m_CD <- merge(link_Ca, link_Da, by=c("from","to"))
#write.csv(m_CD,"CD-common.csv",row.names = F)

# ----
C1 <- as.data.frame(table(link_Ca$from))
D1 <- as.data.frame(table(link_Da$from))
dim(C1)
dim(D1)
colnames(C1) <- c('bacteria','num_C')
colnames(D1) <- c('bacteria','num_D')
CD_b <- inner_join(C1,D1,by='bacteria')
dim(CD_b)
# ----
link_Ca$concat <- paste(link_Ca$from,link_Ca$to, sep=';') 
link_Da$concat <- paste(link_Da$from,link_Da$to, sep=';') 
C_link <- link_Ca[,c(3:4)]
colnames(C_link) <- c("num_C",'concat')
D_link <- link_Da[,c(3:4)]
colnames(D_link) <- c("num_D",'concat')

CD_link <- inner_join(C_link,D_link,by='concat')
CD_link <- CD_link[,c(2,1,3)]
getwd()
#write.csv(CD_link,"CD_link.csv",row.names = F)
split_CDlink <- cSplit(CD_link,"concat",';')
split_CDlink1 <- as.data.frame(table(split_CDlink[,3])) #55 bacteria in commom of CD

# check number of the bacteria respectively in C and D
CD_link1 <- full_join(C_link,D_link,by='concat')
CD_link1 <- CD_link1[,c(2,1,3)]
head(CD_link1)
C_link1 <- CD_link1[,c(1:2)]
D_link1 <- CD_link1[,c(1,3)]
#write.csv(CD_link,"CD_link.csv",row.names = F)
split_C_link1 <- cSplit(C_link1,"concat",';')
split_D_link1 <- cSplit(D_link1,"concat",';')
split_C_link1 <- as.data.frame(table(split_C_link1[,2])) #115 bacteria in commom of C
split_D_link1 <- as.data.frame(table(split_D_link1[,2])) #115 bacteria in commom of D

# C and D 
head(link_Ca)
#-----------16s_contigs_persistent tax-----------------------------
library(readxl)
library(splitstackshape)
X16s_contigs_persistent <- read_excel("G:/My Drive/research/shogun_metagenome/16s_contigs_persistent.xlsx")
#contig_16s <- cSplit(X16s_contigs_persistent,"tax",sep=';')
contig_16s <- as.data.frame(X16s_contigs_persistent) 
typeof(contig_16s)
contig_16s <- as.data.frame(contig_16s[,-c(1:2)])
table(contig_16s$tax_03)
is.data.frame(contig_16s)
dim(contig_16s)
head(contig_16s)
str(contig_16s)
# if row1 is ended with group, then row2 is not ended with group => row1=row2 k=k+1
# if row1 is ended with group, row2 is also ended with group => row1=row3 k=k+2
for (i in 1:nrow(contig_16s)){
  if (grepl('group$',contig_16s[i,1]) & !(grepl('group$',contig_16s[i,2]))){
    for (k in 2:ncol(contig_16s)) {
      contig_16s[i,k-1] <- contig_16s[i,k]}
  }
  }
contig_16s[is.na(contig_16s)] <- 'unknown'

for (i in 1:nrow(contig_16s)){
  if (grepl('group$',contig_16s[i,1]) & (grepl('group$',contig_16s[i,2]))){
    for (k in 3:ncol(contig_16s)) {
      contig_16s[i,k-2] <- contig_16s[i,k]}
  }
}

for (i in 1:nrow(contig_16s)){
  if (grepl('subdivisions$',contig_16s[i,2])){
    for (k in 3:ncol(contig_16s)) {
      contig_16s[i,k-1] <- contig_16s[i,k]}
  }
}

colnames(contig_16s) <- c('phylum','class','order','family','genus','species')
table(contig_16s$phylum)
table(contig_16s$class)
table(contig_16s$order)
table(contig_16s$family)
genus_t <- as.data.frame(table(contig_16s$genus)) # top persistent bacteria
genus_t1 <- genus_t[order(genus_t$Freq),]

