library(readxl)
require(reshape)
library(dplyr)
library(tidyr)
setwd("/Users/luyang/Downloads/annotation_results_megahit_C")
D_arg_Prodigal <- read_excel("matched_ARG_orf_D.fa.emapper.annotations.xlsx")
C_arg_Prodigal <- read_excel("matched_ARG_orf_C.fa.emapper.annotations.xlsx")
head(C_arg_Prodigal)
C31 <- C_arg_Prodigal[,c(1,12)]  #KEGG 7,COG 11
D31 <- D_arg_Prodigal[,c(1,12)] 
colnames(C31) <- c("contig","COG")
colnames(D31) <- c("contig","COG")
C31_a <- C31[,2]
D31_a <- D31[,2]
C <- as.data.frame(table(C31_a))
D <- as.data.frame(table(D31_a))
colnames(C) <- c("COG","C")
colnames(D) <- c("COG","D")
CD <- full_join(C,D,'COG')
m_CD <- reshape::melt(CD,"COG")
library(ggplot2)
library(scales)
ggplot(m_CD, aes(factor(COG), value, fill = variable)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1")+
  labs(x = 'COG Category',y='Genes')+
  theme(text = element_text(size = 12, family = "Times New Roman",color = "black"),
        axis.ticks.length=unit(0.1, "cm"),
        axis.text.x = element_text(color="black", size = 8,face = 'bold'),
        axis.text.y = element_text(color="black", size = 10,face = 'bold'),
        legend.text=element_text(size=10,face = 'italic'),
        legend.title=element_text(size=10,face = 'bold'),
        legend.position = c(.25, .85))+
  guides(fill = guide_legend(ncol =1,title = "Treatment"))


#---------KEGG------------------
ko_list <- C_arg_Prodigal[,7]
ko_list$name <- NA
ko_list$pathway <- NA
head(ko_list)
unique_ko_list <- unique(ko_list)
nrow(unique_ko_list)
unique_ko_list <- filter(unique_ko_list,!is.na(unique_ko_list$`KEGG KO`))
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("KEGGREST")
library(KEGGREST)
for (i in 1:nrow(unique_ko_list)){
  tryCatch({
    x <- as.character(unique_ko_list[i,1])
    y <- unlist(strsplit(x,','))
    if (length(y)==1) {
      query <- keggGet(y)
      unique_ko_list[i,2] <- toString(query[[1]]$NAME)
      unique_ko_list[i,3] <- toString(query[[1]]$PATHWAY)
    }else {
      name_list <- list()
      path_list <- list()
      query <- keggGet(y)
      for (j in 1:length(query)){
        name_list[[j]] <- query[[j]]$NAME
        path_list[[j]] <- query[[j]]$PATHWAY
      }
      unique_ko_list[i,2] <- paste(unlist(name_list), collapse = ';')
      unique_ko_list[i,3] <- paste(unlist(path_list), collapse = ";")
    }}, error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
}
#write.csv(unique_ko_list,"unique_ko_list_C.csv")

ko_list <- D_arg_Prodigal[,7]
ko_list$name <- NA
ko_list$pathway <- NA
head(ko_list)
unique_ko_list <- unique(ko_list)
nrow(unique_ko_list)
unique_ko_list <- filter(unique_ko_list,!is.na(unique_ko_list$`KEGG KO`))
unique_ko_list <- read_excel("C:/Users/Lu Yang/Desktop/KEGG_C.xlsx")


for (i in 1:nrow(unique_ko_list)){
  tryCatch({
    x <- as.character(unique_ko_list[i,1])
    y <- unlist(strsplit(x,','))
    if (length(y)==1) {
      query <- keggGet(y)
      unique_ko_list[i,2] <- toString(query[[1]]$NAME)
      unique_ko_list[i,3] <- toString(query[[1]]$PATHWAY)
    }else {
      name_list <- list()
      path_list <- list()
      query <- keggGet(y)
      for (j in 1:length(query)){
        name_list[[j]] <- query[[j]]$NAME
        path_list[[j]] <- query[[j]]$PATHWAY
      }
      unique_ko_list[i,2] <- paste(unlist(name_list), collapse = ';')
      unique_ko_list[i,3] <- paste(unlist(path_list), collapse = ";")
    }}, error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
}

#write.csv(unique_ko_list,"unique_ko_list_D.csv",row.names=F)
#write.csv(unique_ko_list,"KEGG_C.csv",row.names=F)

unique_ko_list_D <- unique_ko_list
D1 <- cSplit(unique_ko_list_D, "pathway", ";")


setwd("G:/My Drive/research/shogun_metagenome/")
eggnog_C <- read_excel("eggnog_C.xlsx")
eggnog_D <- read_excel("eggnog_D.xlsx")

eggnog <- full_join(eggnog_C,eggnog_D,by="eggnog") 
eggnog$eggnog <- gsub('Multidrug','multidrug',eggnog$eggnog)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
eggnog$eggnog <- firstup(eggnog$eggnog)
eggnog[is.na(eggnog)] <- 0 
eggnog$eggnog <- gsub('Ec 1.1.1.22','UDP-glucose 6-dehydrogenase',eggnog$eggnog)
eggnog$eggnog <- gsub('Ec 2.3.1.81','Aminoglycoside N(3)-acetyltransferase',eggnog$eggnog)
eggnog$eggnog <- gsub('Ec 2.7.7.47','Streptomycin 3-adenylyltransferase',eggnog$eggnog)
eggnog$eggnog <- gsub('Inherit from bactNOG:','',eggnog$eggnog)
eggnog$eggnog <- gsub('Inherit from COG:','',eggnog$eggnog)
eggnog$eggnog <- gsub('Inherit from firmNOG:','',eggnog$eggnog)
eggnog$eggnog <- gsub('Major facilitator superfamily MFS_1','Major facilitator superfamily',eggnog$eggnog)
eggnog$eggnog <- gsub('Major Facilitator superfamily','Major facilitator superfamily',eggnog$eggnog)
eggnog$eggnog <- gsub('Major facilitator superfamily.*','Major facilitator superfamily',eggnog$eggnog)
eggnog$eggnog <- gsub('Penicillin-binding protein 2','Penicillin-binding protein',eggnog$eggnog)
eggnog$eggnog <- gsub('Resistance','Resistance protein',eggnog$eggnog)
eggnog$eggnog <- gsub('Resistance protein to gentamicin, tobramycin.*','Resistance protein',eggnog$eggnog)
eggnog$eggnog <- gsub('Translation elongation factor g','Translation elongation factor G',eggnog$eggnog)
eggnog$eggnog <- gsub('RRNA adenine N-6-methyltransferase','RRNA (Adenine-N(6)-)-methyltransferase',eggnog$eggnog)
eggnog$eggnog <- gsub('Aminoglycoside 6-adenylyltransferase.*','Aminoglycoside 6-adenylyltransferase',eggnog$eggnog)
eggnog$eggnog <- gsub(')-reductase','Reductase',eggnog$eggnog)
eggnog$eggnog <- gsub('transcriptional','Transcriptional',eggnog$eggnog)
eggnog$eggnog <- gsub('Cell wall formation.*','Cell wall formation',eggnog$eggnog)
eggnog$eggnog <- gsub('Major Facilitator','Major facilitator superfamily',eggnog$eggnog)
eggnog$eggnog <- gsub('Methionine sulfoxide reductase.*','Methionine sulfoxide reductase',eggnog$eggnog)
eggnog$eggnog <- gsub('Resistance protein protein','Resistance protein',eggnog$eggnog)

eggnog$eggnog <- gsub('Transcriptional.*','Transcriptional regulator',eggnog$eggnog)
eggnog$eggnog <- gsub('Penicillin binding protein.*','Penicillin-binding protein',eggnog$eggnog)
eggnog$eggnog <- gsub('RND efflux system.*','RND efflux system',eggnog$eggnog)
eggnog$eggnog <- gsub('Streptomycin 3-adenylyltransferase','Streptomycin 3"-adenylyltransferase',eggnog$eggnog)
eggnog$eggnog <- gsub('Udp-glucose 6-dehydrogenase','UDP-glucose 6-dehydrogenase',eggnog$eggnog)
eggnog$eggnog <- gsub('Provides the sole de novo source of dTMP for DNA biosynthesis.*','Thymidylate synthase',eggnog$eggnog)
eggnog$eggnog <- gsub('This enzyme is an effector of chloramphenicol resistance in bacteria.*','Chloramphenicol resistant',eggnog$eggnog)
eggnog$eggnog <- gsub('Streptomycin 3.*','Streptomycin 3"-adenylyltransferase',eggnog$eggnog)
eggnog$eggnog <- gsub('Phosphotransferase enzyme family','Phosphotransferase',eggnog$eggnog)
eggnog$eggnog <- gsub('Chloramphenicol O-acetyltransferase','Chloramphenicol acetyltransferase',eggnog$eggnog)


eggnog <- aggregate(.~eggnog, data=eggnog, FUN=sum)
is.data.frame(eggnog)
