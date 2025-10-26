library("biomaRt")

set <- c("may2015","oct2014","feb2014")

HOST<-"https://may2015.archive.ensembl.org"

datasets <- listDatasets(ensembl)0

hs <- 'hsapiens_gene_ensembl'
MART<-"ENSEMBL_MART_ENSEMBL"
ensembl=useMart(MART, dataset = hs)


#LOADING GENELIST

filename <- 'BLUE_MYC_CLUST'
ext<-'.txt'
table <- paste0("C:/Users/newma/OneDrive/Desktop/Placement/",filename,ext)
GeneList <- read.table(table)


#MAPPING GENE SYMBOLS TO ENSEMBL

hgnc_symbol_list<- GeneList$V1

ensembl_list <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','external_gene_name'), 
                      filters='hgnc_symbol',
                      values=hgnc_symbol_list, 
                      mart=ensembl)

unchanged <- intersect(GeneList$V1,ensembl_list$ensembl_gene_id)

title<-paste0(filename,'_list.txt',sep='')
write.table(ensembl_list,title,sep="\t",row.names=FALSE,quote=FALSE) 



