Sys.setenv(PATH = paste(Sys.getenv("PAH"),"/home/jerome/anaconda3/envs/bacterio/bin/",sep=':'))
system2(command='conda',args = 'list')

library(dplyr)
library(Biostrings)
library(syngulon2)


species= read.csv("/media/jerome/Seagate Expansion Drive//Key/01-selected.species/bacteria.n.csv")
species = species$Organism[species$n>20]
species2 = gsub("_"," ",species)
dir.create("01-accession-list")


####### download accession list
for (i in 1:length(species)) {
  
  demo.search <- esearch(term = paste0(species2[i],"[orgn] and plasmid[title]"), db = 'nuccore', usehistory = TRUE) #search
  accessions <- efetch(demo.search, rettype = "acc",retmode = "text",outfile= paste0('01-accession-list/',species[i],".csv"))#fetch accessions
  accession  <- read.csv(paste0('01-accession-list/',species[i],".csv"),header=F,stringsAsFactors = F) 
  accession$V1 <- unlist(lapply(strsplit(accession$V1,split='\\.'),function(x) x[[1]]))
  write.csv(accession,paste0('01-accession-list/',species[i],".csv"),row.names = F)
  print(i)
  
}

###########Download genome et annotation
nombre = 1000
for (j in 1:length(species)) {
  dir.create(paste0("02-annotation/",species[j]))
  dir.create(paste0("03-genome/",species[j]))
  accession <- read.csv(paste0('01-accession-list/',species[j],".csv"),stringsAsFactors = F)
  accession  <- accession$V1
  if (length(accession)>nombre) {
    accession <- accession[sample(1:length(accession), nombre, replace=F)]
  }
  N.accession <- length(accession)
  for(i in 1:N.accession)
  {
    seq <- try(read.GenBank(accession[i]),silent = T)
    current.annotation <- try(getAnnotationsGenBank(access.nb=accession[i], quiet = TRUE))
    if (class(seq)=="DNAbin" & class(current.annotation)=="data.frame"& (all(is.element(c('start','end','type','product'),colnames(current.annotation))))) {
      write.FASTA(seq,paste0('03-genome/',species[j],"/",accession[i],'.fasta'))
      seq <- readDNAStringSet(paste0('03-genome/',species[j],"/",accession[i],'.fasta'))
      writeXStringSet(seq,paste0('03-genome/',species[j],"/",accession[i],'.fasta'))
      current.annotation <- tibble(current.annotation)
      current.annotation <- current.annotation %>% select(start,end,type,gene,product)
      write.csv(current.annotation,paste0('02-annotation/',species[j],"/",accession[i],'.csv'),row.names = F)
    }
    
    print(i)
  }
}




############screenblast sur les génomes

collicin = read.csv("Colicin_nts.csv",sep = ";")
sequences = collicin$Nucleotide.Sequence
sequences = DNAStringSet(sequences)
names(sequences) = collicin$genename
sequences = sequences[1:17]
writeXStringSet(sequences,"sequences.fasta")
ngenes <- length(sequences)
genomelist <- list.files('03-genome/',full.names = T,recursive = T)
genomelist<- genomelist[grep('.fasta',genomelist)]
ngenomes <- length(genomelist)
result.blast <- matrix(nrow=ngenomes,ncol=ngenes)
rownames(result.blast) <- gsub(gsub(genomelist,pattern = '03-genome/',replacement = ''),pattern = ".fasta",replacement = "" )
genes <- names(sequences)
gene.levels <- levels(factor(genes))
colnames(result.blast) <- gene.levels


unlink("04-extracted-sequences-from-blast",recursive = T)
dir.create("04-extracted-sequences-from-blast")
for(i in 1:ngenomes)
{
  result.blast[i,] <- screenBlastkey(reference =  'sequences.fasta',querry = genomelist[i],min.pc.ident = 80 ,min.pc.length=50,geneDir = "04-extracted-sequences-from-blast/")
  species = strsplit(rownames(result.blast)[i],"/")[[1]][2]
  taille = readDNAStringSet(paste0("03-genome/",species,"/",basename(genomelist[i])))
  taille = length(taille[[1]])
  rownames(result.blast)[i]=paste0(taille,rownames(result.blast)[i])
  print(i)
}
write.csv(result.blast,'99-results/result.blast.csv')


########## screenblast sur les annotations
collicin = read.csv("Colicin_nts.csv",sep = ";")
sequences = collicin$Nucleotide.Sequence
sequences = DNAStringSet(sequences)
names(sequences) = collicin$genename
sequences = sequences[1:17]
writeXStringSet(sequences,"sequences.fasta")
ngenes <- length(sequences)
genomelist <- list.files('05-extracted-sequences-from-annotation/',full.names = T,recursive = T)
genomelist<- genomelist[grep('.fasta',genomelist)]
ngenomes <- length(genomelist)
result.blast.annotation <- matrix(nrow=ngenomes,ncol=ngenes)
rownames(result.blast.annotation) <- gsub(gsub(genomelist,pattern = '05-extracted-sequences-from-annotation/',replacement = ''),pattern = ".fasta",replacement = "")
genes <-names(sequences)
gene.levels <- levels(factor(genes))
colnames(result.blast.annotation) <- gene.levels


for(i in 1:ngenomes)
{
  result.blast.annotation[i,] <- screenBlast(reference =  'sequences.fasta',querry = genomelist[i],min.pc.ident = 80 ,min.pc.length=50)
  species = strsplit(rownames(result.blast.annotation)[i],"/")[[1]][2]
  taille = readDNAStringSet(paste0("03-genome/",species,"/",basename(genomelist[i])))
  taille = length(taille[[1]])
  rownames(result.blast.annotation)[i]=paste0(taille,rownames(result.blast.annotation)[i])
  print(i)
}
write.csv(result.blast.annotation,'99-results/result.blast.annotation.csv')


########## au propre pour otus
espece = list.dirs("04-extracted-sequences-from-blast/",full.names = F,recursive = F)
for (i in 1:length(espece)) {
  fasta.list.newspecies <- list.files(paste0("04-extracted-sequences-from-blast/",espece[i],'/'),full.names = T,recursive = T)
  fasta = list.files(paste0("04-extracted-sequences-from-blast/",espece[i],'/'),full.names = T)
  fasta.list.newspecies <- fasta.list.newspecies[grep('fasta',fasta.list.newspecies)]
  genename <- unique(basename(fasta.list.newspecies))
  for(j in 1:length(genename))
  {
    current.sequence <- readDNAStringSet(fasta.list.newspecies[grep(genename[j],fasta.list.newspecies)])
    writeXStringSet(current.sequence,paste0("04-extracted-sequences-from-blast/",espece[i],"/",genename[j]))
  }
  diri = list.dirs(paste0("04-extracted-sequences-from-blast/",espece[i],'/'),full.names = T)
  diri = diri[2:length(diri)]
  unlink(diri,recursive = T)
}
fasta.list=c()
for (i in 1:length(espece)) {
  fasta.list.newspecies <- list.files(paste0("04-extracted-sequences-from-blast/",espece[i],'/'),full.names = T)
  fasta.list.newspecies <- fasta.list.newspecies[grep('fasta',fasta.list.newspecies)]
  fasta.list <- c(fasta.list,fasta.list.newspecies)
}
fileinfo <- file.info(fasta.list)
fasta.list <- fasta.list[fileinfo$size>0]
genename <- unique(basename(fasta.list))
for(i in 1:length(genename))
{
  current.sequence <- readDNAStringSet(fasta.list[grep(genename[i],fasta.list)])
  writeXStringSet(current.sequence,paste0("04-extracted-sequences-from-blast/",genename[i]))
}




######otus
espece = basename(list.dirs("04-extracted-sequences-from-blast",recursive = F))
geneDir <- "04-extracted-sequences-from-blast"
outDir <- "05-otus/"
collicin <- read.csv('Colicin_nts.csv',sep = ';')
otus(species = espece,collicin = collicin,geneDir = geneDir,outDir = outDir)

#########plot otus phylo
bacteria.table = read.csv("/Key/01-selected.species/bacteria.table.csv")
phylo.from.otus(otusDir = "04-extracted-sequences-from-blast/",bacteria.table = bacteria.table)





######heat map blast
result.blast = read.csv("99-results/result.blast.csv",row.names = 1)
result.blast = data.frame(result.blast)
result.blast = filter_all(result.blast, any_vars(. !=""))
presenceblast <- data.frame(apply(result.blast,2,function(x) as.numeric(x!='')))
presenceblast$cna = rep(0,length(presenceblast$cna))
rownames(presenceblast)=rownames(result.blast)
library(gplots)
matrix.to.plot <- as.matrix(presenceblast)
my_palette <- colorRampPalette(c("white", "gray", "black"))(n = 299)
pdf('99-results/heatmap.pdf',width=18,height = 30)
heatmap.2(matrix.to.plot,dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,150),sepwidth=c(0.015,0.015),
          labRow =rownames(rownames(result.blast)),margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5,breaks=seq(0,1,length.out = 300) )
dev.off()


#######heat map blast annotation
result.blast.annotation = read.csv("99-results/result.blast.annotation.csv",row.names = 1)
result.blast.annotation = data.frame(result.blast.annotation)
result.blast.annotation = filter_all(result.blast.annotation, any_vars(. !=""))
presenceblast <- apply(result.blast.annotation,2,function(x) as.numeric(x!=''))
rownames(presenceblast)=rownames(result.blast.annotation)
library(gplots)
matrix.to.plot <- as.matrix(presenceblast)
my_palette <- colorRampPalette(c("white", "gray", "black"))(n = 299)
pdf('99-results/heatmap.annotation.pdf',width=18,height = 30)
heatmap.2(matrix.to.plot,dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,150),sepwidth=c(0.015,0.015),
          labRow =rownames(rownames(result.blast.annotation)),margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5,breaks=seq(0,1,length.out = 300) )
dev.off()