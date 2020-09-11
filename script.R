Sys.setenv(PATH = paste(Sys.getenv("PAH"),"/home/jerome/anaconda3/envs/bacterio/bin/",sep=':'))
system2(command='conda',args = 'list')

library(dplyr)
library(Biostrings)
library(syngulon)





####### download accession list
species= read.csv("/home/jerome/Documents/1-fichier/37-CTMA-SINGULON-AMARTIN/Analysis/01-selected.species/bacteria.n.csv")
species = species$Organism[species$n>20]
dir.create("01-accession-list/")
download.accession.NCBI(species = species,title = "plasmid",accessionDir = "01-accession-list/",index = 102)

###########Download genome et annotation
species= read.csv("/home/jerome/Documents/1-fichier/37-CTMA-SINGULON-AMARTIN/Analysis/01-selected.species/bacteria.n.csv")
species = species$Organism[species$n>20]
dir.create("02-annotation/")
dir.create("03-genome/")
dl.annot.genome(species = species,NmaxPlasmid = 10000,annotationDir = "02-annotation/",genomeDir = "03-genome/",accessionDir = "01-accession-list/",index=1)


######extract sequences from annotation
species= read.csv("/home/jerome/Documents/1-fichier/37-CTMA-SINGULON-AMARTIN/Analysis/01-selected.species/bacteria.n.csv")
species = species$Organism[species$n>20]
dir.create("06-extracted-sequences-annotation-notrans")
currentannotationfull = extract.sequences.annotation(species = species,annotationDir = "02-annotation/",genomeDir = "03-genome/",outDir = "06-extracted-sequences-annotation-notrans/")

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
dir.create("99-results/")
write.csv(result.blast.annotation,'99-results/result.blast.annotation.csv')





########## au propre pour otus
espece = list.dirs("08-colicinV/",full.names = F,recursive = F)
for (i in 1:length(espece)) {
  fasta.list.newspecies <- list.files(paste0("08-colicinV/",espece[i],'/'),full.names = T,recursive = T)
  fasta = list.files(paste0("08-colicinV/",espece[i],'/'),full.names = T)
  fasta.list.newspecies <- fasta.list.newspecies[grep('fasta',fasta.list.newspecies)]
  genename <- unique(basename(fasta.list.newspecies))
  for(j in 1:length(genename))
  {
    current.sequence <- readDNAStringSet(fasta.list.newspecies[grep(genename[j],fasta.list.newspecies)])
    writeXStringSet(current.sequence,paste0("08-colicinV/",espece[i],"/",genename[j]))
  }
  diri = list.dirs(paste0("08-colicinV/",espece[i],'/'),full.names = T)
  diri = diri[2:length(diri)]
  unlink(diri,recursive = T)
}
fasta.list=c()
for (i in 1:length(espece)) {
  fasta.list.newspecies <- list.files(paste0("08-colicinV/",espece[i],'/'),full.names = T)
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
espece = basename(list.dirs("08-colicinV/",recursive = F))
geneDir <- "08-colicinV/"
dir.create("09-otus")
outDir <- "09-otus/"
collicin <- read.csv('Colicin_nts.csv',sep = ';')
otus(species = espece,collicin = collicin,geneDir = geneDir,outDir = outDir)

#########plot otus phylo
bacteria.table = read.csv("/Key/01-selected.species/bacteria.table.csv")
phylo.from.otus(otusDir = "09-otus/",bacteria.table = bacteria.table)





######heat map blast
result.blast = read.csv("99-results/result.blast.csv",row.names = 1)
result.blast = data.frame(result.blast)
result.blast = filter_all(result.blast, any_vars(. !=""))
presencegenome <- data.frame(apply(result.blast,2,function(x) as.numeric(x!='')))
presencegenome$bss = rep(0,length(presencegenome$bss))#probleme avec les cna donc mis tout a zero car il n'y en a pas
rownames(presencegenome)=rownames(result.blast)
library(gplots)
matrix.to.plot <- as.matrix(presencegenome)
my_palette <- colorRampPalette(c("white", "gray", "black"))(n = 299)
pdf('99-results/heatmap.pdf',width=18,height = 100)
heatmap.2(matrix.to.plot,dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,508),sepwidth=c(0.015,0.015),
          labRow =rownames(rownames(result.blast)),margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5,breaks=seq(0,1,length.out = 300) )
dev.off()
hitgenome = rep(0,length(names(presencegenome)))
for (i in 1:length(names(presencegenome))) {
  hitgenome[i]=sum(presencegenome[i])
}


###calcul de presence des genes
coli = presencegenome%>%filter(.,grepl("Citrobacter",rownames(presencegenome)))
hitcoli = rep(0,length(names(coli)))
for (i in 1:length(names(coli))) {
  hitcoli[i]=sum(coli[i])
}
hitcoli
  #######heat map blast annotation
result.blast.annotation = read.csv("99-results/result.blast.annotation.csv",row.names = 1)
depart = dim(result.blast.annotation)[1]
result.blast.annotation = data.frame(result.blast.annotation)
result.blast.annotation = filter_all(result.blast.annotation, any_vars(. !=""))
presenceblast <- data.frame(apply(result.blast.annotation,2,function(x) as.numeric(x!='')))
arrive = dim(presenceblast)[1]
presenceblast$bss = rep(0,length(presenceblast$bss))#probleme avec les bss donc mis tout a zero car il n'y en a pas
pourcentagegenome = arrive/depart 
rownames(presenceblast)=rownames(result.blast.annotation)
library(gplots)
matrix.to.plot <- as.matrix(presenceblast)
my_palette <- colorRampPalette(c("white", "gray", "black"))(n = 299)
pdf('99-results/heatmap.annotation.pdf',width=18,height = 100)
heatmap.2(matrix.to.plot,dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,400),sepwidth=c(0.015,0.015),
          labRow =rownames(rownames(result.blast.annotation)),margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5,breaks=seq(0,1,length.out = 300) )
dev.off()

hit = rep(0,length(names(presenceblast)))
for (i in 1:length(names(presenceblast))) {
  hit[i]=sum(presenceblast[i])
}

###calcul de presence des genes
coli = presenceblast%>%filter(.,grepl("coli",rownames(presenceblast)))
hitcoli = rep(0,length(names(coli)))
for (i in 1:length(names(coli))) {
  hitcoli[i]=sum(coli[i])
}
hitcoli


############ calcul du pourcentage entre le nombre d'annotation et le nombre de gene trouver avec blast
nombre =0
for (i in 1:dim(presenceblast)[1]) {
  for (j in 1:dim(presenceblast)[2]) {
    if (presenceblast[i,j]==1) {
      nombre= nombre+1
    }
  }
}

pourcentageannotation = nombre/dim(currentannotationfullnom)[1]


######### heatmap pour les especes
bacteria.table= read.csv("/home/jerome/Documents/1-fichier/37-CTMA-SINGULON-AMARTIN/Analysis/01-selected.species/bacteria.table.csv")
result.blast.annotation = read.csv("99-results/result.blast.annotation.csv",row.names = 1)
result.blast.annotation = data.frame(result.blast.annotation)
result.blast.annotation = filter_all(result.blast.annotation, any_vars(. !=""))
presenceblast <- data.frame(apply(result.blast.annotation,2,function(x) as.numeric(x!='')))
presenceblast$bss = rep(0,length(presenceblast$bss))#probleme avec les bss donc mis tout a zero car il n'y en a pas
rownames(presenceblast)=rownames(result.blast.annotation)


correspondance.organism.subgroup <- bacteria.table%>%group_by(Organism) %>% dplyr::count(Organism, SubGroup) %>% dplyr::slice(which.max(n)) %>% dplyr::rename(species=Organism)
species <- unlist(lapply(strsplit(rownames(presenceblast),split='/'),function(x) x[[2]]))
presence <- data.frame(species,presenceblast)
rownames(presence) <- rownames(presenceblast)
##### summary at bacteria species

presence.summary.n <- presence %>% group_by(species) %>% dplyr::summarise(n=n())
presence.summary.mean <- presence %>% group_by(species) %>% select(-species) %>% dplyr::summarise_all(.funs = list(MEAN = ~ round(mean(x = .,na.rm=T),2)))
presence.summary <- full_join(presence.summary.n,presence.summary.mean,by='species')
presence.summary <- presence.summary %>% rename_all(funs(gsub("_MEAN", "", .)))

write.csv(presence.summary,'99-results/presence.summary.csv')





######### arbre phylogenetic
bacteria.table= read.csv("/home/jerome/Documents/1-fichier/37-CTMA-SINGULON-AMARTIN/Analysis/01-selected.species/bacteria.table.csv")

phylo.all(bacteria.table = bacteria.table,Dir = "06-extracted-sequences-annotation-notrans/")

######## distance 2D
fileDir = "06-extracted-sequences-annotation-notrans/"
filelist = list.files(fileDir,full.names = T,recursive = T)
ADN = readDNAStringSet(filelist)
distance <- stringDist(ADN)
pcoa <- cmdscale(distance,k = 2)
pcoa = jitter(pcoa,factor = 200)
mycolor <- rep('black',length(rownames(pcoa)))
presence = rownames(presenceblastnom)
presence = unlist(lapply(strsplit(presence,split='/'),function(x) x[3]))
for (i in 1:length(presence)) {
  ligne = presenceblastnom%>% filter(.,grepl(presence[i],rownames(presenceblastnom)))
  if (ligne$E2 == 1 |ligne$E3 == 1 | ligne$E5 == 1 |ligne$E6 == 1 |ligne$E7 == 1 |ligne$E8 == 1 |ligne$E9 == 1  ) {
    mycolor[rownames(pcoa)==presence[i]] = "red"
  }
  if (ligne$cma == 1 ) {
    mycolor[rownames(pcoa)==presence[i]] = "yellow"
  }
  if (ligne$cta == 1 |ligne$C5 == 1 | ligne$cka == 1  ) {
    mycolor[rownames(pcoa)==presence[i]] = "blue"
  }
  if (ligne$cea == 1  ) {
    mycolor[rownames(pcoa)==presence[i]] = "darkorange"
  }
  if (ligne$cna == 1  ) {
    mycolor[rownames(pcoa)==presence[i]] = "azure4"
  }
  if (ligne$csa == 1  ) {
    mycolor[rownames(pcoa)==presence[i]] = "coral4"
  }
}
i=1
while (i<=length(rownames(pcoa))) {
  current = currentannotationfull%>% filter(.,grepl(strsplit(rownames(pcoa)[i],'-')[[1]][1],currentannotationfull$name))
  j=1
  while (j<=length(rownames(current))) {
    currentnow = current%>%slice(j)
    #if (!isEmpty(grep("colicin V",currentnow$product))) {
     # mycolor[i+j] = "blueviolet"
    #}
    if ("colicin V"==currentnow$product | "colicin V family bacteriocin" ==currentnow$product) {
    mycolor[i+j] = "blueviolet"
    }
    
    j=j+1
  }
  i = i+length(rownames(current))
  print(i)
}
  
  



x <- pcoa[, 1]
y <- -pcoa[, 2] # reflect so North is at the top
## note asp = 1, to ensure Euclidean distances are represented correctly
plot(x, y,  asp = 1, main = "cmdscale(sequences)",col = mycolor)
legend('topright',legend = c('E','cma','cta|C5|cka','cea'),fill=c('red','yellow','blue','darkorange'))
text(x, y, rownames(pcoa), cex = 0.6)


##########analyser les produits
tableannotation = currentannotationfull%>% group_by(product)%>% summarise(n())%>%rename(nombre = `n()`)

tableannotation = tableannotation%>%filter(tableannotation$nombre>10)












######pour les noms
species= read.csv("/home/jerome/Documents/1-fichier/37-CTMA-SINGULON-AMARTIN/Analysis/01-selected.species/bacteria.n.csv")
species = species$Organism[species$n>20]
dir.create("07-extracted-sequences-annotation-notrans-nom")
currentannotationfullnom = extract.sequences.annotation.nom(species = species,annotationDir = "02-annotation/",genomeDir = "03-genome/",outDir = "07-extracted-sequences-annotation-notrans-nom/")


collicin = read.csv("Colicin_nts.csv",sep = ";")
sequences = collicin$Nucleotide.Sequence
sequences = DNAStringSet(sequences)
names(sequences) = collicin$genename
sequences = sequences[1:17]
writeXStringSet(sequences,"sequences.fasta")
ngenes <- length(sequences)
genomelist <- list.files('07-extracted-sequences-annotation-notrans-nom/',full.names = T,recursive = T)
genomelist<- genomelist[grep('.fasta',genomelist)]
ngenomes <- length(genomelist)
result.blast.annotationom <- matrix(nrow=ngenomes,ncol=ngenes)
rownames(result.blast.annotationom) <- gsub(gsub(genomelist,pattern = '07-extracted-sequences-annotation-notrans-nom/',replacement = ''),pattern = ".fasta",replacement = "")
genes <-names(sequences)
gene.levels <- levels(factor(genes))
colnames(result.blast.annotationom) <- gene.levels


for(i in 1:ngenomes)
{
  result.blast.annotationom[i,] <- screenBlast(reference =  'sequences.fasta',querry = genomelist[i],min.pc.ident = 80 ,min.pc.length=50)

  print(i)
}
dir.create("99-results/")
write.csv(result.blast.annotationom,'99-results/result.blast.annotationom.csv')

result.blast.annotationom = read.csv("99-results/result.blast.annotationom.csv",row.names = 1)
depart = dim(result.blast.annotationom)[1]
result.blast.annotationom = data.frame(result.blast.annotationom)
result.blast.annotationom = filter_all(result.blast.annotationom, any_vars(. !=""))
presenceblastnom <- data.frame(apply(result.blast.annotationom,2,function(x) as.numeric(x!='')))
arrive = dim(presenceblastnom)[1]
presenceblastnom$bss = rep(0,length(presenceblastnom$bss))#probleme avec les bss donc mis tout a zero car il n'y en a pas
pourcentagegenome = arrive/depart 
rownames(presenceblastnom)=rownames(result.blast.annotationom)
library(gplots)
matrix.to.plot <- as.matrix(presenceblastnom)
my_palette <- colorRampPalette(c("white", "gray", "black"))(n = 299)
pdf('99-results/heatmap.annotationom.pdf',width=18,height = 100)
heatmap.2(matrix.to.plot,dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,400),sepwidth=c(0.015,0.015),
          labRow =rownames(rownames(result.blast.annotationom)),margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5,breaks=seq(0,1,length.out = 300) )
dev.off()
