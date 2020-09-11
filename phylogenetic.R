library(dplyr)
library(Biostrings)
library(muscle)
library(phangorn)
library(ggtree)
filelist <- list.files(otusDir,full.names = T,recursive = T)
filelist <- filelist[grep(pattern = '.fasta',x =filelist )]

sequences <- readDNAStringSet(filelist)

    names(sequences) <- newname
    maxlength <- max(width(sequences))
    
    align.muscle <- muscle::muscle(sequences)
    dist1 <- stringDist(as(align.muscle,"DNAStringSet"), method="hamming")
    dist1 <- 100*dist1/maxlength
    mytree1 <- upgma(dist1)
    
    nseqeuences <- length(sequences)
    
    
    
    tiplabel <- mytree1$tip.label
    mycolor <- rep('black',length(tiplabel))
    presence = rownames(presenceblastnom)
    presence = unlist(lapply(strsplit(presence,split='/'),function(x) x[3]))
    mycolor[is.element(tiplabel,presence)] = "green"
    
    pdf('99-results/phylogenetic.pdf',width=10,height = 60*log10(nseqeuences))
    p <- ggtree(mytree1)
    p <- p  + geom_tiplab(color=mycolor,offset=0) + xlim(NA, 150) +geom_treescale(0.05,-2.5,width=10,fontsize = 2,linesize = 0.5)
    plot(p)
    dev.off()
    
    