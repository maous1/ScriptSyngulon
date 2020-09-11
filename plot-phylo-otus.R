#' A function to plot phylogenetic tree of all genes based on fasta file in the otus directory
#'
#' @param otusDir the directory with the otus
#' @param bacteria.table bacteria table
#'
#' @return
#' @export
#'
#' @examples
phylo.from.otus <- function(otusDir,bacteria.table)
{
  library(dplyr)
  library(Biostrings)
  library(muscle)
  library(phangorn)
  library(ggtree)
  correspondance.organism.subgroup <- bacteria.table%>%group_by(Organism) %>% dplyr::count(Organism, SubGroup) %>% dplyr::slice(which.max(n)) %>% dplyr::rename(species=Organism)
  filelist <- list.files(otusDir,full.names = T)
  filelist <- filelist[grep(pattern = '.fasta',x =filelist )]
  
  genename <- gsub(basename(filelist),pattern = '.fasta',replacement = '')
  ngenes <- length(genename)
  for(i in 1:ngenes)
  {
    sequences <- readDNAStringSet(filelist[i])
    if (length(sequences)>1) {
      sequences <- sequences[sort.list(names(sequences))]
      currentname <- names(sequences)
      newname <- currentname
      newname.part1 <- paste(unlist(lapply(strsplit(newname,split=':'),function(x) x[[3]])),unlist(lapply(strsplit(newname,split=':'),function(x) x[[2]])),sep='_')
      
      newname.part2 <- unlist(lapply(strsplit(newname,split=':'),function(x) x[[1]]))
      newname.part2 <- data.frame(species=newname.part2)
      newname.part2 <- left_join(newname.part2,correspondance.organism.subgroup,by='species')
      newname.part2 <- paste(newname.part2$SubGroup,newname.part2$species,sep='_')
      newname <- paste(newname.part1,newname.part2,sep='_')
      
      names(sequences) <- newname
      maxlength <- max(width(sequences))
      
      align.muscle <- muscle::muscle(sequences)
      dist1 <- stringDist(as(align.muscle,"DNAStringSet"), method="hamming")
      dist1 <- 100*dist1/maxlength
      mytree1 <- upgma(dist1)
      
      nseqeuences <- length(sequences)
      
      pdf(paste0('99-results/phylogenetic.otus.',genename[i],'.pdf'),width=10,height = 10*log10(nseqeuences))
      p <- ggtree(mytree1)
      p <- p  + geom_tiplab(offset=0) + xlim(NA, 150) +geom_treescale(0.05,-2.5,width=10,fontsize = 2,linesize = 0.5)
      plot(p)
      dev.off()
    }
  }
}
