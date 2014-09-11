setwd("~/Desktop/CRISPR-CRUNCHER")

options(stringsAsFactors = F)

# Analyze ABI sequencing file to determine whether a deletion

# install bioconductor packages
if(!all(c("sangerseqR", "seqinr") %in% rownames(installed.packages()))){
  source("http://bioconductor.org/biocLite.R")
  biocLite("sangerseqR")
  biocLite("seqinr")
  biocLite("org.Hs.eg.db")
}

library(sangerseqR)
library(seqinr)
library(org.Hs.eg.db)

# Human coding sequences taken from Ensembl CDS dump http://useast.ensembl.org/info/data/ftp/index.html

#human_sequences <- read.fasta('humanTranscript200bpFlanking.fasta')
#seq_attrs <- t(sapply(1:length(human_sequences), function(x){
#  seq_attr <- strsplit(getAnnot(human_sequences[[x]]), split = '[>\\|]')[[1]][-1]
#  data.frame(index = x, ensemblGene = seq_attr[1], common = seq_attr[3])
#}, simplify = T))
#seq_attrs <- as.data.frame(apply(seq_attrs, c(1,2), function(x){x[[1]]}))
#save(human_sequences, seq_attrs, file = "humanTranscript200bpFlanking.Rdata")
load("humanTranscript200bpFlanking.Rdata")

# IDconv <- read.delim('ensemblIDconv.txt') # ensembl to common taken from ensembl biomart
IUPAC <- read.delim('IUPAC.txt', header = F)
rownames(IUPAC) <- IUPAC[,1]; IUPAC <- IUPAC[,-1]; colnames(IUPAC) <- c("A", "C", "G", "T")
# Identify which sequences need to be deconvoluted - ABIfiles

ABIfiles <- list.files("ABIfiles")
an_ABIfile <- ABIfiles[1]
for(an_ABIfile in ABIfiles){
  
  # Find the sequence of the wild-type gene
  gene <- strsplit(an_ABIfile, split = '[-_]')[[1]][2]
  gene_attr <- seq_attrs[seq_attrs$common == gene,]
  
  if(nrow(gene_attr) == 0){
    gene_sequence <- NA
    warning("Did not find gene: ", gene)
  }else if(nrow(gene_attr) == 1){
    gene_sequence <- human_sequences[[as.numeric(gene_attr$index)]]
  }else{
    match_lengths <- sapply(as.numeric(gene_attr$index), function(i){
      summary(human_sequences[[as.numeric(gene_attr$index)]])$length
    })
    gene_sequence <- human_sequences[[as.numeric(gene_attr$index)[which.max(match_lengths)]]]
    warning(gene, " is degenerate with ", nrow(gene_attr), " sequence matches (taking the longest)")
  }
  gene_sequence <- DNAString(paste(getSequence(gene_sequence), collapse = ""))
  
  # Load sanger file and find regions of homo/heterozygocity
  ABIfile <- readsangerseq(paste0("ABIfiles/", an_ABIfile))
  hetcalls <- makeBaseCalls(ABIfile, ratio = 0.33)
  hetCode <- convert2IUPAC(hetcalls, returnType = "code") # convert sanger sequence to indicate degenerate bases
  
  # find a stretch of homozygocity matching WT DNA sequence strongly
  homozMat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -5, baseOnly = F)
  homozMatch <- pairwiseAlignment(primarySeq(hetcalls), gene_sequence, substitutionMatrix = homozMat, type = "local", gapOpening=-15, gapExtension=-15)
  # writePairwiseAlignments(homozMatch)
  
  if(homozMatch@pattern@range@width < 50){
    warning("longest match is only ", homozMatch@pattern@range@width, " nucleotides long")
  }
  
  # align the sanger sequence and WT sequence at the point of matching identity
  # trim sanger and gene sequence to the end of alignment
  hetCode <- DNAString(hetCode, start = homozMatch@pattern@range@start + homozMatch@pattern@range@width)
  WTaligned <- DNAString(gene_sequence, start = homozMatch@subject@range@start + homozMatch@subject@range@width)
  
  hetMat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -5, baseOnly = F)
  for(hetNuc in rownames(IUPAC)[rowSums(IUPAC) == 2]){
    hetMat[,hetNuc][names(IUPAC[rownames(IUPAC) == hetNuc,])[IUPAC[rownames(IUPAC) == hetNuc,] == 1]] <- 2
    hetMat[hetNuc,][names(IUPAC[rownames(IUPAC) == hetNuc,])[IUPAC[rownames(IUPAC) == hetNuc,] == 1]] <- 2
  }
  
  # Reconstruct one WT copy by finding indels which are necessary to achieve a heterozygous genotype that is consistent with the WT sequence
  # Map the WT sequence to the full sanger read and infer the complementary base of heterozygous regions.
  # Map this offset sequence back to the gene using a standard substitution matrix
  
  GlobalDel <- pairwiseAlignment(gene_sequence, DNAString(hetCode), type = "local-global", substitutionMatrix = hetMat, gapExtension = -8, gapOpening = -50)
  # call deletions on copy A
  GlobalDel@pattern@indel
  GlobalDel@subject@indel
  inferredMap
  
  compMat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -8, baseOnly = F) 
  
  inferredComplement <- haplotypeSubtract(consensusString(GlobalDel), hetCode)
  inferredMap <- pairwiseAlignment(gene_sequence, inferredComplement, type = "local-global", substitutionMatrix = compMat, gapExtension = -3, gapOpening = -50)
  inferredMap@pattern@indel
  inferredMap@subject@indel
  
  consensusString(inferredMap)
  
  hetMatch  <- pairwiseAlignment(hetCode, WTaligned, type = "local", substitutionMatrix = hetMat, gapExtension = -0.1)
  pairwiseAlignment(hetCode, WTaligned, type = "local", substitutionMatrix = hetMat)
  
  chromatogram(hetcalls, showcalls = "both")
  
  writePairwiseAlignments(inferredMap)
  # take the early stretch of homozygocity and align to the correct gene (local-global)
  # read forward to find stretches of WT that 
  
  # align @ homozygocity
  # identify interuption of heterozygocity breakpoint
  # use n+1 breakpoint to identify the nth indel
  
  # if gap at beginning of sanger - insertion
  # if gap at beginning of genomic - deletion
  
  
  
  
  
}

convert2IUPAC <- function(hetcalls, returnType = "code"){
  # returnType is:
  # matrix: either return a matrix with 4 columns indicating which bases are in agreement with the sanger call
  # code: IUPAC 1-letter summary of degenerate matches
  
  baseCalls <- data.frame(P = strsplit(primarySeq(hetcalls, string = T), split = '')[[1]], S = strsplit(secondarySeq(hetcalls, string = T), split = '')[[1]])

  if(returnType == "code"){
    
    IUPACcall <- sapply(1:nrow(baseCalls), function(i){
      totalCall <- IUPAC[rownames(IUPAC) == baseCalls$P[i],] + IUPAC[rownames(IUPAC) == baseCalls$S[i],]
      totalCall[totalCall > 1] <- 1
      rownames(IUPAC)[apply(IUPAC, 1, function(x){all(x == totalCall)})]
    })
    
    return(paste(IUPACcall, collapse = ""))
    
  }else if(returnType == "matrix"){
    IUPACcall <- t(sapply(1:nrow(baseCalls), function(i){
      totalCall <- IUPAC[rownames(IUPAC) == baseCalls$P[i],] + IUPAC[rownames(IUPAC) == baseCalls$S[i],]
      totalCall[totalCall > 1] <- 1
      totalCall
    }))
    
    return(IUPACcall)
    
  }else{
  stop("invalid arguement for returnType")
  }
  
}
  
haplotypeSubtract <- function(WT, Het){
  WT = strsplit(WT, split = "")[[1]]
  Het = strsplit(Het, split = "")[[1]]
  
  haplo <- mapply(aWT = WT, aHet = Het, function(aWT, aHet){
    if(aWT == "-"){aHet}else{
      if(all(IUPAC[rownames(IUPAC) == aWT,] == IUPAC[rownames(IUPAC) == aHet,])){
        rownames(IUPAC)[rownames(IUPAC) == aWT]
      }else{
        base_diff <- IUPAC[rownames(IUPAC) == aHet,] - IUPAC[rownames(IUPAC) == aWT,]
        base_diff[base_diff < 0] <- 0
        rownames(IUPAC)[apply(IUPAC, 1, function(z){
          all(z == base_diff)
        })]
      }
    }
  })
  
  return(DNAString(paste(haplo, collapse = "")))
}





http://imperialis.inhs.illinois.edu/dmitriev/indel.asp