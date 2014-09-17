setwd("~/Desktop/Misc-Code/CRISPR-Caller")

options(stringsAsFactors = F)

# Analyze ABI sequencing file to determine whether a deletion

# install bioconductor packages
if(!all(c("sangerseqR", "seqinr") %in% rownames(installed.packages()))){
  source("http://bioconductor.org/biocLite.R")
  biocLite("sangerseqR")
  biocLite("seqinr")
  biocLite("org.Hs.eg.db")
  biocLite("biomaRt")

}

library(sangerseqR)
library(seqinr)
library(org.Hs.eg.db)

library(reshape2)
library(zoo)
library(data.table)

###### Functions ########

convert2IUPAC <- function(hetcalls, returnType = "code"){
  
  # Convert sangerseq objects into IUPAC DNA code
  
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
  
  # Take a Heterozygous genotype which is presumably a mixture of a WT haplotype (with indels) and an unknown haplotype
  # base-by-base subtract the WT genotype from the heterozygote to leave a predicted complementary haplotype
  
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

###### Prep ######

# Human coding sequences taken from Ensembl CDS dump http://useast.ensembl.org/info/data/ftp/index.html
# By default saved as .tsv, just altering extension works because it is formatted as a fasta

if(!file.exists("humanTranscript200bpFlanking.Rdata")){
  
  # If not previously performed then turn the fasta file into a list and indexing table so that future loading is faster
  # This is still slow because of the large object which must be loaded
  
  human_sequences <- read.fasta('humanTranscript200bpFlanking.fasta')
  seq_attrs <- t(sapply(1:length(human_sequences), function(x){
    seq_attr <- strsplit(getAnnot(human_sequences[[x]]), split = '[>\\|]')[[1]][-1]
    data.frame(index = x, ensemblGene = seq_attr[1], common = seq_attr[3])
  }, simplify = T))
  seq_attrs <- as.data.frame(apply(seq_attrs, c(1,2), function(x){x[[1]]}))
  save(human_sequences, seq_attrs, file = "humanTranscript200bpFlanking.Rdata")
  
}else{
  load("humanTranscript200bpFlanking.Rdata")
}

library(biomaRt)
ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
getGene("ENSG00000000005", type = "entrezgene", mart = ensembl)

listAttributes(ensembl,page="sequences")
genes = getBM(
  attributes=c("gene_exon_intron","ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
  filters=c("chromosome_name","start","end"),
  values=list("1",231500000,231800000), mart=ensembl
)

# Format IUPAC map from Biostrings

IUPAC <- melt(lapply(IUPAC_CODE_MAP, function(x){strsplit(x, split = '')}))
colnames(IUPAC) <- c("Nucleotide", "Level", "IUPAC_base")
IUPAC <- dcast(IUPAC, IUPAC_base ~ Nucleotide, value.var = "Level", fill = 0)
rownames(IUPAC) <- IUPAC$IUPAC_base; IUPAC <- IUPAC[,-1]

# substitution cost matrices

# HetMat = matrix which allows alignment of heterozygous sequences to reference - e.g. A-A = M-A = M-C
hetMat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -5, baseOnly = F)
for(hetNuc in rownames(IUPAC)[rowSums(IUPAC) == 2]){
  hetMat[,hetNuc][names(IUPAC[rownames(IUPAC) == hetNuc,])[IUPAC[rownames(IUPAC) == hetNuc,] == 1]] <- 2
  hetMat[hetNuc,][names(IUPAC[rownames(IUPAC) == hetNuc,])[IUPAC[rownames(IUPAC) == hetNuc,] == 1]] <- 2
}

# Standard cost of match and mismatch
compMat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -8, baseOnly = F) 

test_gene <- paste(c(rep("C", times = 50), rep("T", times = 50), rep("C", times = 50), rep("T", times = 50), rep("C", times = 50)), collapse = "")
test_query <- paste(c(rep("C", times = 75), rep("A", times = 50), rep("C", times = 75)), collapse = "")
test_align <- pairwiseAlignment(test_gene, test_query, type = "global")



test_align <- pairwiseAlignment(test_query, test_gene, type = "global")
consensusString(test_align)
writePairwiseAlignments(test_align)



##### Iterate through ABIfiles and call indels #####

ABIfiles <- list.files("ABIfiles")
an_ABIfile <- ABIfiles[3]
for(an_ABIfile in ABIfiles){
  
  ### Find the sequence of the wild-type gene ###
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
  
  ### Load sanger file and find regions of homo/heterozygocity ###
  ABIfile <- readsangerseq(paste0("ABIfiles/", an_ABIfile))
  
  hetcalls <- makeBaseCalls(ABIfile, ratio = 0.33)
  hetCode <- convert2IUPAC(hetcalls, returnType = "code") # convert sanger sequence to indicate degenerate bases
  
  # Reconstruct one WT copy by finding indels which are necessary to achieve a heterozygous genotype that is consistent with the WT sequence
  # Map the WT sequence to the full sanger read and infer the complementary base of heterozygous regions.
  # Map this offset sequence back to the gene using a standard substitution matrix
  
  # Determine quality of sanger calls
  peak_position_matrix <- melt(hetcalls@peakPosMatrix)
  colnames(peak_position_matrix) <- c("Position", "Base", "Trace_position")
  peak_position_matrix <- data.table(peak_position_matrix)
  peak_position_matrix$I <- sapply(1:nrow(peak_position_matrix), function(i){
    if(is.na(peak_position_matrix$Trace_position[i])){
      peak_region <- peak_position_matrix$Trace_position[peak_position_matrix$Position == peak_position_matrix$Position[i]]
      peak_region <- peak_region[!is.na(peak_region)]
      max(hetcalls@traceMatrix[min(peak_region):max(peak_region),peak_position_matrix$Base[i]])
    }else{
      hetcalls@traceMatrix[peak_position_matrix$Trace_position[i], peak_position_matrix$Base[i]]
    }
    
  })
  fluoroCounts <- acast(peak_position_matrix, Position ~ Base, value.var = "I")
  # As a quality score - Sum the two greatest counts and divide by the sum of the lower two counts
  fluoroCounts <- fluoroCounts + 1
  sQS <- apply(fluoroCounts, 1, function(i){
    sum(sort(i)[3:4]) / sum(sort(i)[1:2])
  })
  plot(log2(sQS))
  
  # Call the first haplotype
  
  GlobalDel <- pairwiseAlignment(gene_sequence, DNAString(hetCode), type = "local-global", substitutionMatrix = hetMat, gapExtension = -8, gapOpening = -50)
  # call deletions on copy A
  GlobalDel@pattern@indel[[1]]
  GlobalDel@subject@indel[[1]]
  
  #alleleAseq <- pairwiseAlignment(DNAString(consensusString(GlobalDel)), gene_sequence, type = "local", gapExtension = -3, gapOpening = -50)
  #writePairwiseAlignments(alleleAseq)
  
  # Call the second allele
  if(nchar(consensusString(GlobalDel)) != nchar(hetCode)){
    stop("Sequence is misaligned to sanger")
  }
  
  inferredComplement <- haplotypeSubtract(consensusString(GlobalDel), hetCode)
  inferredMap <- pairwiseAlignment(gene_sequence, inferredComplement, type = "local-global", substitutionMatrix = compMat, gapExtension = -8, gapOpening = -50)
  
  inferredMap@pattern@indel[[1]]
  inferredMap@subject@indel[[1]]
  
  #alleleBseq <- pairwiseAlignment(gene_sequence, DNAString(consensusString(inferredMap)), type = "local", gapExtension = -3, gapOpening = -50)
  #writePairwiseAlignments(alleleBseq)
  
  stringSummary <- DNAStringSet(c(consensusString(GlobalDel), consensusString(inferredMap)))
  consensusToGene <- pairwiseAlignment(stringSummary, gene_sequence, type = "local", gapExtension = -3, gapOpening = -50)
  writePairwiseAlignments(consensusToGene)
  
  compareStrings(tmp2)
  
  # Call both mutant alleles and then align each to WT
  # Determine quality of sanger sequence 
  # Flag if one allele is WT
  an_Alignment <- GlobalDel
  an_Alignment <- inferredMap
  QS = sQS
  
  
  
  filter_indels <- function(an_Alignment, QS){
    # insertion and deletion are reversed because functions are w.r.t. gene sequence
    allele_deletions <- insertion(an_Alignment)[[1]]
    allele_insertions <- deletion(an_Alignment)[[1]]
    
    allele_called_indels <- NULL
    if(length(allele_deletions@start) != 0){
      allele_called_indels <- rbind(allele_called_indels, data.frame(type = "deletion", allele_deletions))
    }
    if(length(allele_insertions@start) != 0){
      allele_called_indels <- rbind(allele_called_indels, data.frame(type = "insertion", allele_insertions))
    }
    allele_called_indels <- data.table(allele_called_indels)
    
    allele_called_indels[flankingQS := max(1, start - 2):min(nchar(an_Alignment), end + 2)]
    
    
  }
  
  
  
  
  
  hetMatch  <- pairwiseAlignment(hetCode, WTaligned, type = "local", substitutionMatrix = hetMat, gapExtension = -0.1)
  pairwiseAlignment(hetCode, WTaligned, type = "local", substitutionMatrix = hetMat)
  
  chromatogram(hetcalls, showcalls = "both")
  
  
  
}

