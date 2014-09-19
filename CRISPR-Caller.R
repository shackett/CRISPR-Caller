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
library(biomaRt)

library(ggplot2)
library(reshape2)
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

test_genes <- c("Foxp2", "tp53", "ENSG00000139618", "made_up")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

query_gene <- test_genes[4]
find_gene_sequence(query_gene)

find_gene_sequence <- function(query_gene){
  
  require(biomaRt)
  
  if(length(query_gene) != 1){stop("Query a single gene")}
  
  genes = getBM(
    attributes=list("gene_exon_intron", "ensembl_gene_id","hgnc_symbol"),
    filters="hgnc_symbol",
    values=query_gene, mart=ensembl
  )
  if(nrow(genes) != 0){
    return(genes)
  }
  
  genes = getBM(
    attributes=list("gene_exon_intron", "ensembl_gene_id","hgnc_symbol"),
    filters="ensembl_gene_id",
    values=query_gene, mart=ensembl
  )
  if(nrow(genes) != 0){
    return(genes)
  }
  
  warning(query_gene, " does not match any human gene symbol or ensembl gene ID")
  return(NULL)
}

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
  # nucleotide numbering will correspond directly to position in QS vector
  
  indel_summary <- allele_called_indels[,list(localQC = median(QS[c(max(1, start - 2):min(an_Alignment@subject@range@width, end + 2))])),by = c("type", "start", "end", "width")]
  indel_summary[,verdict := ifelse(localQC < 30 | width < 5, "Misalignment", "True Indel")]
  
  if(!any(indel_summary$verdict == "True Indel")){
    warning("No probable indel located")
  }
  return(indel_summary)
}

###### Prep ######

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

#test_gene <- paste(c(rep("C", times = 50), rep("T", times = 50), rep("C", times = 50), rep("T", times = 50), rep("C", times = 50)), collapse = "")
#test_query <- paste(c(rep("C", times = 75), rep("A", times = 50), rep("C", times = 75)), collapse = "")
#test_align <- pairwiseAlignment(test_gene, test_query, type = "global")
#test_align <- pairwiseAlignment(test_query, test_gene, type = "global")
#consensusString(test_align)
#writePairwiseAlignments(test_align)



##### Iterate through ABIfiles and call indels #####

callAlleles <- function(inputFilePath, outputDirectoryPath, outputFile = NULL){
  
  
  
  
  
  
  
  
}

ABIfiles <- list.files("ABIfiles")
for(an_ABIfile in ABIfiles){
  
  ### Find the sequence of the wild-type gene ###
  gene <- strsplit(an_ABIfile, split = '[-_]')[[1]][2]
  gene_sequence_matches <- find_gene_sequence(gene)
  if(is.null(gene_sequence_matches)){
    warning("A gene sequence could not be located for", an_ABIfile)
    next
  }
  
  if(nrow(gene_sequence_matches) == 1){
    gene_sequence <- gene_sequence_matches$gene_exon_intron
  }else{
    gene_sequence <- gene_sequence_matches$gene_exon_intron[which.max(nchar(gene_sequence_matches$gene_exon_intron))]
    warning(gene, " is degenerate with ", nrow(gene_attr), " sequence matches (taking the longest)")
  }
  gene_sequence <- DNAString(gene_sequence)
  
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
  print(qplot(y = log10(sQS), x = 1:length(sQS)) + geom_hline(y = log10(30), size = 3, col = "RED"))
  
  # Call the first allele
  
  GlobalDel <- pairwiseAlignment(gene_sequence, DNAString(hetCode), type = "local-global", substitutionMatrix = hetMat, gapExtension = -8, gapOpening = -50)
  
  # Call the second allele
  if(nchar(consensusString(GlobalDel)) != nchar(hetCode)){
    stop("Sequence is misaligned to sanger")
  }
  
  inferredComplement <- haplotypeSubtract(consensusString(GlobalDel), hetCode)
  inferredMap <- pairwiseAlignment(gene_sequence, inferredComplement, type = "local-global", substitutionMatrix = compMat, gapExtension = -8, gapOpening = -50)
  
  stringSummary <- DNAStringSet(c(consensusString(GlobalDel), consensusString(inferredMap)))
  consensusToGene <- pairwiseAlignment(stringSummary, gene_sequence, type = "local", gapExtension = -3, gapOpening = -50)
  #writePairwiseAlignments(consensusToGene, file = "")
  
  cat(an_ABIfile)
  cat("\n------------------------------------------\n")
  print(filter_indels(GlobalDel, sQS), digits = 2)
  cat("------------------------------------------\n")
  print(filter_indels(inferredMap, sQS), digits = 2)
  cat('==========================================\n\n')
  
}

