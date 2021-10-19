#' @param seqs (Required). A character vector of the sequences to be assigned, or an object 
#' coercible by \code{\link{getUniques}}.
#'   
#' @param refFasta (Required). The path to the reference fasta file, or an 
#' R connection Can be compressed.
#' This reference fasta file should be formatted so that the id lines correspond to the
#' taxonomy (or classification) of the associated sequence, and each taxonomic level is 
#' separated by a semicolon. Eg.
#' 
#'  >Kingom;Phylum;Class;Order;Family;Genus;   
#'  ACGAATGTGAAGTAA......   
#' 
#' @param minBoot (Optional). Default 50. 
#' The minimum bootstrap confidence for assigning a taxonomic level.
#'   
#' @param tryRC (Optional). Default FALSE. 
#' If TRUE, the reverse-complement of each sequences will be used for classification if it is a better match to the reference
#' sequences than the forward sequence.
#'   
#' @param outputBootstraps (Optional). Default FALSE.
#'  If TRUE, bootstrap values will be retained in an integer matrix. A named list containing the assigned taxonomies (named "taxa") 
#'  and the bootstrap values (named "boot") will be returned. Minimum bootstrap confidence filtering still takes place,
#'  to see full taxonomy set minBoot=0
#'   
#' @param taxLevels (Optional). Default is c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species").
#' The taxonomic levels being assigned. Truncates if deeper levels not present in
#' training fasta.
#'   
#' @param multithread (Optional). Default is FALSE.
#'  If TRUE, multithreading is enabled and the number of available threads is automatically determined.   
#'  If an integer is provided, the number of threads to use is set by passing the argument on to
#'  \code{\link{setThreadOptions}}.
#'   
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, print status to standard output.
#'   
#' @return A character matrix of assigned taxonomies exceeding the minBoot level of
#'   bootstrapping confidence. Rows correspond to the provided sequences, columns to the
#'   taxonomic levels. NA indicates that the sequence was not consistently classified at
#'   that level at the minBoot threshhold.
#'   
#'   If outputBootstraps is TRUE, a named list containing the assigned taxonomies (named "taxa") 
#'   and the bootstrap values (named "boot") will be returned.
#' 
#' @export
#' 
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead id
#' 
#' @examples
#' seqs <- getSequences(system.file("extdata", "example_seqs.fa", package="dada2"))
#' training_fasta <- system.file("extdata", "example_train_set.fa.gz", package="dada2")
#' taxa <- assignTaxonomy(seqs, training_fasta)
#' taxa80 <- assignTaxonomy(seqs, training_fasta, minBoot=80, multithread=2)
#' 
assignTaxonomy <- function(seqs, refFasta, minBoot=50, tryRC=FALSE, outputBootstraps=FALSE,
                           taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                           multithread=FALSE, verbose=FALSE) {
  MIN_REF_LEN <- 20 # Enforced minimum length of reference seqs. Must be bigger than the kmer-size used (8).
  MIN_TAX_LEN <- 50 # Minimum length of input sequences to get a taxonomic assignment
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  if(min(nchar(seqs)) < MIN_TAX_LEN) {
    warning("Some sequences were shorter than ", MIN_TAX_LEN, " nts and will not receive a taxonomic classification.")
  }
  # Read in the reference fasta
  refsr <- readFasta(refFasta)
  lens <- width(sread(refsr))
  if(any(lens<MIN_REF_LEN)) {
    refsr <- refsr[lens>=MIN_REF_LEN]
    warning(paste0("Some reference sequences were too short (<", MIN_REF_LEN, "nts) and were excluded."))
  }
  refs <- as.character(sread(refsr))
  tax <- as.character(id(refsr))
  tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x)) # Remove leading/trailing whitespace
  # Sniff and parse UNITE fasta format
  UNITE <- FALSE
  if(all(grepl("FU\\|re[pf]s", tax[1:10]))) {
    UNITE <- TRUE
    cat("UNITE fungal taxonomic reference detected.\n")
    tax <- sapply(strsplit(tax, "\\|"), `[`, 5)
    tax <- gsub("[pcofg]__unidentified;", "_DADA2_UNSPECIFIED;", tax)
    tax <- gsub(";s__(\\w+)_", ";s__", tax)
    tax <- gsub(";s__sp$", ";_DADA2_UNSPECIFIED", tax)
  }
  # Crude format check
  if(!grepl(";", tax[[1]])) {
    if(length(unlist(strsplit(tax[[1]], "\\s")))==3) {
      stop("Incorrect reference file format for assignTaxonomy (this looks like a file formatted for assignSpecies).")
    } else {
      stop("Incorrect reference file format for assignTaxonomy.")
    }
  }
  # Parse the taxonomies from the id string
  tax.depth <- sapply(strsplit(tax, ";"), length)
  td <- max(tax.depth)
  for(i in seq(length(tax))) {
    if(tax.depth[[i]] < td) {
      for(j in seq(td - tax.depth[[i]])) {
        tax[[i]] <- paste0(tax[[i]], "_DADA2_UNSPECIFIED;")
      }
    }
  }
  # Create the integer maps from reference to type ("genus") and for each tax level
  genus.unq <- unique(tax)
  ref.to.genus <- match(tax, genus.unq)
  tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol=td, byrow=TRUE)
  tax.df <- as.data.frame(tax.mat)
  for(i in seq(ncol(tax.df))) {
    tax.df[,i] <- factor(tax.df[,i])
    tax.df[,i] <- as.integer(tax.df[,i])
  }
  tax.mat.int <- as.matrix(tax.df)
  ### Assign
  # Parse multithreading argument
  if(is.logical(multithread)) {
    if(multithread==TRUE) { RcppParallel::setThreadOptions(numThreads = "auto") }
    else { RcppParallel::setThreadOptions(numThreads = 1) }
  } else if(is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
  } else {
    warning("Invalid multithread parameter. Running as a single thread.")
    RcppParallel::setThreadOptions(numThreads = 1)
  }
  # Run C assignemnt code
  assignment <- C_assign_taxonomy2(seqs, rc(seqs), refs, ref.to.genus, tax.mat.int, tryRC, verbose)
  # Parse results and return tax consistent with minBoot
  bestHit <- genus.unq[assignment$tax]
  boots <- assignment$boot
  taxes <- strsplit(bestHit, ";")
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,]>=minBoot])
  # Convert to character matrix
  tax.out <- matrix(NA_character_, nrow=length(seqs), ncol=td)
  for(i in seq(length(seqs))) {
    if(length(taxes[[i]]) > 0) {
      tax.out[i,1:length(taxes[[i]])] <- taxes[[i]]
    }
  }
  rownames(tax.out) <- seqs
  colnames(tax.out) <- taxLevels[1:ncol(tax.out)]
  tax.out[tax.out=="_DADA2_UNSPECIFIED"] <- NA_character_
  if(outputBootstraps){
      # Convert boots to integer matrix
      boots.out <- matrix(boots, nrow=length(seqs), ncol=td)
      rownames(boots.out) <- seqs
      colnames(boots.out) <- taxLevels[1:ncol(boots.out)]
      list(tax=tax.out, boot=boots.out)
  } else {
    tax.out
  }
}

# Helper function for assignSpecies
mapHits <- function(x, refs, keep, sep="/") {
  hits <- refs[x]
  hits[grepl("Escherichia", hits, fixed=TRUE) | grepl("Shigella", hits, fixed=TRUE)] <- "Escherichia/Shigella"
  if(length(unique(hits))<=keep) {
    rval <- do.call(paste, c(as.list(sort(unique(hits))), sep=sep))
  } else { rval <- NA_character_ }
  if(length(rval)==0) rval <- NA_character_
  rval
}

# Match curated genus names to binomial genus names
# Handles Clostridium groups and split genera names
matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
  if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
  if((gen.tax==gen.binom) || 
     grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) || 
     grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
