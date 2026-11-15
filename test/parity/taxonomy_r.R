# R half of the taxonomy parity check.
#   Rscript taxonomy_r.R <outdir> <train_set.fa.gz> <species.fa.gz> <queries.tsv (column "sequence")>
# Then: julia --project=. test/parity/taxonomy_julia.jl <outdir> <train_set> <species>
suppressMessages(library(dada2))
a <- commandArgs(TRUE)
out <- a[1]; ref <- a[2]; spref <- a[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f <- function(x) file.path(out, x)
w <- function(m, name) write.table(m, f(name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(Sys.getlocale("LC_COLLATE"), f("locale.txt"))

seqs <- read.delim(a[4])$sequence
s22 <- seqs[22]
seqs <- c(seqs, rc(seqs[1:20]), substr(seqs[21], 1, 45),
          paste0(substr(s22, 1, 100), "NNNN", substr(s22, 105, nchar(s22))))
writeLines(seqs, f("queries.txt"))
ext <- system.file("extdata", package = "dada2")
writeLines(c(file.path(ext, "example_train_set.fa.gz"), file.path(ext, "example_species_assignment.fa.gz")),
           f("extdata.txt"))

runs <- list(default = list(ref = ref, tryRC = FALSE, seed = 100, minBoot = 50),
             rc      = list(ref = ref, tryRC = TRUE,  seed = 100, minBoot = 50),
             boot80  = list(ref = ref, tryRC = FALSE, seed = 7,   minBoot = 80),
             example = list(ref = file.path(ext, "example_train_set.fa.gz"), tryRC = FALSE, seed = 1, minBoot = 50))
tt <- NULL
for (nm in names(runs)) {
  r <- runs[[nm]]
  set.seed(r$seed)
  t0 <- Sys.time()
  t <- assignTaxonomy(seqs, r$ref, minBoot = r$minBoot, tryRC = r$tryRC, outputBootstraps = TRUE, multithread = FALSE)
  if (nm == "default") writeLines(sprintf("%.1f", as.numeric(Sys.time() - t0, units = "secs")), f("time_default.txt"))
  w(t$tax, paste0("tax_", nm, ".tsv")); w(t$boot, paste0("boot_", nm, ".tsv"))
  if (nm == "default") tt <- t$tax
}

acgt <- dada2:::C_isACGT(seqs)
writeLines(as.character(as.integer(acgt)), f("acgt.txt"))
q <- seqs[acgt]
w(assignSpecies(q, spref), "sp_single.tsv")
w(assignSpecies(q, spref, allowMultiple = TRUE, tryRC = TRUE), "sp_multi_rc.tsv")
w(assignSpecies(q, spref, allowMultiple = 3), "sp_multi3.tsv")
w(assignSpecies(q, file.path(ext, "example_species_assignment.fa.gz"), allowMultiple = TRUE), "sp_example.tsv")
w(addSpecies(tt[acgt, ], spref), "addsp.tsv")
w(addSpecies(tt[acgt, ], spref, allowMultiple = TRUE, tryRC = TRUE), "addsp_multi_rc.tsv")
cat("R done\n")
