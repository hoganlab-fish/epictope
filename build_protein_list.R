#!/usr/bin/env Rscript
# Rebuild zebrafish_proteins.csv (the app's dropdown list) from UniProt's
# zebrafish reference proteome: one entry per gene (the "one sequence per
# gene" set, reviewed Swiss-Prot and unreviewed TrEMBL entries alike). That's
# the same gene set as AlphaFold DB's zebrafish proteome archive, so most
# listed proteins can have a structure on disk.
# usage: Rscript build_protein_list.R [output.csv]   (run from the repo root)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "zebrafish_proteins.csv"
proteome <- "UP000000437_7955"
options(timeout = max(600, getOption("timeout")))

# Full UniProt flat-file records for the one-per-gene set (~43 MB download)
dat_file <- tempfile(fileext = ".dat.gz")
download.file(paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/",
                     "reference_proteomes/Eukaryota/", sub("_.*", "", proteome), "/", proteome, ".dat.gz"),
              dat_file, mode = "wb", quiet = TRUE)
dat <- readLines(gzfile(dat_file))

# Each record starts with an ID line; number records so fields from other
# lines can be matched back to them.
is_id <- startsWith(dat, "ID   ")
record <- cumsum(is_id)
n <- sum(is_id)

# First value of a field per record ("" where the record has none)
first_per_record <- function(keep, pattern) {
  vals <- sub(pattern, "\\1", dat[keep])
  first <- !duplicated(record[keep])
  res <- character(n)
  res[record[keep][first]] <- trimws(vals[first])
  res
}

# Protein names in UniProt's display style, e.g. "Name (EC 1.2.3.4) (Alt name)
# (Short name)": every Full=/Short=/EC= value on the record's DE lines, up to
# any Contains:/Includes: section (names of cleaved chains/domains).
de <- which(startsWith(dat, "DE   "))
de_rec <- record[de]
de_txt <- dat[de]
in_main <- ave(grepl("^DE   (Contains|Includes):", de_txt), de_rec, FUN = cumsum) == 0
parts <- regmatches(de_txt, regexec("(Full|Short|EC)=([^;{]+)", de_txt))
key <- vapply(parts, function(x) if (length(x)) x[2] else NA_character_, "")
val <- vapply(parts, function(x) if (length(x)) trimws(x[3]) else NA_character_, "")
ok <- !is.na(key) & in_main
val <- ifelse(key == "EC", paste("EC", val), val)
names_by_rec <- split(val[ok], de_rec[ok])
protein_names <- character(n)
protein_names[as.integer(names(names_by_rec))] <- vapply(names_by_rec, function(v)
  paste0(v[1], if (length(v) > 1) paste0(" (", v[-1], ")", collapse = "")), "")

id_lines <- dat[is_id]
res <- data.frame(
  # primary accession: first one on the first AC line
  uniprot_id   = first_per_record(startsWith(dat, "AC   "), "^AC   ([^;]+);.*"),
  gene_name    = first_per_record(startsWith(dat, "GN   Name="), "^GN   Name=([^;{ ]+).*"),
  protein_name = protein_names,
  length       = as.integer(sub(".*;\\s+([0-9]+) AA\\.$", "\\1", id_lines)),
  alphafold_id = first_per_record(startsWith(dat, "DR   AlphaFoldDB;"), "^DR   AlphaFoldDB; ([^;]+);.*")
)
res <- res[order(res$uniprot_id), ]
# quote only protein_name (it can contain commas), keeping the list's
# shell-friendly layout (plain header, unquoted IDs)
writeLines(paste(names(res), collapse = ","), out)
write.table(res, out, sep = ",", quote = 3, qmethod = "double", na = "",
            row.names = FALSE, col.names = FALSE, append = TRUE)

message(nrow(res), " proteins written to ", out, " (", sum(grepl(" Reviewed;", id_lines)),
        " reviewed; ", sum(nzchar(res$alphafold_id)), " with an AlphaFold DB entry; ",
        sum(!nzchar(res$gene_name)), " without a gene name)")
