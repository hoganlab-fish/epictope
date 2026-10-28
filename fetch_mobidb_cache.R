#!/usr/bin/env Rscript
# One-off: download MobiDB entries (including ANCHOR binding regions) for every
# protein in the protein list into data/mobidb/<accession>.json, so the IUPred2A
# fallback in app.R / single_score.R reads them from disk instead of the API.
# Safe to re-run: accessions already on disk are skipped. For the full
# zebrafish list (~26.6k proteins, 50 per request, 10-30 s apart) a first run
# takes about 3 hours.
# usage: Rscript fetch_mobidb_cache.R [proteins.csv]   (run from the repo root)

args <- commandArgs(trailingOnly = TRUE)
csv <- if (length(args) >= 1) args[1] else "zebrafish_proteins.csv"
out_dir <- file.path("data", "mobidb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

accs <- unique(trimws(read.csv(csv, stringsAsFactors = FALSE)$uniprot_id))
accs <- accs[!file.exists(file.path(out_dir, paste0(accs, ".json")))]
batches <- split(accs, ceiling(seq_along(accs) / 50))

for (i in seq_along(batches)) {
  # random 10-30 s gap between requests to go easy on the MobiDB server
  if (i > 1) Sys.sleep(runif(1, 10, 30))
  url <- paste0("https://mobidb.org/api/download?acc=",
                paste(batches[[i]], collapse = ","), "&format=json")
  # one JSON entry per line; accessions MobiDB doesn't have are left out
  lines <- readLines(url, warn = FALSE)
  lines <- lines[nzchar(lines)]
  for (line in lines) {
    writeLines(line, file.path(out_dir, paste0(jsonlite::fromJSON(line)$acc, ".json")))
  }
  message("batch ", i, "/", length(batches), ": ", length(lines), " of ",
          length(batches[[i]]), " found")
}

missing <- accs[!file.exists(file.path(out_dir, paste0(accs, ".json")))]
if (length(missing)) message("Not in MobiDB: ", paste(missing, collapse = ", "))
