# Load libraries once
library(dplyr)
library(magrittr)
library(tibble)

read_circtools_output <- function(circ_output_dir, samples = NULL, chrMT = TRUE, save_csv = TRUE) {
  
  # 1. Read CircCoordinates
  coords_file <- file.path(circ_output_dir, "CircCoordinates")
  if(!file.exists(coords_file)) stop("CircCoordinates file not found in ", circ_output_dir)
  
  circ_coords <- read.table(coords_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  message("CircCoordinates rows read: ", nrow(circ_coords))
  
  circ_coords$bsj_id <- paste0(circ_coords$Chr, ":", circ_coords$Start, "|", circ_coords$End, ":", circ_coords$Strand)
  circ_coords$length <- circ_coords$End - circ_coords$Start + 1
  
  circ_df <- circ_coords %>%
    dplyr::rename(
      chr = Chr,
      start = Start,
      end = End,
      strand = Strand,
      gene_id = Gene,
      isoform_cirexon = `Start-End Region`
    ) %>%
    dplyr::mutate(seq = NA, samples = NA)
  
  # 2. Read CircRNACount
  counts_file <- file.path(circ_output_dir, "CircRNACount")
  if(!file.exists(counts_file)) stop("CircRNACount file not found in ", circ_output_dir)
  
  circ_counts <- read.table(counts_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  message("CircRNACount rows read: ", nrow(circ_counts))
  
  if(!is.null(samples)) {
    colnames(circ_counts)[4:ncol(circ_counts)] <- samples
  }
  
  circ_counts$bsj_id <- paste0(circ_counts$Chr, ":", circ_counts$Start, "|", circ_counts$End, ":", "+")
  
  # Map strand from circ_df
  strand_map <- circ_df %>% select(bsj_id, strand)
  circ_counts <- circ_counts %>% left_join(strand_map, by = "bsj_id")
  circ_counts$bsj_id <- paste0(circ_counts$Chr, ":", circ_counts$Start, "|", circ_counts$End, ":", circ_counts$strand)
  
  # Build circ_mtx
  circ_mtx <- circ_counts %>%
    dplyr::select(bsj_id, dplyr::everything(), -Chr, -Start, -End, -strand) %>%
    tibble::column_to_rownames("bsj_id")
  
  # 3. Remove MT if requested
  if(!chrMT) {
    circ_df <- circ_df %>% filter(chr != "MT")
    circ_mtx <- circ_mtx[rownames(circ_mtx) %in% circ_df$bsj_id, , drop = FALSE]
  }
  
  # 4. Save CSV
  if(save_csv) {
    write.csv(circ_df, file.path(circ_output_dir, "circ_df.csv"), row.names = FALSE)
    write.csv(circ_mtx, file.path(circ_output_dir, "circ_mtx.csv"), row.names = TRUE)
    message("Saved circ_df.csv and circ_mtx.csv in ", circ_output_dir)
  }
  
  return(list(circ_df = circ_df, circ_mtx = circ_mtx))
}


