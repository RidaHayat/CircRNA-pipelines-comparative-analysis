
library(dplyr)

generate_circexplorer2_output_files <- function(circexplorer2_outdir,
                                                outdir,
                                                chrMT = TRUE,
                                                samples = NULL,
                                                save_rds = TRUE,
                                                save_csv = TRUE) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message("Created output folder: ", outdir)
  }
  
  # List sample directories
  sample_dirs <- list.dirs(circexplorer2_outdir, full.names = TRUE, recursive = FALSE)
  if (!is.null(samples)) {
    sample_dirs <- sample_dirs[sample_dirs %in% file.path(circexplorer2_outdir, samples)]
  }
  
  circ_df_list <- list()
  circ_mtx_list <- list()
  
  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    message("\nProcessing sample: ", sample_name)
    
    circ_bsj_file  <- file.path(sample_dir, "back_spliced_junction.bed")
    circ_anno_file <- file.path(sample_dir, "circularRNA_known.txt")
    
    if (!file.exists(circ_bsj_file))  stop("back_spliced_junction.bed not found for ", sample_name)
    if (!file.exists(circ_anno_file)) stop("circularRNA_known.txt not found for ", sample_name)
    
    # ----------------------------
    # 1) Read BSJ bed
    # ----------------------------
    bsj <- read.table(circ_bsj_file,
                      header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE,
                      quote = "", comment.char = "")
    
    if (ncol(bsj) < 3) stop("Unexpected BED format (<3 columns) in: ", circ_bsj_file)
    
    # standardize columns
    chr   <- bsj[[1]]
    start <- as.integer(bsj[[2]])
    end   <- as.integer(bsj[[3]])
    
    # score/count: BED col5 if exists and numeric, else 1
    count <- rep(1L, length(chr))
    if (ncol(bsj) >= 5) {
      suppressWarnings(tmp <- as.numeric(bsj[[5]]))
      if (!all(is.na(tmp))) count <- ifelse(is.na(tmp), 1, tmp)
    }
    
    # strand: BED col6 if exists else "*"
    strand <- rep("*", length(chr))
    if (ncol(bsj) >= 6) strand <- bsj[[6]]
    
    bsj_core <- paste0(chr, ":", start, "|", end)
    bsj_id   <- paste0(bsj_core, ":", strand)
    
    bsj_df <- data.frame(
      chr = chr,
      start = start,
      end = end,
      strand = strand,
      bsj_core = bsj_core,
      bsj_id = bsj_id,
      length = (end - start + 1),
      stringsAsFactors = FALSE
    )
    
    # ----------------------------
    # 2) Read annotation (circularRNA_known.txt)
    # ----------------------------
    anno <- read.table(circ_anno_file,
                       header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE,
                       quote = "", comment.char = "",
                       check.names = FALSE)
    
    # Renaming the columns to match our expected format
    colnames(anno) <- c("chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
                        "itemRgb", "exonCount", "exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    
    # Remove unnecessary columns
    anno <- anno[, c("chrom", "start", "end", "name", "strand", "readNumber", "geneName")]
    
    # Add BSJ identifier
    anno$bsj_core <- paste0(anno$chrom, ":", anno$start, "|", anno$end)
    
    # Merge BSJ info with annotation data
    circ_df <- merge(bsj_df, anno, by = "bsj_core", all.x = TRUE)
    
    # Add the sample name to circ_df
    circ_df$sample <- sample_name
    
    # ----------------------------
    # 3) Build per-sample count matrix
    # ----------------------------
    circ_mtx <- data.frame(bsj_id = bsj_id, count = count, stringsAsFactors = FALSE) %>%
      group_by(bsj_id) %>%
      summarise(!!sample_name := sum(count, na.rm = TRUE), .groups = "drop")
    
    circ_df_list[[sample_name]]  <- circ_df
    circ_mtx_list[[sample_name]] <- circ_mtx
  }
  
  # ----------------------------
  # 4) Combine circ_df and circ_mtx across samples
  # ----------------------------
  circ_df_all <- bind_rows(circ_df_list)
  
  circ_mtx_all <- circ_mtx_list[[1]]
  if (length(circ_mtx_list) > 1) {
    for (i in 2:length(circ_mtx_list)) {
      circ_mtx_all <- full_join(circ_mtx_all, circ_mtx_list[[i]], by = "bsj_id")
    }
  }
  circ_mtx_all[is.na(circ_mtx_all)] <- 0
  
  # set rownames like CIRI-full style
  circ_mtx_mat <- as.data.frame(circ_mtx_all)
  rownames(circ_mtx_mat) <- circ_mtx_mat$bsj_id
  circ_mtx_mat$bsj_id <- NULL
  
  # ----------------------------
  # 5) Optionally remove mitochondrial circRNAs
  # ----------------------------
  if (!chrMT) {
    mt_vals <- c("MT", "chrM", "M", "chrMT")
    circ_df_all <- circ_df_all %>% filter(!(chrom %in% mt_vals))
    circ_mtx_mat <- circ_mtx_mat[rownames(circ_mtx_mat) %in% circ_df_all$bsj_id, , drop = FALSE]
  }
  
  # ----------------------------
  # 6) Save the final results
  # ----------------------------
  if (save_rds) {
    saveRDS(circ_df_all, file = file.path(outdir, "circ_df_all.rds"))
    saveRDS(circ_mtx_mat, file = file.path(outdir, "circ_mtx_all.rds"))
  }
  if (save_csv) {
    write.csv(circ_df_all, file = file.path(outdir, "circ_df_all.csv"), row.names = FALSE)
    write.csv(circ_mtx_mat, file = file.path(outdir, "circ_mtx_all.csv"), row.names = TRUE)
  }
  
  return(list(circ_df = circ_df_all, circ_mtx = circ_mtx_mat))
}
