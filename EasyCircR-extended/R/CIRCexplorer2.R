# ============================================================
# CIRCexplorer2 Pipeline Wrapper
# ============================================================

read_samplefile <- function(samples_file) {
  
  if (!file.exists(samples_file)) {
    stop("Sample file not found: ", samples_file)
  }
  
  samples <- read.table(
    samples_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_cols <- c("FileName1", "FileName2", "SampleName")
  
  if (!all(required_cols %in% colnames(samples))) {
    stop("Sample file must contain columns: FileName1, FileName2, SampleName")
  }
  
  return(samples)
}



run_circexplorer2 <- function(samples_file,
                              genome_fa,
                              ref_txt,
                              trim_reads_length = NULL,
                              force = FALSE,
                              threads = parallel::detectCores(),
                              remove_temporary_files = TRUE) {
  
  samples <- read_samplefile(samples_file)
  
  base_outdir <- "EasyCirc/circRNA"
  trimdir     <- file.path(base_outdir, "trimmed_fastq")
  circ_dir    <- file.path(base_outdir, "CIRCexplorer2")
  
  dir.create(base_outdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(circ_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------------------------------------------------
  # Check BWA index
  # ------------------------------------------------------------
  if (!file.exists(paste0(genome_fa, ".bwt"))) {
    cat("BWA index not found. Building index...\n")
    system(paste("bwa index", genome_fa))
  } else {
    cat("BWA index found. Skipping indexing.\n")
  }
  
  # ------------------------------------------------------------
  # Loop over samples
  # ------------------------------------------------------------
  for (i in 1:nrow(samples)) {
    
    samplename <- samples$SampleName[i]
    fq1 <- samples$FileName1[i]
    fq2 <- samples$FileName2[i]
    
    sample_outdir <- file.path(circ_dir, samplename)
    
    if (force && dir.exists(sample_outdir)) {
      unlink(sample_outdir, recursive = TRUE)
    }
    
    # ------------------------------------------------------------
    # Optional trimming
    # ------------------------------------------------------------
    if (!is.null(trim_reads_length)) {
      
      .trim_reads(
        fq1 = fq1,
        fq2 = fq2,
        trim_reads_length = trim_reads_length,
        outdir = trimdir,
        samplename = samplename,
        force = force,
        n_core = threads
      )
      
      fq1 <- file.path(trimdir, paste0(samplename, "_R1.trimmed.fq.gz"))
      fq2 <- file.path(trimdir, paste0(samplename, "_R2.trimmed.fq.gz"))
    }
    
    # ------------------------------------------------------------
    # Run CIRCexplorer2
    # ------------------------------------------------------------
    if (!dir.exists(sample_outdir)) {
      
      cat("------------------------------------------------------------\n")
      cat("Running CIRCexplorer2:", samplename, "\n")
      
      .run_circexplorer2(
        fq1 = fq1,
        fq2 = fq2,
        genome_fa = genome_fa,
        ref_txt = ref_txt,
        outdir = circ_dir,
        samplename = samplename,
        threads = threads
      )
      
    } else {
      cat(samplename, "already processed ... skipping\n")
    }
  }
  
  # ------------------------------------------------------------
  # Clean intermediate files
  # ------------------------------------------------------------
  if (remove_temporary_files) {
    clean_circexplorer2_folder(samples_file)
  }
  
  cat("\n✅ All samples processed successfully!\n")
}



# ============================================================
# Internal Function: Run CIRCexplorer2 (Single Sample)
# ============================================================

.run_circexplorer2 <- function(fq1, fq2,
                               genome_fa,
                               ref_txt,
                               outdir,
                               samplename,
                               threads = 10) {
  
  cat("------------------------------------------------------------\n")
  cat("CIRCexplorer2:", samplename, "\n")
  cat("------------------------------------------------------------\n")
  
  sample_dir <- file.path(outdir, samplename)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  sam_file  <- file.path(sample_dir, paste0(samplename, ".sam"))
  parse_log <- file.path(sample_dir, "parse.log")
  anno_log  <- file.path(sample_dir, "annotate.log")
  
  # ------------------------------------------------------------
  # 1. BWA Alignment
  # ------------------------------------------------------------
  bwa_cmd <- paste(
    "bash -c 'bwa mem -t", threads,
    genome_fa,
    "<(zcat", fq1, ")",
    "<(zcat", fq2, ")",
    ">", sam_file, "'"
  )
  
  cat("Running BWA...\n")
  system(bwa_cmd)
  
  # ------------------------------------------------------------
  # 2. Parse
  # ------------------------------------------------------------
  parse_cmd <- paste(
    "cd", sample_dir, "&&",
    "CIRCexplorer2 parse -t BWA",
    basename(sam_file),
    ">", basename(parse_log)
  )
  
  cat("Running parse...\n")
  system(parse_cmd)
  
  # ------------------------------------------------------------
  # 3. Annotate
  # ------------------------------------------------------------
  annotate_cmd <- paste(
    "cd", sample_dir, "&&",
    "CIRCexplorer2 annotate",
    "-r", ref_txt,
    "-g", genome_fa,
    "-b back_spliced_junction.bed",
    "-o circularRNA_known.txt",
    ">", basename(anno_log)
  )
  
  cat("Running annotate...\n")
  system(annotate_cmd)
  
  cat("CIRCexplorer2 finished:", samplename, "\n")
}



# ============================================================
# Trim Reads (Optional)
# ============================================================

.trim_reads <- function(fq1, fq2,
                        trim_reads_length,
                        outdir,
                        samplename,
                        force = FALSE,
                        n_core = parallel::detectCores()) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  fq1_paired   <- file.path(outdir, paste0(samplename,"_R1.trimmed.fq.gz"))
  fq1_unpaired <- file.path(outdir, paste0(samplename,"_R1un.trimmed.fq.gz"))
  fq2_paired   <- file.path(outdir, paste0(samplename,"_R2.trimmed.fq.gz"))
  fq2_unpaired <- file.path(outdir, paste0(samplename,"_R2un.trimmed.fq.gz"))
  
  if (!force && file.exists(fq1_paired) && file.exists(fq2_paired)) {
    cat("Already trimmed ... skipping\n")
  } else {
    
    trimmomatic_jar <- system.file("java","trimmomatic-0.39.jar", package="EasyCircR")
    cutfilter <- paste0("CROP:",trim_reads_length," MINLEN:",trim_reads_length)
    
    command <- paste("java -jar", trimmomatic_jar,
                     "PE -phred33 -threads", n_core,
                     fq1, fq2,
                     fq1_paired, fq1_unpaired,
                     fq2_paired, fq2_unpaired,
                     cutfilter)
    
    cat("Running Trimmomatic...\n")
    system(command)
  }
}



# ============================================================
# Clean Intermediate Files
# ============================================================

clean_circexplorer2_folder <- function(samples_file) {
  
  samples <- read_samplefile(samples_file)
  base_folder <- "EasyCirc/circRNA/CIRCexplorer2"
  
  for (i in 1:nrow(samples)) {
    
    samplename <- samples$SampleName[i]
    sample_path <- file.path(base_folder, samplename)
    
    if (dir.exists(sample_path)) {
      
      cat("Cleaning intermediate files in:", sample_path, "\n")
      
      unlink(file.path(sample_path, paste0(samplename, ".sam")))
      unlink(file.path(sample_path, "parse.log"))
      unlink(file.path(sample_path, "annotate.log"))
    }
  }
}