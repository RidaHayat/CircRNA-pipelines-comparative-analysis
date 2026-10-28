# =========================
# Helper: consistent sample naming
# =========================

get_sample_name <- function(fastq_file) {
  sample_name <- basename(fastq_file)
  # Remove typical FASTQ suffixes and extra temporary suffixes
  sample_name <- sub("(_R?1|_R?2|_1|_2|_1d|\\.rrnaols_detect)?\\.fastq(\\.gz)?$", "", sample_name)
  # Replace any weird characters with underscore
  sample_name <- gsub("[^A-Za-z0-9_-]", "_", sample_name)
  return(sample_name)
}
# =========================
# Read sample file
# =========================
read_samplefile <- function(samples_file) {
  if(!file.exists(samples_file)) stop("Samples file not found: ", samples_file)
  
  samples <- read.table(
    samples_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required_cols <- c("SampleName", "FileName1", "FileName2")
  if(!all(required_cols %in% colnames(samples))) stop("Sample file must contain columns: SampleName, FileName1, FileName2")
  
  return(samples)
}

# =========================
# Step1: Adapter removal & quality trimming with Flexbar
# =========================
run_flexbar_trim <- function(fastq_dir = "data/fastq",
                             output_dir = "flexbar",
                             test_sample = NULL,
                             flexbar_path = "flexbar",
                             min_read_length = 30) {
  
  if(!dir.exists(fastq_dir)) stop("FASTQ directory does not exist: ", fastq_dir)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  r1_files <- list.files(fastq_dir, pattern = "_R1\\.fastq(\\.gz)?$", full.names = TRUE)
  r2_files <- list.files(fastq_dir, pattern = "_R2\\.fastq(\\.gz)?$", full.names = TRUE)
  
  if(!is.null(test_sample)) {
    r1_files <- r1_files[grepl(test_sample, r1_files)]
    r2_files <- r2_files[grepl(test_sample, r2_files)]
    if(length(r1_files) == 0 || length(r2_files) == 0) stop("Test sample FASTQ not found: ", test_sample)
  }
  
  if(length(r1_files) != length(r2_files)) stop("Number of R1 and R2 files does not match!")
  r1_files <- sort(r1_files)
  r2_files <- sort(r2_files)
  
  for(i in seq_along(r1_files)) {
    fq1 <- r1_files[i]
    fq2 <- r2_files[i]
    sample_name <- get_sample_name(fq1)
    
    out_sample_dir <- file.path(output_dir, sample_name)
    if(!dir.exists(out_sample_dir)) dir.create(out_sample_dir, recursive = TRUE)
    
    out_prefix <- file.path(out_sample_dir, sample_name)
    
    cmd <- sprintf("%s -r %s -p %s -t %s --min-read-length %s",
                   flexbar_path, fq1, fq2, out_prefix, min_read_length)
    
    cat("\nTrimming sample:", sample_name, "\n")
    cat("Command:\n", cmd, "\n")
    system(cmd)
  }
  
  cat("\nAll trimming finished. Trimmed files are in:", output_dir, "\n")
  invisible(TRUE)
}

# =========================
# Step 2: rRNA removal with Bowtie2
# =========================
run_rrna_filter <- function(flexbar_dir = "flexbar",
                            rRNA_fasta = "data/human_rRNA.fasta",
                            output_dir = "rrna",
                            test_sample = NULL,
                            threads = 1,
                            bowtie2_path = "bowtie2") {
  
  if(!file.exists(bowtie2_path)) stop("Bowtie2 not found at: ", bowtie2_path)
  if(!dir.exists(flexbar_dir)) stop("Flexbar directory not found: ", flexbar_dir)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  index_dir <- "rrna_index"
  if(!dir.exists(index_dir)) dir.create(index_dir)
  index_prefix <- file.path(index_dir, "human_rRNA")
  
  if(!file.exists(paste0(index_prefix, ".1.bt2"))) {
    cat("Building Bowtie2 index from", rRNA_fasta, "...\n")
    system(sprintf("%s-build %s %s", bowtie2_path, rRNA_fasta, index_prefix))
  }
  
  r1_files <- list.files(flexbar_dir, pattern = "_1\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  r2_files <- list.files(flexbar_dir, pattern = "_2\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  
  if(!is.null(test_sample)) {
    r1_files <- r1_files[grepl(test_sample, r1_files)]
    r2_files <- r2_files[grepl(test_sample, r2_files)]
    if(length(r1_files) == 0 || length(r2_files) == 0) stop("Test sample not found: ", test_sample)
  }
  
  r1_files <- sort(r1_files)
  r2_files <- sort(r2_files)
  if(length(r1_files) != length(r2_files)) stop("Number of R1 and R2 files does not match!")
  
  for(i in seq_along(r1_files)) {
    sample_name <- get_sample_name(r1_files[i])
    cat("\nFiltering rRNA for sample:", sample_name, "\n")
    
    sample_out <- file.path(output_dir, sample_name)
    if(!dir.exists(sample_out)) dir.create(sample_out, recursive = TRUE)
    
    cmd <- sprintf("%s -x %s -1 %s -2 %s --very-sensitive -p %d --un-conc %s > /dev/null 2> %s.log",
                   bowtie2_path,
                   index_prefix,
                   r1_files[i],
                   r2_files[i],
                   threads,
                   file.path(sample_out, paste0(sample_name, "_%d.rrna.fastq")),
                   file.path(sample_out, sample_name))
    
    cat("Command:\n", cmd, "\n")
    system(cmd)
  }
  
  cat("\nrRNA removal finished. Filtered reads are in:", output_dir, "\n")
  invisible(TRUE)
}

# =========================
# Step 3: STAR genome index
# =========================
generate_star_index <- function(genome_fasta,
                                gtf_file,
                                output_dir,
                                threads = 10,
                                sjdb_overhang = 100,
                                star_exe,
                                force = FALSE) {
  
  if(dir.exists(output_dir) && !force) {
    message("STAR index already exists — skipping.")
    return(invisible(TRUE))
  }
  
  if(dir.exists(output_dir) && force) unlink(output_dir, recursive = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cmd <- sprintf("%s --runMode genomeGenerate --runThreadN %d --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang %d",
                 star_exe, threads, output_dir, genome_fasta, gtf_file, sjdb_overhang)
  
  message("Building STAR index...")
  system(cmd)
  
  invisible(TRUE)
}

# =========================
# Step 4: STAR joint mapping
# =========================
run_star_joint <- function(fastq_dir,
                           genome_dir,
                           output_dir = "star_quant",
                           gtf_file,
                           threads = 10,
                           star_exe,
                           test_sample = NULL) {
  
  if(!dir.exists(fastq_dir)) stop("FASTQ directory not found: ", fastq_dir)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  r1_files <- list.files(fastq_dir, pattern = "_1.*\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  r2_files <- list.files(fastq_dir, pattern = "_2.*\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  
  if(!is.null(test_sample)) {
    r1_files <- r1_files[grepl(test_sample, r1_files)]
    r2_files <- r2_files[grepl(test_sample, r2_files)]
    if(length(r1_files) == 0 || length(r2_files) == 0) stop("Test sample FASTQ not found: ", test_sample)
  }
  
  r1_files <- sort(r1_files)
  r2_files <- sort(r2_files)
  if(length(r1_files) != length(r2_files)) stop("Number of R1 and R2 files does not match!")
  
  for(i in seq_along(r1_files)) {
    sample_name <- get_sample_name(r1_files[i])
    sample_out <- file.path(output_dir, sample_name)
    if(!dir.exists(sample_out)) dir.create(sample_out, recursive = TRUE)
    
    out_prefix <- file.path(sample_out, paste0(sample_name, "_"))
    read_cmd <- ifelse(grepl("\\.gz$", r1_files[i]), " --readFilesCommand zcat \\\n", " \\\n")
    
    cmd <- paste0(
      star_exe, " \\\n",
      "  --runThreadN ", threads, " \\\n",
      "  --genomeDir ", genome_dir, " \\\n",
      "  --genomeLoad NoSharedMemory \\\n",
      "  --readFilesIn ", r1_files[i], " ", r2_files[i], read_cmd,
      "  --outFileNamePrefix ", out_prefix, " \\\n",
      "  --outReadsUnmapped Fastx \\\n",
      "  --outSAMattributes NH HI AS nM NM MD jM jI XS \\\n",
      "  --outSJfilterOverhangMin 15 15 15 15 \\\n",
      "  --outFilterMultimapNmax 20 \\\n",
      "  --outFilterScoreMin 1 \\\n",
      "  --outFilterMatchNminOverLread 0.7 \\\n",
      "  --outFilterMismatchNmax 999 \\\n",
      "  --outFilterMismatchNoverLmax 0.05 \\\n",
      "  --alignIntronMin 20 \\\n",
      "  --alignIntronMax 1000000 \\\n",
      "  --alignMatesGapMax 1000000 \\\n",
      "  --alignSJoverhangMin 15 \\\n",
      "  --alignSJDBoverhangMin 10 \\\n",
      "  --alignSoftClipAtReferenceEnds No \\\n",
      "  --chimSegmentMin 15 \\\n",
      "  --chimScoreMin 15 \\\n",
      "  --chimScoreSeparation 10 \\\n",
      "  --chimJunctionOverhangMin 15 \\\n",
      "  --sjdbGTFfile ", gtf_file, " \\\n",
      "  --quantMode GeneCounts \\\n",
      "  --twopassMode Basic \\\n",
      "  --chimOutType Junctions SeparateSAMold"
    )
    
    cat("\nMapping sample:", sample_name, "\n")
    cat("Command:\n", cmd, "\n")
    system(cmd)
  }
  
  cat("\nSTAR quantification finished. Output is in:", output_dir, "\n")
  invisible(TRUE)
}

# =========================
# Step 5: STAR mate-specific mapping
# =========================
run_star_mate_separate <- function(fastq_dir,
                                   genome_dir,
                                   output_dir = "star_quant",
                                   gtf_file,
                                   threads = 10,
                                   star_exe,
                                   test_sample = NULL) {
  
  if(!dir.exists(fastq_dir)) stop("FASTQ directory not found: ", fastq_dir)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  r1_files <- list.files(fastq_dir, pattern = "_1.*\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  r2_files <- list.files(fastq_dir, pattern = "_2.*\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  
  if(!is.null(test_sample)) {
    r1_files <- r1_files[grepl(test_sample, r1_files)]
    r2_files <- r2_files[grepl(test_sample, r2_files)]
    if(length(r1_files) == 0 || length(r2_files) == 0) stop("Test sample FASTQ not found: ", test_sample)
  }
  
  r1_files <- sort(r1_files)
  r2_files <- sort(r2_files)
  if(length(r1_files) != length(r2_files)) stop("Number of R1 and R2 files does not match!")
  
  for(i in seq_along(r1_files)) {
    sample_name <- get_sample_name(r1_files[i])
    sample_out <- file.path(output_dir, sample_name)
    if(!dir.exists(sample_out)) dir.create(sample_out, recursive = TRUE)
    
    read_cmd <- ifelse(grepl("\\.gz$", r1_files[i]), " --readFilesCommand zcat \\\n", " \\\n")
    
    # Mate1
    mate1_out <- file.path(sample_out, "mate1")
    if(!dir.exists(mate1_out)) dir.create(mate1_out)
    out_prefix1 <- file.path(mate1_out, paste0(sample_name, "_mate1_"))
    
    cmd1 <- paste0(
      star_exe, " \\\n",
      "  --runThreadN ", threads, " \\\n",
      "  --genomeDir ", genome_dir, " \\\n",
      "  --genomeLoad NoSharedMemory \\\n",
      "  --readFilesIn ", r1_files[i], read_cmd,
      "  --outFileNamePrefix ", out_prefix1, " \\\n",
      "  --outReadsUnmapped Fastx \\\n",
      "  --outSAMattributes NH HI AS nM NM MD jM jI XS \\\n",
      "  --outSJfilterOverhangMin 15 15 15 15 \\\n",
      "  --outFilterMultimapNmax 20 \\\n",
      "  --outFilterScoreMin 1 \\\n",
      "  --outFilterMatchNminOverLread 0.7 \\\n",
      "  --outFilterMismatchNmax 999 \\\n",
      "  --outFilterMismatchNoverLmax 0.05 \\\n",
      "  --alignIntronMin 20 \\\n",
      "  --alignIntronMax 1000000 \\\n",
      "  --alignMatesGapMax 1000000 \\\n",
      "  --alignSJoverhangMin 15 \\\n",
      "  --alignSJDBoverhangMin 10 \\\n",
      "  --alignSoftClipAtReferenceEnds No \\\n",
      "  --chimSegmentMin 15 \\\n",
      "  --chimScoreMin 15 \\\n",
      "  --chimScoreSeparation 10 \\\n",
      "  --chimJunctionOverhangMin 15 \\\n",
      "  --sjdbGTFfile ", gtf_file, " \\\n",
      "  --quantMode GeneCounts \\\n",
      "  --twopassMode Basic \\\n",
      "  --chimOutType Junctions SeparateSAMold"
    )
    
    cat("\nMapping mate1 for:", sample_name, "\n")
    cat(cmd1, "\n")
    system(cmd1)
    
    # Mate2
    mate2_out <- file.path(sample_out, "mate2")
    if(!dir.exists(mate2_out)) dir.create(mate2_out)
    cmd2 <- sub(r1_files[i], r2_files[i], cmd1)
    cmd2 <- sub("mate1", "mate2", cmd2)
    cmd2 <- sub("_mate1_", "_mate2_", cmd2)
    
    cat("\nMapping mate2 for:", sample_name, "\n")
    cat(cmd2, "\n")
    system(cmd2)
  }
  
  cat("\nSTAR mate-specific mapping finished. Output:", output_dir, "\n")
  invisible(TRUE)
}

# =========================
# Step 6: Convert SAM -> BAM
# =========================
convert_sam_to_bam <- function(star_dir, mates = c("mate1", "mate2"), convert_joint = TRUE) {
  if(!dir.exists(star_dir)) stop("STAR output directory not found: ", star_dir)
  
  samples <- list.dirs(star_dir, full.names = TRUE, recursive = FALSE)
  
  for(sample_path in samples) {
    sample_name <- basename(sample_path)
    
    if(convert_joint) {
      joint_sam <- list.files(sample_path, pattern = "Aligned.out.sam$", full.names = TRUE)
      if(length(joint_sam) > 0) {
        joint_bam <- file.path(sample_path, "Aligned.noS.bam")
        temp_bam <- file.path(sample_path, "temp.bam")
        cat("\nConverting joint SAM -> BAM for sample:", sample_name, "\n")
        system2("samtools", args = c("view", "-bS", joint_sam, "-o", temp_bam))
        system2("samtools", args = c("sort", "-o", joint_bam, temp_bam))
        file.remove(temp_bam)
        system2("samtools", args = c("index", joint_bam))
      } else message("No joint Aligned.out.sam found for sample: ", sample_name)
    }
    
    for(mate in mates) {
      mate_path <- file.path(sample_path, mate)
      if(!dir.exists(mate_path)) next
      sam_file <- list.files(mate_path, pattern = "Aligned.out.sam$", full.names = TRUE)
      if(length(sam_file) == 0) next
      bam_file <- file.path(mate_path, "Aligned.noS.bam")
      temp_bam <- file.path(mate_path, "temp.bam")
      
      cat("\nConverting", mate, "SAM -> BAM for sample:", sample_name, "\n")
      system2("samtools", args = c("view", "-bS", sam_file, "-o", temp_bam))
      system2("samtools", args = c("sort", "-o", bam_file, temp_bam))
      file.remove(temp_bam)
      system2("samtools", args = c("index", bam_file))
    }
  }
  
  cat("\nSAM -> BAM conversion finished.\n")
  invisible(TRUE)
}

# =========================
# Step 7: Generate Circtools input files
# =========================
generate_circtools_files <- function(base_dir,
                                     output_dir = "circtools_detect") {
  
  base_dir <- path.expand(base_dir)
  output_dir <- path.expand(output_dir)
  
  if(!dir.exists(base_dir)) stop("Base directory not found: ", base_dir)
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Find all 'star' directories recursively
  star_dirs <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
  star_dirs <- star_dirs[basename(star_dirs) == "star"]
  
  if(length(star_dirs) == 0) stop("No STAR output directories found in: ", base_dir)
  
  samplesheet <- c()
  mate1_files <- c()
  mate2_files <- c()
  bam_files <- c()
  
  for(star_dir in star_dirs) {
    sample_path <- star_dir
    sample_name <- basename(dirname(star_dir))  # parent folder = sample name
    
    # Joint Chimeric files
    joint_chimeric <- list.files(sample_path, pattern = "_Chimeric.out.junction$", full.names = TRUE)
    if(length(joint_chimeric) > 0) samplesheet <- c(samplesheet, joint_chimeric[1])
    
    # Joint BAM files
    joint_bam <- list.files(sample_path, pattern = "Aligned.noS.bam$", full.names = TRUE)
    if(length(joint_bam) > 0) bam_files <- c(bam_files, joint_bam[1])
    
    # Mate1 Chimeric files
    mate1_path <- file.path(sample_path, "mate1")
    if(dir.exists(mate1_path)) {
      mate1_chimeric <- list.files(mate1_path, pattern = "_mate1_Chimeric.out.junction$", full.names = TRUE)
      if(length(mate1_chimeric) > 0) mate1_files <- c(mate1_files, mate1_chimeric)
    }
    
    # Mate2 Chimeric files
    mate2_path <- file.path(sample_path, "mate2")
    if(dir.exists(mate2_path)) {
      mate2_chimeric <- list.files(mate2_path, pattern = "_mate2_Chimeric.out.junction$", full.names = TRUE)
      if(length(mate2_chimeric) > 0) mate2_files <- c(mate2_files, mate2_chimeric)
    }
  }
  
  # Write output files
  writeLines(as.character(samplesheet), file.path(output_dir, "samplesheet"))
  writeLines(as.character(mate1_files), file.path(output_dir, "mate1"))
  writeLines(as.character(mate2_files), file.path(output_dir, "mate2"))
  writeLines(as.character(bam_files), file.path(output_dir, "bam_files.txt"))
  
  cat("Generated Circtools input files in:", output_dir, "\n")
  cat("Samples detected:", length(star_dirs), "\n")
  
  invisible(list(samplesheet = samplesheet,
                 mate1 = mate1_files,
                 mate2 = mate2_files,
                 bam_files = bam_files))
}


.run_circtools <- function(sample_row,
                           base_dir,
                           genome_fasta,
                           gtf_file,
                           rRNA_fasta,
                           threads = 10,
                           flexbar_path,
                           bowtie2_path,
                           star_exe,
                           force = FALSE) {
  
  # ✅ Get clean sample name
  fq1 <- sample_row$FileName1
  fq2 <- sample_row$FileName2
  sample_name <- get_sample_name(fq1)                  # derive from FASTQ
  sample_name <- gsub("[^A-Za-z0-9_-]", "_", sample_name)  # replace weird chars
  
  # Create sample directory
  sample_dir <- file.path(base_dir, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  # STAR index folder (shared for all samples)
  star_index_dir <- file.path(base_dir, "star_index")
  
  # Build STAR index once
  generate_star_index(
    genome_fasta = genome_fasta,
    gtf_file = gtf_file,
    output_dir = star_index_dir,
    threads = threads,
    star_exe = star_exe,
    force = FALSE
  )
  
  # Step 1 — Flexbar
  run_flexbar_trim(
    fastq_dir = dirname(fq1),
    output_dir = file.path(sample_dir, "flexbar"),
    test_sample = sample_name,
    flexbar_path = flexbar_path
  )
  
  # Step 2 — rRNA removal
  run_rrna_filter(
    flexbar_dir = file.path(sample_dir, "flexbar"),
    rRNA_fasta = rRNA_fasta,
    output_dir = file.path(sample_dir, "rrna"),
    threads = threads,
    bowtie2_path = bowtie2_path
  )
  
  # Step 3 — STAR joint mapping
  run_star_joint(
    fastq_dir = file.path(sample_dir, "rrna"),
    genome_dir = star_index_dir,
    output_dir = file.path(sample_dir, "star"),
    gtf_file = gtf_file,
    threads = threads,
    star_exe = star_exe
  )
  
  # Step 4 — STAR mate-specific mapping
  run_star_mate_separate(
    fastq_dir = file.path(sample_dir, "rrna"),
    genome_dir = star_index_dir,
    output_dir = file.path(sample_dir, "star"),
    gtf_file = gtf_file,
    threads = threads,
    star_exe = star_exe
  )
  
  # Step 5 — Convert SAM → BAM
  convert_sam_to_bam(file.path(sample_dir, "star"))
  
  message("Finished sample: ", sample_name)
}


run_circtools <- function(sample_file,
                          base_dir,
                          genome_fasta,
                          gtf_file,
                          rRNA_fasta,
                          threads = 4,
                          force = FALSE,
                          flexbar_path = "flexbar",
                          bowtie2_path = "bowtie2",
                          star_exe = "STAR",
                          generate_detect_files = TRUE) {
  
  samples <- read_samplefile(sample_file)
  
  for (i in seq_len(nrow(samples))) {
    
    tryCatch({
      .run_circtools(
        sample_row = samples[i, ],
        base_dir = base_dir,
        genome_fasta = genome_fasta,
        gtf_file = gtf_file,
        rRNA_fasta = rRNA_fasta,
        threads = threads,
        flexbar_path = flexbar_path,
        bowtie2_path = bowtie2_path,
        star_exe = star_exe,
        force = force
      )
    },
    error = function(e) {
      message("ERROR in sample: ", samples$SampleName[i])
      message("Reason: ", e$message)
    })
  }
  
  # ✅ Generate detect input ONCE (after loop)
  if (generate_detect_files) {
    generate_circtools_files(
      base_dir = base_dir,
      output_dir = file.path(base_dir, "circtools_detect")
    )
  }
}