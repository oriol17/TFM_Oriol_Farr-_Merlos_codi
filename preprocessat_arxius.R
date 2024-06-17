library(ShortRead)
library(Biostrings)
library(parallel)

# Directori a definir on guardarem els arxius MODIFICAR!!! Ull viu, no posar / al final, que ja n'afegeix una file.path
fastq_directory <- "/home/jaume/Documentos/uoc_prova"
fastq_files = list.files(pattern = ".fastq")

# Llista d'arxius FASTQ
#fastq_files <- list("SRR18283144_1.fastq","SRR18283144_2.fastq","SRR18283143_1.fastq","SRR18283143_2.fastq","SRR18283142_1.fastq","SRR18283142_2.fastq","SRR18283141_1.fastq","SRR18283141_2.fastq","SRR18283149_1.fastq","SRR18283149_2.fastq","SRR18283148_1.fastq","SRR18283148_2.fastq","SRR18283147_1.fastq","SRR18283147_2.fastq","SRR18283146_1.fastq","SRR18283146_2.fastq","SRR18283157_1.fastq","SRR18283157_2.fastq","SRR18283156_1.fastq","SRR18283156_2.fastq","SRR18283155_1.fastq","SRR18283155_2.fastq","SRR18283154_1.fastq","SRR18283154_2.fastq")

# Ajuntem path amb nom dels arxius fastq:

fastq_paths <- file.path(fastq_directory, fastq_files)


# 1- Primer Trim
trim_sequences <- function(fastq_file, cutoff = "5", minlen = 50) {
  trimmed_fastq <- trimTails(fastq_file, k = 2, a = cutoff)
  keep <- width(sread(trimmed_fastq)) >= minlen
  keep[is.na(keep)] <- FALSE  # Manejar NA
  trimmed_fastq <- trimmed_fastq[keep]
  return(trimmed_fastq)
}

# 2- Funció per eliminar adaptadors
adapter_sequence_left <- "AGATCGGAAGAGCACACGTCTGAAC"
adapter_sequence_right <- "AGATCGGAAGAGCGTCGTGTAGGGA"

trim_adapters <- function(fastq_data, adapter_sequence_left, adapter_sequence_right) {
  adapter_left <- DNAString(adapter_sequence_left)
  adapter_right <- DNAString(adapter_sequence_right)
  trimmed <- trimLRPatterns(Lpattern = adapter_left, Rpattern = adapter_right, subject = sread(fastq_data))
  ids <- id(fastq_data)
  original_quality <- quality(fastq_data)
  quality <- narrow(original_quality, start=1, end=width(trimmed))
  trimmed_fastq <- ShortReadQ(sread = trimmed, quality = quality, id = ids)
  return(trimmed_fastq)
}


# 3- Funció per filtrar seqüències de baixa qualitat

filter_low_quality <- function(fastq_data, min_quality = 20) {
  quality_filter <- apply(as(quality(fastq_data), "matrix"), 1, function(x) all(x >= min_quality))
  quality_filter[is.na(quality_filter)] <- FALSE  # Manejar NA
  filtered <- fastq_data[quality_filter]
  return(filtered)
}


# 4- Funció per eliminar duplicats

remove_duplicates <- function(fastq_data) {
  sequences <- sread(fastq_data)
  unique_indices <- !duplicated(as.character(sequences))
  unique_fastq <- fastq_data[unique_indices]
  return(unique_fastq)
}


# 5- Funció per normalitzar el contingut de GC 

normalize_gc_content <- function(fastq_data) { 
  gc_content <- letterFrequency(sread(fastq_data), "GC", as.prob = TRUE) 
  gc_mean <- mean(gc_content) 
  gc_sd <- sd(gc_content) 
  gc_filter <- gc_content >= (gc_mean - 3 * gc_sd) & gc_content <= (gc_mean + 3 * gc_sd) 
  gc_filter[is.na(gc_filter)] <- FALSE  # Manejar NA
  normalized_fastq <- fastq_data[gc_filter] 
  return(normalized_fastq) 
}


# Funció per processar els arxius
process_fastq_file <- function(fastq_file) {
  #fastq_data <- readFastq(fastq_file)
  # Aplicar funcions de processat
  trimmed_fastq <- trim_sequences(fastq_file)
  trimmed_adapters_fastq <- trim_adapters(trimmed_fastq, adapter_sequence_left, adapter_sequence_right)
  filtered_fastq <- filter_low_quality(trimmed_adapters_fastq)
  final_trimmed_fastq <- trim_sequences(filtered_fastq, cutoff = "5", minlen = 50)
  unique_fastq <- remove_duplicates(final_trimmed_fastq)
  normalized_fastq <- normalize_gc_content(unique_fastq)
}

# Processar tots els arxius de la llista

read_write_fastq = function(fastq_file){
fastq_data <- FastqStreamer(fastq_file, n = 1e6)  # Leer en fragmentos de 1 millón de lecturas
output_file <- paste0(tools::file_path_sans_ext(fastq_file), "_processed.fastq")
while (length(vector <- yield(fastq_data)) > 0) {
  normalized_fastq = process_fastq_file(vector)
  writeFastq(normalized_fastq, file = output_file, mode = "a", compress = TRUE)
  }
close(fastq_data)
}

#unlist(lapply(fastq_paths,process_fastq_file))
mc.cores = detectCores()
mclapply(fastq_paths, read_write_fastq, mc.cores = 16)



