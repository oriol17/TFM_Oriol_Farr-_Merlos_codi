if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("anilchalisey/rseqR")
library(rseqR)

devtools::install_github("anilchalisey/rseqR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
library(ShortRead)
library(Biostrings)

fastq_file1 <- "D:/Master/TFM/SRAs/SRR18283141_1.fastq"
fastq_file2 <- "D:/Master/TFM/SRAs/SRR18283141_2.fastq"

fastq_data1 <- readFastq(fastq_file1)
fastq_data2 <- readFastq(fastq_file2)

trim_sequences <- function(fastq_data, cutoff = "5", minlen = 50) {
  trimmed_fastq <- trimTails(fastq_data, k = 2, a = cutoff)
  keep <- width(sread(trimmed_fastq)) >= minlen
  trimmed_fastq <- trimmed_fastq[keep]
  return(trimmed_fastq)
}


# 2- Funció per eliminarlos
adapter_sequence_left <- "AGATCGGAAGAGCACACGTCTGAAC"
adapter_sequence_right <- "AGATCGGAAGAGCGTCGTGTAGGGA"

trim_adapters <- function(fastq_data, adapter_sequence_left, adapter_sequence_right) {
  adapter_left <- DNAString(adapter_sequence_left)
  adapter_right <- DNAString(adapter_sequence_right)
  trimmed <- trimLRPatterns(Lpattern = adapter_left, Rpattern = adapter_right, subject = sread(fastq_data))
  ids <- id(fastq_data)
  #quality <- quality(fastq_data)
  original_quality <- quality(fastq_data)
  quality <- narrow(original_quality, start=1, end=width(trimmed))
  trimmed_fastq <- ShortReadQ(sread = trimmed, quality = quality, id = ids)
  return(trimmed_fastq)
}

# Funció per filtrar seqüències de baixa qualitat
filter_low_quality <- function(fastq_data, min_quality = 20) {
  quality_filter <- apply(as(quality(fastq_data), "matrix"), 1, function(x) all(x >= min_quality))
  filtered <- fastq_data[quality_filter]
  return(filtered)
}

#filter_low_quality <- function(fastq_data, min_quality = 20) {
#  # Obtener la calidad en formato matrix directamente
#  quality_matrix <- as(quality(fastq_data), "matrix")
#  
#  # Convertir la calidad de formato ASCII a valores Phred
#  quality_numeric <- matrix(as.numeric(charToRaw(as.character(quality_matrix))) - 33, nrow = nrow(quality_matrix), ncol = ncol(quality_matrix))
#  
#  # Filtrar las secuencias basadas en la calidad mínima
#  quality_filter <- apply(quality_numeric, 1, function(x) all(x >= min_quality))
#  
#  # Filtrar el fastq_data usando el filtro de calidad
#  filtered <- fastq_data[quality_filter]
#  return(filtered)
#}

# Funció per eliminar duplicats
remove_duplicates <- function(fastq_data) {
  unique_fastq <- unique(fastq_data)
  return(unique_fastq)
}

remove_duplicates <- function(fastq_data) {
  # Obtenir les seqüències i convertir-les a un caràcter per comparar
  sequences <- sread(fastq_data)
  
  # Trobar posicions úniques basades en les seqüències
  unique_indices <- !duplicated(as.character(sequences))
  
  # Filtrar el fastq_data per mantenir només les seqüències úniques
  unique_fastq <- fastq_data[unique_indices]
  
  return(unique_fastq)
}

#   Obtenció final
trimmed_fastq1 <- trim_sequences(fastq_data1)
trimmed_fastq2 <- trim_sequences(fastq_data2)

# Retalla adaptadors
trimmed_adapters_fastq <- trim_adapters(trimmed_fastq2, adapter_sequence_left, adapter_sequence_right)

# Filtra seqüències de baixa qualitat
filtered_fastq <- filter_low_quality(trimmed_adapters_fastq)

# Retalla les cues de les seqüències amb la funció trim_sequences

#final_trimmed_fastq <- trim_sequences(filtered_fastq, cutoff = "5", minlen = 50)
final_trimmed_fastq <- trim_sequences(trimmed_adapters_fastq, cutoff = "5", minlen = 50)

# Elimina duplicats
unique_fastq <- remove_duplicates(final_trimmed_fastq)

# Escriu l'arxiu FASTQ processat
output_file <- "path/to/your/output_file.fastq"
writeFastq(unique_fastq, file = output_file2, compress = FALSE)

# Comprova la qualitat amb FastQC
system(paste("fastqc", output_file))

# Funció per normalitzar el contingut de GC 
normalize_gc_content <- function(fastq_data) { 
  gc_content <- letterFrequency(sread(fastq_data), "GC", as.prob = TRUE) 
  gc_mean <- mean(gc_content) 
  gc_sd <- sd(gc_content) 
# Filtra seqüències amb un contingut de GC fora de 3 desviacions estàndard de la mitjana 
  gc_filter <- gc_content >= (gc_mean - 3 * gc_sd) & gc_content <= (gc_mean + 3 * gc_sd) 
  normalized_fastq <- fastq_data[gc_filter] 
  return(normalized_fastq) 
}


# Normalitza el contingut de GC 
normalized_fastq <- normalize_gc_content(unique_fastq)

writeFastq(normalized_fastq, file = output_file2, compress = FALSE)
#cutoff <- "5"
#trimmed_fastq <- trimTails(fastq_data2, a = cutoff, k = 2)
#minlen <- 50
#keep <- width(sread(trimmed_fastq)) >= minlen
#trimmed_fastq <- trimmed_fastq[keep]
#
#cutoff <- "20"
#k <- 2
#trimmed_fastq <- trimTailw(fastq_data2, k = k, a = cutoff, halfwidth = 2)
#minlen <- 50
#keep <- width(sread(trimmed_fastq)) >= minlen
#trimmed_fastq <- trimmed_fastq[keep]

trimmed_fastq1 <- trim_sequences(fastq_data1)
trimmed_fastq2 <- trim_sequences(fastq_data2)

output_file1 <- "D:/Master/TFM/SRAs/SRR18283141_1_trimmed_by_r.fastq"
output_file2 <- "D:/Master/TFM/SRAs/SRR18283141_2_trimmed_by_r5.fastq"
writeFastq(trimmed_fastq1, output_file1)
writeFastq(trimmed_fastq2, file = output_file2, compress = FALSE)



## Aliniació de seqüències amb el genoma de referència:
BiocManager::install(c("Rsubread", "methylKit"))
library(Rsubread)
library(methylKit)

# Indexem el genoma de referència
buildindex(basename="hg38_index", reference="D:/Master/TFM/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa", memory="4G")


# Aliniem la secuencia
align(index="hg38_index", readfile1="D:/path/to/sequence.fasta", output_file="aligned.bam")


# Convertim de BAM a BED
system("bedtools bamtobed -i aligned.bam > aligned.bed")

# Crear objecte methylRawList a partir de l'arxiu BED
my_methylation_data <- read(file="aligned.bed", sample.id="sample1", assembly="hg38", treatment=c(0), context="CpG")

# Inspeccionar dades
head(my_methylation_data)


# Anàlisi de metilació:

# Filtrar dades
filtered_data <- filterByCoverage(my_methylation_data, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

# Normalizar dades
normalized_data <- normalizeCoverage(filtered_data)

# Calcular la diferencia de metilación
meth_diff <- calculateDiffMeth(normalized_data, overdispersion="MN", test="Chisq", mc.cores=1)

# Obtenir les regions diferencialment metilades
diff_meth_regions <- getMethylDiff(meth_diff, difference=25, qvalue=0.01)

# Visualitzar les regions diferencialment metilades
head(diff_meth_regions)