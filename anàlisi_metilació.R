## Alineació de seqüències amb el genoma de referència:
library(Rsubread)
library(methylKit)
library(ggplot2)

#Definim el path on trobem el genoma, caldrà canviar-ho!!
genome_path = "D:/Master/TFM/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Indexem el genoma de referència
buildindex(basename="hg38_index", reference=genome_path)

# Llista de seqüències fastq
fastq_files <- list("SRR18283144_1_processed.fastq","SRR18283144_2_processed.fastq","SRR18283143_1_processed.fastq","SRR18283143_2_processed.fastq","SRR18283142_1_processed.fastq","SRR18283142_2_processed.fastq","SRR18283141_1_processed.fastq","SRR18283141_2_processed.fastq","SRR18283149_1_processed.fastq","SRR18283149_2_processed.fastq","SRR18283148_1_processed.fastq","SRR18283148_2_processed.fastq","SRR18283147_1_processed.fastq","SRR18283147_2_processed.fastq","SRR18283146_1_processed.fastq","SRR18283146_2_processed.fastq","SRR18283157_1_processed.fastq","SRR18283157_2_processed.fastq","SRR18283156_1_processed.fastq","SRR18283156_2_processed.fastq","SRR18283155_1_processed.fastq","SRR18283155_2_processed.fastq","SRR18283154_1_processed.fastq","SRR18283154_2_processed.fastq")

# Definim directori on les trobarem:

fastq_directory <- "D:/Master/TFM/SRAs"

fastq_paths <- file.path(fastq_directory, fastq_files)


# Definim índex del genoma
index <- "hg38_index"

# Ruta base on escriurem els arxius, s'ha de definir!!!!
base_path <- "D:/path/to/"


# Creem una llista buida per les dades de metilació
methylation_data_list <- list()

for (fastq_path in fastq_paths) {
  # Obtenir nom de la seqüència sense l'extensió fastq
  sample_id <- sub("_processed.fastq$", "", basename(fastq_path))
  
  # Definir rutes d'arxius
  output_bam <- file.path(base_path, paste0(sample_id, ".bam"))
  output_bed <- file.path(base_path, paste0(sample_id, ".bed"))
  
  # Aliniar la seqüència
  align(index=index, readfile1=fastq_path, output_file=output_bam)
  
  # Convertir de BAM a BED
  system(paste("bedtools bamtobed -i", output_bam, ">", output_bed))
  
  # Crear objecte methylRawList a partir de l'arxiu BED
  my_methylation_data <- methRead(file=output_bed, sample.id=sample_id, assembly="hg38", treatment=c(0), context="CpG")
  
  # Generem una variable per cada element de la llista
  variable_name <- paste0("methyl_", sample_id)
  assign(variable_name, my_methylation_data)
  
  # Afegim objecte a la llista
  methylation_data_list[[sample_id]] <- my_methylation_data
  
  # Fem un print del head
  print(head(my_methylation_data))
}

# Definim grups de mostres

ns_BAV_samples <- c("SRR18283144_1", "SRR18283144_2", "SRR18283143_1", "SRR18283143_2", 
                    "SRR18283142_1", "SRR18283142_2", "SRR18283141_1", "SRR18283141_2")

ts_BAV_samples <- c("SRR18283149_1", "SRR18283149_2", "SRR18283148_1", "SRR18283148_2", 
                    "SRR18283147_1", "SRR18283147_2", "SRR18283146_1", "SRR18283146_2")

ts_TAV_samples <- c("SRR18283157_1", "SRR18283157_2", "SRR18283156_1", "SRR18283156_2", 
                    "SRR18283155_1", "SRR18283155_2", "SRR18283154_1", "SRR18283154_2")



# Creem llistes d'objectes de tipus methylRaw per cada grup:

ns_BAV_list <- lapply(ns_BAV_samples, function(x) get(paste0("methyl_", x)))
ts_BAV_list <- lapply(ts_BAV_samples, function(x) get(paste0("methyl_", x)))
ts_TAV_list <- lapply(ts_TAV_samples, function(x) get(paste0("methyl_", x)))



# Convinem les llistes d'objectes methylRaw en objectes methylRawList
ns_BAV_methylRawList <- new("methylRawList", ns_BAV_list)
ts_BAV_methylRawList <- new("methylRawList", ts_BAV_list)
ts_TAV_methylRawList <- new("methylRawList", ts_TAV_list)


# Unir els 3 grups en una llista
combined_list <- list(ns_BAV_methylRawList, ts_BAV_methylRawList, ts_TAV_methylRawList)

# Unir les dades en un sol objecte
combined_methylation_data <- unite(combined_list)

# Fer la normalització
normalized_data <- normalizeCoverage(combined_methylation_data)

# Analitzar la diferència de metilació entre grups
diff_methylation <- calculateDiffMeth(normalized_data, covariates = NULL, overdispersion = "MN", test = "Chisq")

# Obtenir les regions especialment metilades
methylation_differentially_expressed <- getMethylDiff(diff_methylation, difference = 25, qvalue = 0.01)


# Histograma per veure les diferències de metilació
hist(methylation_differentially_expressed$diff.meth, main = "Distribució de les diferències de metilació", xlab = "Diferència de metilació (%)", breaks = 50)


# Gràfic de distribució de metilació
met_levels <- getMeth(normalized_data, type = "raw")
df <- data.frame(
  sample = rep(c("ns_BAV", "ts_BAV", "ts_TAV"), each = length(met_levels[[1]])),
  methylation = unlist(met_levels)
)

ggplot(df, aes(x = sample, y = methylation)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribució dels nivells de metilació per grup", x = "Grup", y = "Nivell de metilació")
