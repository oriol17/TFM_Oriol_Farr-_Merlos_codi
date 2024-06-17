library(Rsubread)
library(methylKit)
library(ggplot2)


# Aquí poseu el directori on tingueu els BEDs
beds_directory <- "C:/Users/oriol/OneDrive/Documents/Master/TFM/beds/beds"
llista_beds <- c("SRR18283141.bed","SRR18283142.bed","SRR18283143.bed","SRR18283144.bed","SRR18283146.bed","SRR18283147.bed","SRR18283148.bed","SRR18283149.bed","SRR18283154.bed","SRR18283155.bed","SRR18283156.bed","SRR18283157.bed")

beds_paths <- file.path(beds_directory, llista_beds)


methylation_data_list <- list()

for (bed_file in beds_paths) {
  # Asignar nombres a las columnas
  bed_data <- read.table(bed_file, header=FALSE, sep="\t")
  colnames(bed_data) <- c("chromosome", "start", "end", "coverage", "strand")
  
  # Generar valores ficticios para methylated_coverage y unmethylated_coverage
  set.seed(42) # Para reproducibilidad
  bed_data$methylated_coverage <- round(bed_data$coverage * runif(nrow(bed_data), 0.3, 0.7))
  bed_data$unmethylated_coverage <- bed_data$coverage - bed_data$methylated_coverage
  final_bed_data <- bed_data[, c("chromosome", "start", "end", "strand", "coverage", "methylated_coverage", "unmethylated_coverage")]
  output_bed_final <- paste("C:/Users/oriol/OneDrive/Documents/Master/TFM/beds/beds/processed/", basename(bed_file), sep = "")
  write.table(final_bed_data, file=output_bed_final, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  sample_id <- sub(".bed$", "", basename(bed_file))

  # Crear objecte methylRawList a partir de l'arxiu BED
  my_methylation_data <- methRead(location=output_bed_final, sample.id=sample_id, assembly="hg38", treatment=c(0), context="CpG")

  # Generem una variable per cada element de la llista
  variable_name <- paste0("methyl_", sample_id)
  assign(variable_name, my_methylation_data)

  # Afegim objecte a la llista
  methylation_data_list[[sample_id]] <- my_methylation_data

  # Fem un print del head
  print(head(my_methylation_data))
}

# Definim grups de mostres

ns_BAV_samples <- c("SRR18283144", "SRR18283143", "SRR18283142", "SRR18283141")

ts_BAV_samples <- c("SRR18283149", "SRR18283148", "SRR18283147", "SRR18283146")

ts_TAV_samples <- c("SRR18283157", "SRR18283156","SRR18283155", "SRR18283154")

#samples <- c("SRR18283144", "SRR18283143", "SRR18283142", "SRR18283141", "SRR18283149", "SRR18283148", "SRR18283147", "SRR18283146", "SRR18283157", "SRR18283156","SRR18283155", "SRR18283154")

# Creem llistes d'objectes de tipus methylRaw per cada grup:

ns_BAV_list <- lapply(ns_BAV_samples, function(x) get(paste0("methyl_", x)))
ts_BAV_list <- lapply(ts_BAV_samples, function(x) get(paste0("methyl_", x)))
ts_TAV_list <- lapply(ts_TAV_samples, function(x) get(paste0("methyl_", x)))

#methyl_raw_list_samples <- new("methylRawList", methylation_data_list)
#class(methyl_raw_list_samples)

# Combinem les llistes d'objectes methylRaw en objectes methylRawList
ns_BAV_methylRawList <- new("methylRawList", ns_BAV_list)
ts_BAV_methylRawList <- new("methylRawList", ts_BAV_list)
ts_TAV_methylRawList <- new("methylRawList", ts_TAV_list)


# Unir els 3 grups en una llista
#combined_list <- c(ns_BAV_methylRawList, ts_BAV_methylRawList, ts_TAV_methylRawList)
#
#class(methyl_raw_list_samples)
## Unir les dades en un sol objecte
#combined_methylation_data <- unite(methyl_raw_list_samples, destrand=FALSE, min.per.group=2L)
#
#summary(methyl_raw_list_samples)


# Normalitzar les dades
norm_ns_BAV <- normalizeCoverage(ns_BAV_methylRawList)
norm_ts_BAV <- normalizeCoverage(ts_BAV_methylRawList)
norm_ts_TAV <- normalizeCoverage(ts_TAV_methylRawList)

# Filtrar les dades
filtered_ns_BAV <- filterByCoverage(norm_ns_BAV, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
filtered_ts_BAV <- filterByCoverage(norm_ts_BAV, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
filtered_ts_TAV <- filterByCoverage(norm_ts_TAV, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

#unir_un_grup <- unite(filtered_ns_BAV, destrand=FALSE)

# Combinar y crear objecte methylBase
combined_methylRawList <- new("methylRawList", c(filtered_ns_BAV@.Data, filtered_ts_BAV@.Data, filtered_ts_TAV@.Data))

# Asignar tratamientos: 0 para ns_BAV, 1 para ts_BAV, 2 para ts_TAV
treatment <- c(rep(0, length(filtered_ns_BAV@.Data)), rep(1, length(filtered_ts_BAV@.Data)), rep(2, length(filtered_ts_TAV@.Data)))

# Unir les dades per a l'anàlisi diferencial
united_methylation <- unite(combined_methylRawList, destrand=FALSE, min.per.group=2L, allow.cartesian=TRUE, by=.EACHI)


## Fins aquí, la part de dalt és el que falla, el unite no funciona


# Anàlisi diferencial de metilació
diff_meth <- calculateDiffMeth(united_methylation, covariates=data.frame(treatment), overdispersion="MN", test="Chisq")

# Resultats
diff_meth_results <- getMethylDiff(diff_meth)
head(diff_meth_results)

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
