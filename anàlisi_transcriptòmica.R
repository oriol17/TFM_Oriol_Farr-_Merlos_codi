if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligo")
aBiocManager::install("limma")
BiocManager::install("GEOquery")

library(oligo)
library(limma)
library(GEOquery)

# Veure si hi ha diferències entre XM i XP
# Llegim fitxer de Targets
targets_XM_XP <- read.csv("C:/Users/oriol/OneDrive/Documents/Master/TFM/targets_Xm_Xp.csv")
arxiusCel <- paste("C:/Users/oriol/OneDrive/Documents/Master/TFM/Dades transcriptòmica/Dades segon estudi", targets_XM_XP$Nom_Arxiu, sep = "/")
noms_abreviats <- targets_XM_XP$Nom_abreviat
colors_xp_xm <- targets_XM_XP$Color

rawData <- read.celfiles(arxiusCel)

# Boxplot mostres no normalitzades
boxplot(rawData, names=noms_abreviats, col=colors_xp_xm, las = 2, main= "Boxplot de les mostres del cromosoma X d'origen patern (vermell) i matern (verd) sense normalitzar", cex.main=1)

# Agrupació jeràrquica mostres no normalitzades
# Dades d'expressió:
expr_no_norm <- exprs(rawData)
dist_no_norm <- dist(t(expr_no_norm))
hc_no_norm <- hclust(dist_no_norm, method = "average")
plot(hc_no_norm, labels = noms_abreviats, main = "Dendrograma de clustering jeràrquic de M i P", xlab = "", sub = "", cex = 0.8)

# Anàlisi de components principals:
pca_resultat <- prcomp(t(expr_no_norm), scale. = TRUE)

pca_dades <- data.frame(pca_resultat$x[, 1:2])

pca_percentatges <- summary(pca_resultat)

# Extreure els percentatges de variança explicada
percentatge_var_explicada <- pca_percentatges$importance[2, ] * 100

# Veure els percentatges de variança explicada per les dues primeres components principals
percentatge_var_explicada[1:2]


pca_dades$mostres <- noms_abreviats
pca_dades$colors <- colors_xp_xm

library(ggplot2)
ggplot(pca_dades, aes(x = PC1, y = PC2, label = mostres, color = colors)) +
  geom_point(size = 3) +
  geom_text(vjust = 2, hjust = 0.5) +  # Añadir etiquetas
  labs(title = "PCA: primeres dues components principals",
       x = "PC1, 82,17%",
       y = "PC2, 3,4%") +
  theme_minimal()


BiocManager::install("arrayQualityMetrics")
library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = rawData, 
                    outdir = "C:/Users/oriol/OneDrive/Documents/Master/TFM/array_quality_metrics_xp_xm_no_norm_no_log", 
                    force = TRUE, 
                    do.logtransform = FALSE)



# Mateix procés amb les dades normalitzades amb la funció rma:
dades_normalitzades <- rma(rawData)
boxplot(dades_normalitzades, names=noms_abreviats, col=colors_xp_xm, las = 2, main= "Boxplot de les mostres del cromosoma X d'origen patern (vermell) i matern (verd) normalitzades", cex.main=1)

expr_norm <- exprs(dades_normalitzades)
dist_norm <- dist(t(expr_norm))
hc_norm <- hclust(dist_norm, method = "average")
plot(hc_norm, labels = noms_abreviats, main = "Dendrograma de les mostres normalitzades", xlab = "", sub = "", cex = 1)

pca_norm_resultat <- prcomp(t(expr_norm), scale. = TRUE)

pca_dades_norm <- data.frame(pca_norm_resultat$x[, 1:2])

pca_percentatges_norm <- summary(pca_norm_resultat)

percentatge_var_explicada_norm <- pca_percentatges_norm$importance[2, ] * 100

percentatge_var_explicada_norm[1:2]

pca_dades_norm$mostres <- noms_abreviats
pca_dades_norm$colors <- colors_xp_xm

ggplot(pca_dades_norm, aes(x = PC1, y = PC2, label = mostres, color = colors)) +
  geom_point(size = 3) +
  geom_text(vjust = 2, hjust = 0.5) +  # Añadir etiquetas
  labs(title = "PCA: primeres dues components principals amb dades normalitzades",
       x = "PC1, 20%",
       y = "PC2, 12,05%") +
  theme_minimal()

library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = dades_normalitzades, 
                    outdir = "C:/Users/oriol/OneDrive/Documents/Master/TFM/array_quality_metrics_xp_xm_norm", 
                    force = TRUE, 
                    do.logtransform = FALSE)

# Exclusió de la mostra 7:

arxiusCel_sense7 <- arxiusCel[-7]
noms_abreviats_sense7 <- noms_abreviats[-7]
colors_xp_xm_sense7 <- colors_xp_xm[-7]

rawData_sense7 <- read.celfiles(arxiusCel_sense7)

dades_normalitzades_sense7 <- rma(rawData_sense7)
boxplot(dades_normalitzades_sense7, names=noms_abreviats_sense7, col=colors_xp_xm_sense7, las = 2, main= "Boxplot de les mostres del cromosoma X d'origen patern (vermell) i matern (verd) normalitzades, excloent outliers", cex.main=1)

expr_norm_sense7 <- exprs(dades_normalitzades_sense7)
dist_norm_sense7 <- dist(t(expr_norm_sense7))
hc_norm_sense7 <- hclust(dist_norm_sense7, method = "average")
plot(hc_norm_sense7, labels = noms_abreviats_sense7, main = "Dendrograma de les mostres normalitzades sense outliers", xlab = "", sub = "", cex = 1)

pca_norm_resultat_sense7 <- prcomp(t(expr_norm_sense7), scale. = TRUE)

pca_dades_norm_sense7 <- data.frame(pca_norm_resultat_sense7$x[, 1:2])

pca_percentatges_norm_sense7 <- summary(pca_norm_resultat_sense7)

percentatge_var_explicada_norm_sense7 <- pca_percentatges_norm_sense7$importance[2, ] * 100

percentatge_var_explicada_norm_sense7[1:2]

pca_dades_norm_sense7$mostres <- noms_abreviats_sense7
pca_dades_norm_sense7$colors <- colors_xp_xm_sense7

ggplot(pca_dades_norm_sense7, aes(x = PC1, y = PC2, label = mostres, color = colors)) +
  geom_point(size = 3) +
  geom_text(vjust = 2, hjust = 0.5) +  # Añadir etiquetas
  labs(title = "PCA: primeres dues components principals amb dades normalitzades sense outliers",
       x = "PC1, 20,32%",
       y = "PC2, 12,12%") +
  theme_minimal()

library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = dades_normalitzades_sense7, 
                    outdir = "C:/Users/oriol/OneDrive/Documents/Master/TFM/array_quality_metrics_xp_xm_norm_sense_outliers", 
                    force = TRUE, 
                    do.logtransform = FALSE)

# Filtratge no específic:

library(Biobase)
library(genefilter)
library(hgu133plus2.db)

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
BiocManager::install("pd.hg.u133.plus.2.db")

annotation(dades_normalitzades_sense7) <- "hgu133plus2"
exprset_filtrat <- nsFilter(dades_normalitzades_sense7,
                            var.func = IQR,
                            var.cutoff = 0.75,
                            var.filter = TRUE,
                            require.entrez = TRUE,
                            filterByQuantile = TRUE,
                            feature.exclude = "^AFFX")


# Anàlisi d'expressió diferencial:

grups <- targets_XM_XP$Grup[-7]
grups <- factor(grups)
disseny_matriu <- model.matrix(~0 + grups)
colnames(disseny_matriu) <- levels(grups)

matriu_contrastos <- makeContrasts(XP_vs_XM = P - M, levels = disseny_matriu)

fit_disseny <- lmFit(dades_normalitzades_sense7, disseny_matriu)

# Calcular los contrastes
fit_contrastos <- contrasts.fit(fit_disseny, matriu_contrastos)

# Aplicar el método de estimación de varianza bayesiana
fit_bayes <- eBayes(fit_contrastos)

resultats <- topTable(fit_bayes, adjust = "fdr", sort.by = "P", number = Inf)

gens_diferencialment_expressats <- resultats[abs(resultats$logFC) > 1.2 & resultats$adj.P.Val < 0.05, ]




# Segon conjunt de CEL files:

# Llegim fitxer de Targets
targets_complet <- read.csv("C:/Users/oriol/OneDrive/Documents/Master/TFM/targets_complet.csv")
arxiusCel_complet <- paste("C:/Users/oriol/OneDrive/Documents/Master/TFM/Dades transcriptòmica/Dades complet/", targets_complet$Nom_Arxiu, sep = "/")
noms_abreviats_complet <- targets_complet$Nom_abreviat
colors_complet <- targets_complet$Color

rawData_complet <- read.celfiles(arxiusCel_complet)

# Boxplot mostres no normalitzades
boxplot(rawData_complet, names=noms_abreviats_complet, col=colors_complet, las = 2, main= "Boxplot de les mostres XX (vermell) i XO, Síndrome de Turner (verd) sense normalitzar", cex.main=1)

# Agrupació jeràrquica mostres no normalitzades
# Dades d'expressió:
expr_no_norm_complet <- exprs(rawData_complet)
dist_no_norm_complet <- dist(t(expr_no_norm_complet))
hc_no_norm_complet <- hclust(dist_no_norm_complet, method = "average")
plot(hc_no_norm_complet, labels = noms_abreviats_complet, main = "Dendrograma de les mostres en l'estudi que combina GSE46687 y GSE58435", xlab = "", sub = "", cex = 0.8)

# Anàlisi de components principals:
pca_resultat_complet <- prcomp(t(expr_no_norm_complet), scale. = TRUE)

pca_dades_complet <- data.frame(pca_resultat_complet$x[, 1:2])

pca_percentatges_complet <- summary(pca_resultat_complet)

# Extreure percentatges de variança explicada
percentatge_var_explicada_complet <- pca_percentatges_complet$importance[2, ] * 100

# Veure els percentatges de variança explicada pels dos primeres components
percentatge_var_explicada_complet[1:2]


pca_dades_complet$mostres <- noms_abreviats_complet
pca_dades_complet$colors <- colors_complet

library(ggplot2)
ggplot(pca_dades_complet, aes(x = PC1, y = PC2, label = mostres, color = colors)) +
  geom_point(size = 3) +
  geom_text(vjust = 2, hjust = 0.5) +  # Añadir etiquetas
  labs(title = "PCA: primeres dues components principals",
       x = "PC1, 83,22%",
       y = "PC2, 3,73%") +
  theme_minimal()

library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = rawData_complet, 
                    outdir = "C:/Users/oriol/OneDrive/Documents/Master/TFM/array_quality_metrics_complet_no_norm", 
                    force = TRUE, 
                    do.logtransform = FALSE)


# Normalitzem conjunt de dades complet i repetim anàlisi

dades_complet_norm <- rma(rawData_complet)

# Boxplot mostres no normalitzades
boxplot(dades_complet_norm, names=noms_abreviats_complet, col=colors_complet, las = 2, main= "Boxplot de les mostres XX (vermell) i XO, Síndrome de Turner (verd) un cop feta la normalització", cex.main=1)

# Agrupació jeràrquica mostres no normalitzades
# Dades d'expressió:
expr_complet_norm <- exprs(dades_complet_norm)
dist_complet_norm <- dist(t(expr_complet_norm))
hc_complet_norm <- hclust(dist_complet_norm, method = "average")
plot(hc_complet_norm, labels = noms_abreviats_complet, main = "Dendrograma de les mostres normalitzades en l'estudi que combina GSE46687 y GSE58435", xlab = "", sub = "", cex = 0.8)

# Anàlisi de components principals:
pca_resultat_complet_norm <- prcomp(t(expr_complet_norm), scale. = TRUE)

pca_dades_complet_norm <- data.frame(pca_resultat_complet_norm$x[, 1:2])

pca_percentatges_complet_norm <- summary(pca_resultat_complet_norm)

# Extreure percentatges de variança explicada
percentatge_var_explicada_complet_norm <- pca_percentatges_complet_norm$importance[2, ] * 100

# Veure els percentatges de variança explicada pels dos primeres components
percentatge_var_explicada_complet_norm[1:2]


pca_dades_complet_norm$mostres <- noms_abreviats_complet
pca_dades_complet_norm$colors <- colors_complet

library(ggplot2)
ggplot(pca_dades_complet_norm, aes(x = PC1, y = PC2, label = mostres, color = colors)) +
  geom_point(size = 3) +
  geom_text(vjust = 2, hjust = 0.5) +  # Añadir etiquetas
  labs(title = "PCA: primeres dues components principals en dades normalitzades",
       x = "PC1, 42,96%",
       y = "PC2, 7,59%") +
  theme_minimal()

library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = dades_complet_norm, 
                    outdir = "C:/Users/oriol/OneDrive/Documents/Master/TFM/array_quality_metrics_complet_norm", 
                    force = TRUE, 
                    do.logtransform = FALSE)


# Filtratge no específic del dataset complet
library(Biobase)
library(genefilter)
library(hgu133plus2.db)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
BiocManager::install("pd.hg.u133.plus.2.db")

annotation(dades_complet_norm) <- "hgu133plus2"
exprset_filtrat <- nsFilter(dades_complet_norm,
                            var.func = IQR,
                            var.cutoff = 0.75,
                            var.filter = TRUE,
                            require.entrez = TRUE,
                            filterByQuantile = TRUE,
                            feature.exclude = "^AFFX")

# Anàlisi d'expressió diferencial del dataset complet:

grups_complet <- targets_complet$Grup
grups_complet <- factor(grups_complet)
disseny_matriu_complet <- model.matrix(~0 + grups_complet)
colnames(disseny_matriu_complet) <- levels(grups_complet)

matriu_contrastos_complet <- makeContrasts(XO_vs_XX = XO - XX, levels = disseny_matriu_complet)

fit_disseny_complet <- lmFit(dades_complet_norm, disseny_matriu_complet)

# Calcular contrastos
fit_contrastos_complet <- contrasts.fit(fit_disseny_complet, matriu_contrastos_complet)

# Aplicar el mètode d'estimació de variança bayesiana
fit_bayes_complet <- eBayes(fit_contrastos_complet)

resultats_complet <- topTable(fit_bayes_complet, adjust = "fdr", sort.by = "P", number = Inf)

gens_diferencialment_expressats_complet <- resultats_complet[abs(resultats_complet$logFC) > 1.2 & resultats_complet$adj.P.Val < 0.05, ]

nrow(gens_diferencialment_expressats_complet)

library(annotate)
library(hgu133plus2.db)

probes_complet <- row.names(gens_diferencialment_expressats_complet)

# Obtenir Entrez IDs
entrez_ids_complet <- getEG(probes_complet, "hgu133plus2")

# Obtenir noms de gens
gene_names_complet <- getSYMBOL(probes_complet, "hgu133plus2")

# Crear un dataframe amb la informació
df_gens_entrezid_probes <- data.frame(ProbeID = probes_complet, EntrezID = entrez_ids_complet, GeneName = gene_names_complet, p_valor = gens_diferencialment_expressats_complet[4], logFC = gens_diferencialment_expressats_complet[1])

# Mostrar la informació
print(df_gens_entrezid_probes)

resultats_complet_unics <- df_gens_entrezid_probes[!duplicated(df_gens_entrezid_probes$GeneName), ]


# Volcanoplot

genes_mes_expressats <- head(gens_diferencialment_expressats_complet[order(gens_diferencialment_expressats_complet$adj.P.Val), ], 10)

# Obtenir els noms dels 10 gens més diferencialment expressats
noms_gens_mes_expressats <- getSYMBOL(row.names(gens_diferencialment_expressats_complet), "hgu133plus2")

# Generar el volcanoplot y ressaltar els 10 gens més diferencialment expressats
volcanoplot(fit_bayes_complet, highlight = 10, names = noms_gens_mes_expressats, main = paste("Gens diferencialment expressats amb p =< 0.05", colnames(matriu_contrastos_complet), sep = "\n"))

# Afegir les línies de tall per logFC
abline(v = c(-1.2, 1.2))

# Afegir els noms dels 10 gens més diferencialment expressats
with(gens_diferencialment_expressats_complet, text(logFC, -log10(adj.P.Val), labels = noms_gens_mes_expressats, pos = 4, cex = 0.7, col = "red"))


# Heatmaps
install.packages("gplots")
library(gplots)
files_seleccionades <- rownames(expr_norm) %in% rownames(gens_diferencialment_expressats_complet)
dades_seleccionades <- expr_norm[files_seleccionades, ]
coolmap(dades_seleccionades, cluster.by = "de pattern", linkage.row = "complete", linkage.col = "complete", show.dendrogram = "both", main = "Heatmap XO.vs.XX p =< 0.05")


# Anàlisi de significació biològica
if (!requireNamespace("GOstats", quietly = TRUE)) {
  BiocManager::install("GOstats")
}
library(GOstats)
probes_seleccionades <- df_gens_entrezid_probes$ProbeID
entrez_seleccionats <- df_gens_entrezid_probes$EntrezID

#Eliminem duplicats
entrezUniverse <- AnnotationDbi::select(hgu133plus2.db, rownames(expr_norm), "ENTREZID")$ENTREZID

gens_no_duplicats<- entrez_seleccionats[!duplicated(entrez_seleccionats)]

GOparams <- new("GOHyperGParams", geneIds = entrez_seleccionats, universeGeneIds = entrezUniverse,
                 annotation = "hgu133plus2.db", ontology = "BP", pvalueCutoff = 0.05)

# Test de Fisher

hyper_g_test <- hyperGTest(GOparams)

head(summary(hyper_g_test), n = 10)
dim(summary(hyper_g_test))


gens_sign_bio <-resultats_complet["adj.P.Val"]<0.08
IDs_seleccionats <- rownames(resultats_complet)[gens_sign_bio]

EntrezIDs_seleccionats <- AnnotationDbi:::select(hgu133plus2.db, IDs_seleccionats, c("ENTREZID"))

EntrezIDs_seleccionats <- EntrezIDs_seleccionats$ENTREZID

llista_entrezids <- list(EntrezIDs_seleccionats)
nom_taula <- deparse(substitute(resultats_complet))
names(llista_entrezids) <- nom_taula

sapply(llista_entrezids, length)

# Obtenir els gens mapejats a termes GO per homo sapiens
mapped_genes_GO <- mappedkeys(org.Hs.egGO)

# Obtenir els gens mapejats a vies KEGG per Homo sapiens
mapped_genes_KEGG <- mappedkeys(org.Hs.egPATH)

# Combina GO i KEGG
mapped_genes <- union(mapped_genes_GO , mapped_genes_KEGG)


if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  BiocManager::install("ReactomePA")
}

# Carregar ReactomePA
library(ReactomePA)

llista_de_dades <- llista_entrezids
noms_per_comparar <- names(llista_de_dades)
universe <- mapped_genes
for (i in 1:length(llista_de_dades)){
  genesIn <- llista_de_dades[[i]]
  comparacio <- noms_per_comparar[i]
  resultat_enriquit <- enrichPathway(gene = genesIn,
                                  pvalueCutoff = 0.1,
                                  readable = T,
                                  pAdjustMethod = "BH",
                                  organism = "human",
                                  universe = universe)
  head(resultat_enriquit)
}

enrichplot::cnetplot(resultat_enriquit, categorySize = "geneNum", schowCategory = 15, vertex.label.cex = 0.75)
