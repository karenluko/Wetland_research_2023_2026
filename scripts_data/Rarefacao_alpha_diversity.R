setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data")

library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(janitor)
library(tidyverse)
library(ggrepel)

dados <- read.delim("feature-table_Wetland 1.tsv",
                    row.names = 1,
                    header = TRUE,
                    sep = "\t")

dados <- dados %>%
        row_to_names(row_number = 1)

dados <- dados[, 1:22]

dados <- as.data.frame(lapply(dados, as.numeric))

dados <- dados[, colSums(dados) > 0]


dados_matriz <- as.matrix(dados)

## -----------------------------
## Metadados das amostras
## -----------------------------

metadata <- data.frame(
        Sample  = colnames(dados_matriz),
        Wetland = c(rep("Permanent", 6), rep("Seasonal", 6),
                    rep("Permanent", 5), rep("Seasonal", 5)), # S14 e S19 corrigidos
        Depth   = rep(c("0-10 cm", "10-20 cm"), times = 11)
)

## -----------------------------
## Normalização (rarefação a mesma profundidade)
## -----------------------------


lib_sizes <- colSums(dados_matriz)
summary(lib_sizes)

# Profundidade mínima para rarefação
min_depth <- min(lib_sizes)
min_depth

set.seed(123) 

# Rarefaz todas as amostras para min_depth
dados_rarefied <- t(
        rrarefy(t(dados_matriz), sample = min_depth)
)

## -----------------------------
##sequencing depth x ASVs
## -----------------------------

# Sequência de esforços amostrais até a profundidade mínima entre amostras
esforcos <- seq(1, min_depth, length.out = 100)
esforcos <- round(esforcos)

rarefacao_dados <- lapply(1:ncol(dados_matriz), function(i) {
        max_reads <- sum(dados_matriz[, i])  
        esforcos_validos <- esforcos[esforcos <= max_reads] 
        otu_vals <- as.numeric(
                rarefy(dados_matriz[, i], 
                       sample = esforcos_validos,
                       se = FALSE)
        )
        
        data.frame(
                Sample = colnames(dados_matriz)[i],
                Effort = esforcos_validos,
                OTUs   = otu_vals
        ) %>%
                left_join(metadata, by = "Sample")
})


dados_long <- bind_rows(rarefacao_dados) %>% drop_na()

# Gráfico 1 – por tipo de wetland

r1 <- ggplot(dados_long,
             aes(x = Effort, y = OTUs,
                 color = Wetland, group = Sample)) +
        geom_line(linewidth = 1) +
        facet_wrap(~ Depth) +
        theme_light() +
        labs(
                title = "",
                x = "Sequenced reads",
                y = "Number of observed ASVs",
                color = "Wetland type"
        ) +
        scale_color_manual(values = c("Permanent" = "turquoise3",
                                      "Seasonal"  = "violetred3")) +
        theme(
                legend.position = "bottom",
                text = element_text(size = 12)
        )

print(r1)

ggsave("rarefacao_wetland.tiff",
       plot = r1,
       width = 10,
       height = 8,
       dpi = 300,
       units = "in",
       device = "tiff")

# Gráfico 2 – amostras individuais

r2 <- ggplot(dados_long,
             aes(x = Effort, y = OTUs,
                 color = Sample,
                 linetype = Wetland,
                 group = Sample)) +
        geom_line(linewidth = 1) +
        facet_wrap(~ Depth) +
        theme_light() +
        labs(
                title = "",
                x = "Sequenced reads",
                y = "Number of observed ASVs",
                color = "Sample",
                linetype = "Wetland type"
        ) +
        scale_linetype_manual(values = c("Permanent" = "solid",
                                         "Seasonal"  = "dashed")) +
        theme(
                legend.position = "bottom",
                text = element_text(size = 12)
        )

print(r2)

ggsave("rarefacao_samples.tiff",
       plot = r2,
       width = 10,
       height = 8,
       dpi = 300,
       units = "in",
       device = "tiff")
##########################################################

# Índices de diversidade (tabela rarefeita)


library(vegan)  

comm <- t(dados_rarefied)

# Observed richness, Chao1, Shannon e Simpson
Observed <- specnumber(comm)                    
chao_mat <- estimateR(comm)                     
Chao1   <- chao_mat["S.chao1", ]

Shannon  <- diversity(comm, index = "shannon")
Simpson  <- diversity(comm, index = "simpson")

dados_diversidade <- data.frame(
        Sample          = rownames(comm),
        Observed        = Observed,
        Chao1           = Chao1,
        Shannon         = Shannon,
        Simpson         = Simpson
)


amostras <- c("S1", "S10", "S11", "S12", "S13", "S15", "S16", "S17", "S18", 
              "S2", "S20", "S21", "S22", "S23", "S24", "S3", "S4", "S5", 
              "S6", "S7", "S8", "S9")

Wetland <- c("Permanent", rep("Seasonal", 3), rep("Permanent", 6),
             rep("Seasonal", 5), rep("Permanent", 4), rep("Seasonal", 3))

Depth <- c("0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "0-10 cm", "0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "10-20 cm", "10-20 cm", "0-10 cm", "10-20 cm", "0-10 cm",
           "10-20 cm", "0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "0-10 cm", "10-20 cm", "0-10 cm")

Season <- c(rep("Rainy", 4),
            rep("Dry", 5),
            rep("Rainy", 1),
            rep("Dry", 5),
            rep("Rainy", 7))

metadata_alpha <- data.frame(
        Sample  = amostras,
        Wetland = Wetland,
        Depth   = Depth,
        Season  = Season
)

dados_diversidade <- dados_diversidade %>%
        left_join(metadata_alpha, by = "Sample")

dados_diversidade$Sample <- factor(
        dados_diversidade$Sample,
        levels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                   "S10", "S11", "S12", "S13", "S15", "S16", "S17", "S18",
                   "S20", "S21", "S22", "S23", "S24")
)

dados_diversidade$Season <- factor(
        dados_diversidade$Season,
        levels = c("Rainy", "Dry")
)

dados_long_div <- dados_diversidade %>%
        pivot_longer(
                cols = c(Observed, Chao1, Shannon, Simpson),
                names_to = "Index",
                values_to = "Value"
        ) %>%
        mutate(
                Index = factor(
                        Index,
                        levels = c("Observed", "Chao1", "Shannon", "Simpson"),
                        labels = c("Observed richness", "Chao1", "Shannon", "Simpson")
                )
        )

x <- ggplot(dados_long_div,
            aes(x = Depth, y = Value, color = Wetland)) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
        facet_grid(Index ~ Wetland + Season, scales = "free_y") +
        theme_light() +
        labs(
                title = "",
                x = "Depth",
                y = "Index",
                color = "Wetland Type"
        ) +
        scale_color_manual(values = c("Permanent" = "turquoise3",
                                      "Seasonal"  = "violetred3")) +
        theme(
                legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                text = element_text(size = 12)
        )

plot(x)

ggsave("diversidade_indice_rarefied.tiff",
       plot = x,
       width = 10,
       height = 8,
       dpi = 300,
       units = "in",
       device = "tiff")

######################################################
#ANOVA

##Shannon
anova_shannon <- aov(Shannon ~ Wetland + Season + Depth, data = dados_diversidade)
summary(anova_shannon)
TukeyHSD(anova_shannon)

#Observed richness
anova_observed <- aov(Observed ~ Wetland + Season + Depth, data = dados_diversidade)
summary(anova_observed)
TukeyHSD(anova_observed)

#Chao1
anova_chao1 <- aov(Chao1 ~ Wetland + Season + Depth, data = dados_diversidade)
summary(anova_chao1)
TukeyHSD(anova_chao1)

#Simpson
anova_simpson <- aov(Simpson ~ Wetland + Season + Depth, data = dados_diversidade)
summary(anova_simpson)
TukeyHSD(anova_simpson)

write.csv(dados_rarefied, "dados_rarefied_table.csv", row.names = TRUE)

#######################################
#Beta diversity (todos juntos)

library(compositions)
library(vegan)
library(dplyr)
library(ggplot2)

comm <- t(dados_rarefied)  
comm <- as.matrix(comm)

dim(comm)       
rownames(comm)


comm_nz <- comm + 1

# CLR transformation (Aitchison space)
comm_clr <- clr(comm_nz)     

# Distância Euclidiana sobre os dados CLR
dist_clr <- dist(comm_clr, method = "euclidean")

metadata <- dados_diversidade %>%
        select(Sample, Wetland, Season, Depth) %>%
        distinct()

# Garantir que a ordem de Sample bate com a ordem de rownames(comm_clr)
metadata <- metadata %>%
        arrange(match(Sample, rownames(comm_clr)))

all(metadata$Sample == rownames(comm_clr))

# PERMANOVA com Aitchison (CLR + dist Euclidiana)
permanova_clr <- adonis2(
        dist_clr ~ Wetland + Season + Depth,
        data = metadata
)

permanova_clr

#PERMDISP (checar dispersão entre grupos)
bd_wetland <- betadisper(dist_clr, metadata$Wetland)
anova(bd_wetland)

# p > 0.05 → dispersão homogênea

pcoa_clr <- cmdscale(dist_clr, eig = TRUE, k = 2)

scores <- as.data.frame(pcoa_clr$points)
colnames(scores) <- c("PCoA1", "PCoA2")
scores$Sample <- rownames(scores)

scores <- scores %>%
        left_join(metadata, by = "Sample")

ggplot(scores, aes(x = PCoA1, y = PCoA2,
                   color = Wetland,
                   shape = Season)) +
        geom_point(size = 3) +
        theme_light() +
        labs(
                x = "PCoA1",
                y = "PCoA2",
                color = "Wetland type",
                shape = "Season"
        )

ggplot(scores, aes(x = PCoA1, y = PCoA2,
                   color = Wetland,
                   shape = Season)) +
        geom_point(size = 3) +
        facet_wrap(~ Depth) +
        theme_light()
#########################################
#Betadiversity separado bacteria e archeae

library(dplyr)
library(compositions)
library(vegan)
library(readr)


dados <- read.delim("feature-table_Wetland 1.tsv",
                    row.names = 1,
                    header = TRUE,
                    sep = "\t")


dados <- dados %>%
        row_to_names(row_number = 1)


dados <- dados[, 1:22]


dados <- as.data.frame(lapply(dados, as.numeric))


dados <- dados[, colSums(dados) > 0]


dados_matriz <- as.matrix(dados)

# Metadados das amostras

amostras <- c("S1", "S10", "S11", "S12", "S13", "S15", "S16", "S17", "S18", 
              "S2", "S20", "S21", "S22", "S23", "S24", "S3", "S4", "S5", 
              "S6", "S7", "S8", "S9")

Wetland <- c("Permanent", rep("Seasonal", 3), rep("Permanent", 6),
             rep("Seasonal", 5), rep("Permanent", 4), rep("Seasonal", 3))

Depth <- c("0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "0-10 cm", "0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "10-20 cm", "10-20 cm", "0-10 cm", "10-20 cm", "0-10 cm",
           "10-20 cm", "0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "0-10 cm", "10-20 cm", "0-10 cm")

Season <- c(rep("Rainy", 4),
            rep("Dry", 5),
            rep("Rainy", 1),
            rep("Dry", 5),
            rep("Rainy", 7))

metadata <- data.frame(
        Sample  = amostras,
        Wetland = Wetland,
        Depth   = Depth,
        Season  = Season
)


arch <- read_csv2("Class_archeae.csv")
bact <- read_csv2("Class_bacteria.csv")


sample_names <- metadata$Sample 
# Archaea
arch_samples <- arch %>%
        select(Class, all_of(sample_names))

arch_mat <- as.matrix(arch_samples[ , -1])        
rownames(arch_mat) <- arch_samples$Class          

# Bacteria
bact_samples <- bact %>%
        select(Class, all_of(sample_names))

bact_mat <- as.matrix(bact_samples[ , -1])
rownames(bact_mat) <- bact_samples$Class


arch_mat_t <- t(arch_mat)   
bact_mat_t <- t(bact_mat)   


all(rownames(arch_mat_t) == metadata$Sample)
all(rownames(bact_mat_t) == metadata$Sample)

## 3. Pseudocount + CLR

arch_clr <- clr(arch_mat_t + 1)
bact_clr <- clr(bact_mat_t + 1)

## 4. Distâncias (Aitchison = Euclidiana no CLR)

arch_dist <- dist(arch_clr, method = "euclidean")
bact_dist <- dist(bact_clr, method = "euclidean")

## 5. PERMANOVA

adonis2(arch_dist ~ Wetland + Season + Depth, data = metadata)
adonis2(bact_dist ~ Wetland + Season + Depth, data = metadata)

anova(betadisper(arch_dist, metadata$Wetland))
anova(betadisper(bact_dist, metadata$Wetland))

anova(betadisper(arch_dist, metadata$Season))
anova(betadisper(bact_dist, metadata$Season))

anova(betadisper(arch_dist, metadata$Depth))
anova(betadisper(bact_dist, metadata$Depth))

#p > 0.05 → dispersão homogênea

library(dplyr)
library(ggplot2)

## PCoA – Archaea
arch_pcoa <- cmdscale(arch_dist, eig = TRUE, k = 2)

scores_arch <- as.data.frame(arch_pcoa$points)
colnames(scores_arch) <- c("PCoA1", "PCoA2")
scores_arch$Sample <- rownames(scores_arch)

scores_arch <- scores_arch %>%
        left_join(metadata, by = "Sample") %>%
        mutate(Domain = "Archaea")  # para usar depois se quiser combinar

ggplot(scores_arch, aes(x = PCoA1, y = PCoA2,
                        color = Wetland,
                        shape = Season)) +
        geom_point(size = 3) +
        theme_light() +
        labs(
                x = "PCoA1",
                y = "PCoA2",
                color = "Wetland type",
                shape = "Season"
        )

ggplot(scores_arch, aes(x = PCoA1, y = PCoA2,
                        color = Wetland,
                        shape = Season)) +
        geom_point(size = 3) +
        facet_wrap(~ Depth) +
        theme_light() +
        labs(
                x = "PCoA1",
                y = "PCoA2",
                color = "Wetland type",
                shape = "Season"
        )

## PCoA – Bacteria
bact_pcoa <- cmdscale(bact_dist, eig = TRUE, k = 2)

scores_bact <- as.data.frame(bact_pcoa$points)
colnames(scores_bact) <- c("PCoA1", "PCoA2")
scores_bact$Sample <- rownames(scores_bact)

scores_bact <- scores_bact %>%
        left_join(metadata, by = "Sample") %>%
        mutate(Domain = "Bacteria")

ggplot(scores_bact, aes(x = PCoA1, y = PCoA2,
                        color = Wetland,
                        shape = Season)) +
        geom_point(size = 3) +
        facet_wrap(~ Depth) +
        theme_light() +
        labs(
                x = "PCoA1",
                y = "PCoA2",
                color = "Wetland type",
                shape = "Season"
        )

scores_all <- bind_rows(scores_arch, scores_bact)

ggplot(scores_all, aes(x = PCoA1, y = PCoA2,
                       color = Wetland,
                       shape = Season)) +
        geom_point(size = 3) +
        facet_grid(Domain ~ Depth) +
        theme_light() +
        labs(
                x = "PCoA1",
                y = "PCoA2",
                color = "Wetland type",
                shape = "Season"
        )

arch_var <- arch_pcoa$eig / sum(arch_pcoa$eig) * 100
bact_var <- bact_pcoa$eig / sum(bact_pcoa$eig) * 100


ggplot(scores_arch, aes(PCoA1, PCoA2,
                        color = Wetland,
                        shape = Season)) +
        geom_point(size = 3) +
        facet_wrap(~ Depth, ncol = 1) +
        theme_light() +
        labs(
                x = paste0("PCoA1 (", round(arch_var[1], 1), "%)"),
                y = paste0("PCoA2 (", round(arch_var[2], 1), "%)"),
                color = "Wetland type",
                shape = "Season"
        )

ggplot(scores_bact, aes(PCoA1, PCoA2,
                        color = Wetland,
                        shape = Season)) +
        geom_point(size = 3) +
        facet_wrap(~ Depth, ncol = 1) + 
        theme_light() +
        labs(
                x = paste0("PCoA1 (", round(bact_var[1], 1), "%)"),
                y = paste0("PCoA2 (", round(bact_var[2], 1), "%)"),
                color = "Wetland type",
                shape = "Season"
        )


#########################################################
#PCA

library(ggplot2)
library(vegan)
library(dplyr)
library(ggrepel)

amostras <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", 
              "S10", "S11", "S12", "S13", "S15", "S16", "S17", "S18", "S20", 
              "S21", "S22", "S23", "S24")

Wetland <- c(rep("Permanent", 6), rep("Seasonal", 6), rep("Permanent", 5), rep("Seasonal", 5))

Depth <- c("0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "0-10 cm", "10-20 cm", "0-10 cm", "0-10 cm", "10-20 cm", 
           "0-10 cm","10-20 cm", "10-20 cm", "0-10 cm", "10-20 cm",
           "0-10 cm", "10-20 cm")

Season <- c(rep("Rainy", 12),
            rep("Dry", 10))

metadata <- data.frame(
        Sample  = amostras,
        Wetland = Wetland,
        Depth   = Depth,
        Season  = Season
)


arch <- read_csv2("Class_archeae.csv")
bact <- read_csv2("Class_bacteria.csv")

arch <- arch %>% filter(Class != "Unclassified") 
bact <- bact %>% filter(Class != "Unclassified")

sample_names <- metadata$Sample

# Archaea
arch_samples <- arch %>%
        select(Class, all_of(sample_names))

arch_mat <- as.matrix(arch_samples[ , -1])        
rownames(arch_mat) <- arch_samples$Class           

# Bacteria
bact_samples <- bact %>%
        select(Class, all_of(sample_names))

bact_mat <- as.matrix(bact_samples[ , -1])
rownames(bact_mat) <- bact_samples$Class

arch_mat_t <- t(arch_mat) 
bact_mat_t <- t(bact_mat) 

all(rownames(arch_mat_t) == metadata$Sample)
all(rownames(bact_mat_t) == metadata$Sample)

## 3. Pseudocount + CLR

arch_clr <- clr(arch_mat_t + 1)
bact_clr <- clr(bact_mat_t + 1)

pca_arch <- prcomp(arch_mat_t, center = TRUE, scale. = FALSE)

scores_arch <- as.data.frame(pca_arch$x[,1:2])
scores_arch$Sample <- rownames(scores_arch)
scores_arch <- left_join(scores_arch, metadata, by = "Sample")

load_arch <- as.data.frame(pca_arch$rotation[,1:2])
load_arch$Class <- rownames(load_arch)

mult <- min(
        (max(scores_arch$PC1) - min(scores_arch$PC1)) / (max(load_arch$PC1) - min(load_arch$PC1)),
        (max(scores_arch$PC2) - min(scores_arch$PC2)) / (max(load_arch$PC2) - min(load_arch$PC2))
)

load_arch$PC1 <- load_arch$PC1 * mult * 0.7
load_arch$PC2 <- load_arch$PC2 * mult * 0.7


# Biplot
pA<-ggplot() +
        geom_point(data = scores_arch,
                   aes(PC1, PC2, color = Wetland, shape = Season),
                   size = 3) +
        geom_segment(data = load_arch,
                     aes(x = 0.0, y = 0.0, xend = PC1, yend = PC2),
                     color = "gray60") +
        geom_text_repel(data = load_arch,
                        aes(PC1, PC2, label = Class),
                        size = 3,
                        max.overlaps = Inf,
                        box.padding = 0.3,
                        point.padding = 0.2,
                        force = 2,           
                        force_pull = 1) +
        labs(
                x = paste0("PC1 (", round(summary(pca_arch)$importance[2,1] * 100, 1), "%)"),
                y = paste0("PC2 (", round(summary(pca_arch)$importance[2,2] * 100, 1), "%)"),
                color = "Wetland type",
                shape = "Season"
        ) +
        theme_light()

pA <- pA + guides(color = "none", shape = "none")

pca_bact <- prcomp(bact_mat_t, center = T, scale. = F)

scores_bact <- as.data.frame(pca_bact$x[,1:2])
scores_bact$Sample <- rownames(scores_bact)
scores_bact <- left_join(scores_bact, metadata, by = "Sample")

load_bact <- as.data.frame(pca_bact$rotation[,1:2])
load_bact$Class <- rownames(load_bact)

mult_b <- min(
        (max(scores_bact$PC1) - min(scores_bact$PC1)) / (max(load_bact$PC1) - min(load_bact$PC1)),
        (max(scores_bact$PC2) - min(scores_bact$PC2)) / (max(load_bact$PC2) - min(load_bact$PC2))
)


load_bact_plot <- load_bact %>%
        mutate(strength = sqrt(PC1^2 + PC2^2)) %>%
        arrange(desc(strength)) %>%
        slice(1:15)

Lmax <- quantile(load_bact_plot$strength, 0.5)

load_bact_plot <- load_bact_plot %>%
        mutate(
                scale_factor = ifelse(strength > Lmax, Lmax / strength, 1),
                PC1_plot = PC1 * scale_factor,
                PC2_plot = PC2 * scale_factor
        )


pB <- ggplot() +
        geom_point(
                data = scores_bact,
                aes(PC1, PC2, color = Wetland, shape = Season),
                size = 3
        ) +
        geom_segment(
                data = load_bact_plot,
                aes(x = 0, y = 0, xend = PC1_plot, yend = PC2_plot),
                color = "gray60"
        ) +
        geom_text_repel(
                data = load_bact_plot,
                aes(PC1_plot, PC2_plot, label = Class),
                size = 3,
                max.overlaps = Inf,
                segment.color = NA
        ) +
        theme_light() +
        labs(
                x = paste0("PC1 (", round(summary(pca_bact)$importance[2,1] * 100, 1), "%)"),
                y = paste0("PC2 (", round(summary(pca_bact)$importance[2,2] * 100, 1), "%)"),
                color = "Wetland type",
                shape = "Season"
        )

pB 

fig_pca <- pA + pB +
        plot_annotation(tag_levels = "A")

ggsave("Figure_PCA_AB.png",
        fig_pca,
        width = 16,
        height = 6,
        dpi = 300)
