setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data")
dados.depth <-read.csv2("Dados_16S_depth.csv", stringsAsFactors = T, na.strings = "NA")

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(cowplot)

dados.depth <- dados.depth %>%
        rename(`Soil moisture` = Humidity)

vars_env <- c("pH", "Eh", "Soil moisture", "Csoil", "N", "P", "K")
vars_taxa <- c(
        "CandidatusNitrososphaera","Methanobacterium","Methanocella","Methanomassiliicoccus","Methermicoccus","Nitrosotalea",
        "Pedosphaerae","Saprospirae","Spartobacteria","Acidobacteria","Acidobacteriia","Actinobacteria","Alphaproteobacteria",
        "Anaerolineae","Bacilli","Bacteroidia","Betaproteobacteria","Clostridia","DA052","Deltaproteobacteria","Ellin6529",
        "Flavobacteriia","Gammaproteobacteria","Gemmatimonadetes","Ktedonobacteria","MB.A2.108","Nitrospira","Planctomycetia",
        "SJA.176","Solibacteres","Thermoleophilia","TK10","TK17","Zetaproteobacteria"
)

archaea_taxa <- c("CandidatusNitrososphaera","Methanobacterium","Methanocella",
                  "Methanomassiliicoccus","Methermicoccus","Nitrosotalea")
bacteria_taxa <- setdiff(vars_taxa, archaea_taxa)


if (!exists("dados.depth")) stop("Objeto 'dados.depth' não existe. Carrega o CSV primeiro.")
missing_taxa <- setdiff(vars_taxa, colnames(dados.depth))
if (length(missing_taxa) > 0) stop("Faltam colunas de taxa em 'dados.depth': ", paste(missing_taxa, collapse = ", "))
env_cols <- intersect(vars_env, colnames(dados.depth))
if (length(env_cols) == 0) stop("Nenhuma das variáveis ambientais listadas existe em 'dados.depth'.")
meta_cols <- c("Wetland.Type","Season","Depth","Sampling.point") %>% intersect(colnames(dados.depth))

taxa_df <- dados.depth %>% select(all_of(vars_taxa))
taxa_df[is.na(taxa_df)] <- NA

# CLR by group
clr_by_group <- function(df_group, pseudocount = 1e-6) {
        mat <- as.matrix(df_group)
        mat[mat == 0 & !is.na(mat)] <- pseudocount
        clr_mat <- t(apply(mat, 1, function(x) {
                if (all(is.na(x))) return(rep(NA_real_, length(x)))
                lx <- log(x)
                lx - mean(lx, na.rm = TRUE)
        }))
        colnames(clr_mat) <- colnames(mat)
        as.data.frame(clr_mat)
}

taxa_archaea_df <- taxa_df %>% select(all_of(intersect(archaea_taxa, colnames(taxa_df))))
taxa_bacteria_df <- taxa_df %>% select(all_of(intersect(bacteria_taxa, colnames(taxa_df))))

archaea_clr_df <- if (ncol(taxa_archaea_df) > 0) clr_by_group(taxa_archaea_df) else NULL
bacteria_clr_df <- if (ncol(taxa_bacteria_df) > 0) clr_by_group(taxa_bacteria_df) else NULL

n_samples <- nrow(dados.depth)
taxa_clr_df <- as.data.frame(matrix(NA_real_, nrow = n_samples, ncol = length(vars_taxa)))
colnames(taxa_clr_df) <- vars_taxa
rownames(taxa_clr_df) <- rownames(dados.depth)

if (!is.null(archaea_clr_df)) taxa_clr_df[colnames(archaea_clr_df)] <- archaea_clr_df
if (!is.null(bacteria_clr_df)) taxa_clr_df[colnames(bacteria_clr_df)] <- bacteria_clr_df

clr_data <- bind_cols(
        dados.depth %>% select(all_of(meta_cols)),
        dados.depth %>% select(all_of(env_cols)),
        taxa_clr_df
)

clr_data <- clr_data %>%
        mutate(Depth = recode(Depth,
                              "0-10" = "0-10 cm",
                              "10-20" = "10-20 cm"))

#heatmaps
domain_map <- tibble(Taxon = vars_taxa,
                     Domain = ifelse(Taxon %in% archaea_taxa, "Archaea", "Bacteria"))

safe_cor_test <- function(x, y) {
        ok <- complete.cases(x, y)
        x <- x[ok]; y <- y[ok]
        if (length(x) < 3 || var(x, na.rm = TRUE) == 0 || var(y, na.rm = TRUE) == 0) {
                return(list(r = NA_real_, p = NA_real_))
        }
        tt <- cor.test(x, y, method = "pearson")
        list(r = as.numeric(tt$estimate), p = tt$p.value)
}

wetland_types <- unique(clr_data$Wetland.Type)
depth_levels <- unique(clr_data$Depth)

cor_results <- map_dfr(wetland_types, function(w) {
        map_dfr(depth_levels, function(d) {
                df_sub <- clr_data %>% filter(Wetland.Type == w, Depth == d)
                combos <- expand.grid(Taxon = vars_taxa, Variable = env_cols, stringsAsFactors = FALSE) %>%
                        mutate(Wetland.Type = w, Depth = d)
                if (nrow(df_sub) < 3) {
                        combos$R <- NA_real_; combos$P <- NA_real_
                        return(combos %>% select(Wetland.Type, Depth, Taxon, Variable, R, P))
                }
                combos <- combos %>%
                        mutate(
                                R = map2_dbl(Taxon, Variable, ~ safe_cor_test(df_sub[[.x]], df_sub[[.y]])$r),
                                P = map2_dbl(Taxon, Variable, ~ safe_cor_test(df_sub[[.x]], df_sub[[.y]])$p)
                        ) %>% select(Wetland.Type, Depth, Taxon, Variable, R, P)
                combos
        })
})

cor_results <- cor_results %>%
        left_join(domain_map, by = "Taxon") %>%
        mutate(Significance = case_when(
                !is.na(P) & P <= 0.001 ~ "***",
                !is.na(P) & P <= 0.01  ~ "**",
                !is.na(P) & P <= 0.05  ~ "*",
                TRUE ~ ""
        ))

taxon_order_df <- cor_results %>%
        group_by(Domain, Taxon) %>%
        summarize(mean_abs_R = mean(abs(R), na.rm = TRUE), .groups = "drop") %>%
        arrange(Domain, desc(mean_abs_R))

taxon_levels_bac <- taxon_order_df %>% filter(Domain == "Bacteria") %>% pull(Taxon)
taxon_levels_arch <- taxon_order_df %>% filter(Domain == "Archaea") %>% pull(Taxon)

cor_results <- cor_results %>%
        mutate(
                Taxon = factor(Taxon, levels = c(rev(taxon_levels_bac), rev(taxon_levels_arch))),
                Variable = factor(Variable, levels = env_cols),
                Depth = factor(Depth, levels = depth_levels),
                Wetland.Type = factor(Wetland.Type, levels = wetland_types)
        )

cor_results <- cor_results %>%
        mutate(
                Taxon_italics = gsub("_", " ", Taxon),  # só garante que não haja underscores
                Taxon_italics = paste0("italic('", Taxon_italics, "')")
        )

plot_domain_heatmap <- function(df, domain_name, title = domain_name) {
        df_sub <- df %>% filter(Domain == domain_name)
        ggplot(df_sub, aes(x = Variable, y = Taxon_italics, fill = R)) +
                geom_tile(color = "white", na.rm = TRUE) +
                geom_text(aes(label = Significance), size = 3, na.rm = TRUE) +
                scale_fill_gradient2(low = "orchid2", mid = "turquoise", high = "royalblue",
                                     midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
                facet_grid(Depth ~ Wetland.Type, scales = "free_y", space = "free_y") +
                labs(x = "Environmental variables", y = paste0("Taxa (CLR)"),
                     title = paste0(title)) +
                scale_y_discrete(labels = function(x) parse(text = x))+
                theme_minimal(base_family = "Helvetica") +
                theme(
                        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        strip.text = element_text(face = "bold"),
                        panel.grid = element_blank()
                )
}

plot_domain_heatmap(cor_results, "Bacteria", "Bacteria")
p_arch<-plot_domain_heatmap(cor_results, "Archaea", "Archaea")

ggsave("heatmap_archaea.png", plot = p_arch, width = 10, height = 6, dpi = 300)

###################
cor_results <- cor_results %>%
        mutate(
                R2 = R^2,
                R2_label = ifelse(!is.na(R2), sprintf("%.2f", R2), ""),
                Label = paste0(R2_label, Significance)
        )

plot_domain_heatmap <- function(df, domain_name, title = domain_name) {
        df_sub <- df %>% filter(Domain == domain_name)
        ggplot(df_sub, aes(x = Variable, y = Taxon_italics, fill = R)) +
                geom_tile(color = "white", na.rm = TRUE) +
                geom_text(aes(label = Label), size = 3, na.rm = TRUE)+
                scale_fill_gradient2(low = "orchid2", mid = "turquoise", high = "royalblue",
                                     midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
                facet_grid(Depth ~ Wetland.Type, scales = "free_y", space = "free_y") +
                labs(x = "Environmental variables", y = paste0("Taxa (CLR)"),
                     title = paste0(title)) +
                scale_y_discrete(labels = function(x) parse(text = x))+
                theme_minimal(base_family = "Helvetica") +
                theme(
                        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        strip.text = element_text(face = "bold"),
                        panel.grid = element_blank()
                )
}

plot_domain_heatmap(cor_results, "Bacteria", "Bacteria")
plot_domain_heatmap(cor_results, "Archaea", "Archaea")

