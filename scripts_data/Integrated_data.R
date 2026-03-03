setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data")
df <- read.csv2("tabela_mestre.csv", check.names = FALSE)

library(tidyverse)
library(writexl)
library(broom)
library(agricolae)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(FSA)

df<-df %>% rename(
        "TCsoil" = CT,
        "Soil moisture" = umidade,
        "pH" = pH,
        "Eh" = Eh,
        "NO₃⁻" = Nitrato,
        "SO₄²⁻" = Sulfato,
        "Fe²⁺" = `Fe(II)`,
        "Fe³⁺" = `Fe(III)`,
        "CO₂" = CO2,
        "CH₄" = CH4,
        "SUVA₂₅₄" = SUVA,
        "E₂/E₃" = E2_E3,
        "E₄/E₅" = E4_E5,
        "E₄/E₆" = E4_E6,
        "SR" = SR,
        "Carboxylic acids" = carboxilicos,
        "Phenolic hydroxyls" = fenolicos,
        "Total acidity" = total,
        "Functional aromaticity" = aromaticidade_func,
        "Spectroscopic aromaticity" = aromaticidade_esp,
        "ΔABS(nat-red)" = EAC
)

vars_quim <- df %>%
        dplyr::select(
                TCsoil,
                `Soil moisture`,
                pH,
                Eh,
                `NO₃⁻`,
                `SO₄²⁻`,
                `Fe²⁺`,
                `Fe³⁺`,
                `CO₂`,
                `CH₄`,
                `SUVA₂₅₄`,
                `E₂/E₃`,
                `E₄/E₅`,
                `E₄/E₆`,
                SR,
                `Carboxylic acids`,
                `Phenolic hydroxyls`,
                `Total acidity`,
                `Functional aromaticity`,
                `Spectroscopic aromaticity`,
                `ΔABS(nat-red)`
        )

vars_quim <- names(vars_quim)

#Normality
normality_results <- df %>%
        pivot_longer(
                cols = all_of(vars_quim),
                names_to = "variable",
                values_to = "value"
        ) %>%
        mutate(value = as.numeric(value)) %>%
        group_by(matriz, tipo, variable) %>%
        summarise(
                p_value = ifelse(
                        sum(!is.na(value)) >= 3 & sd(value, na.rm = TRUE) > 0,
                        shapiro.test(value[!is.na(value)])$p.value,
                        NA_real_
                ),
                normal = ifelse(!is.na(p_value) & p_value > 0.05, "Sim", "Não"),
                .groups = "drop"
        )

#write_xlsx(normality_results, "resultados_normalidade.xlsx")


# Kruskal-Wallis+Dunn


vars_solo <- df %>%
        dplyr::select(
                TCsoil,
                `Soil moisture`,
                pH,
                Eh,
                `NO₃⁻`,
                `SO₄²⁻`,
                `Fe²⁺`,
                `Fe³⁺`,
                `CO₂`,
                `CH₄`
        )

vars_solo <- names(vars_solo)

vars_humicas <- df %>%
        dplyr::select(`SUVA₂₅₄`,
                      `E₂/E₃`,
                      `E₄/E₅`,
                      `E₄/E₆`,
                      SR,
                      `Carboxylic acids`,
                      `Phenolic hydroxyls`,
                      `Total acidity`,
                      `Functional aromaticity`,
                      `Spectroscopic aromaticity`,
                      `ΔABS(nat-red)`
        )

vars_humicas <- names(vars_humicas)

kruskal_lista <- list()
dunn_lista <- list()

# =========================================================
#  SOLO (sem matriz)
# =========================================================

df_solo <- df %>%
        mutate(grupo = interaction(tipo, estacao, drop = TRUE))

for(v in vars_solo){
        
        formula_str <- as.formula(paste0("`", v, "` ~ grupo"))
        
        kw <- tryCatch(
                kruskal.test(formula_str, data = df_solo),
                error = function(e) NULL
        )
        
        if(!is.null(kw)){
                
                kruskal_lista[[v]] <- data.frame(
                        Variavel = v,
                        Grupo = "solo",
                        Chi_square = kw$statistic,
                        df = kw$parameter,
                        p_value = kw$p.value
                )
                
                if(kw$p.value < 0.05){
                        
                        dunn <- FSA::dunnTest(formula_str,
                                              data = df_solo,
                                              method = "bh")
                        
                        tmp <- dunn$res
                        tmp$Variavel <- v
                        tmp$Grupo <- "solo"
                        
                        dunn_lista[[paste0(v,"_solo")]] <- tmp
                }
        }
}

# Humic fractions

df_hum <- df %>%
        mutate(grupo = interaction(tipo, estacao, matriz, drop = TRUE))

for(v in vars_humicas){
        
        formula_str <- as.formula(paste0("`", v, "` ~ grupo"))
        
        kw <- tryCatch(
                kruskal.test(formula_str, data = df_hum),
                error = function(e) NULL
        )
        
        if(!is.null(kw)){
                
                kruskal_lista[[v]] <- data.frame(
                        Variavel = v,
                        Grupo = "humicas",
                        Chi_square = kw$statistic,
                        df = kw$parameter,
                        p_value = kw$p.value
                )
                
                if(kw$p.value < 0.05){
                        
                        dunn <- FSA::dunnTest(formula_str,
                                              data = df_hum,
                                              method = "bh")
                        
                        tmp <- dunn$res
                        tmp$Variavel <- v
                        tmp$Grupo <- "humicas"
                        
                        dunn_lista[[paste0(v,"_hum")]] <- tmp
                }
        }
}


final_kw <- bind_rows(kruskal_lista)
final_dunn <- bind_rows(dunn_lista)

write_xlsx(
        list(
                Kruskal = final_kw,
                Dunn = final_dunn
        ),
        "Resultados_Kruskal_Dunn.xlsx"
)

# Table mean +sd


resumo_media_dp <- function(data, variaveis, grupos){
        
        data %>%
                pivot_longer(
                        cols = all_of(variaveis),
                        names_to = "Variavel",
                        values_to = "Valor"
                ) %>%
                group_by(across(all_of(grupos)), Variavel) %>%
                summarise(
                        N = sum(!is.na(Valor)),
                        Media = mean(Valor, na.rm = TRUE),
                        DP = sd(Valor, na.rm = TRUE),
                        Resultado = sprintf("%.2f ± %.2f", Media, DP),
                        .groups = "drop"
                )
}


#Soi: season and wetland type


tabela_solo <- resumo_media_dp(
        data = df,
        variaveis = vars_solo,
        grupos = c("tipo", "estacao")
)


# =========================================================
# HÚMICAS — estação, tipo e matriz
# =========================================================

tabela_humicas <- resumo_media_dp(
        data = df,
        variaveis = vars_humicas,
        grupos = c("tipo", "estacao", "matriz")
)


write_xlsx(
        list(
                Solo_resumo = tabela_solo,
                Solo_estacao = compar_estacao_solo,
                Solo_tipo = compar_tipo_solo,
                Humicas_resumo = tabela_humicas,
                Humicas_estacao = compar_estacao_hum,
                Humicas_tipo = compar_tipo_hum
        ),
        "Tabela_Media_DP.xlsx"
)


# CORRELAÇÕES

library(tidyverse)
library(ggcorrplot)


df_solo_cor <- df %>%
        select(all_of(vars_solo)) %>%
        mutate(across(everything(), as.numeric)) %>%
        drop_na()


cor_solo <- cor(df_solo_cor, method = "spearman")


library(tidyverse)
library(ggcorrplot)
library(Hmisc)


df_solo_cor <- df %>%
        select(all_of(vars_solo)) %>%
        mutate(across(everything(), as.numeric)) %>%
        drop_na()


cor_res <- rcorr(as.matrix(df_solo_cor), type = "spearman")

cor_mat <- cor_res$r  
p_mat   <- cor_res$P  

p_cor_solo <- ggcorrplot(
        cor_mat,
        type = "lower",
        lab = TRUE,
        p.mat = p_mat,
        sig.level = 0.05,
        insig = "blank",      #esconde não-significativos
        colors = c("#67001F", "white", "#053061")) +
        theme_light(base_family = "Helvetica")+
        theme(plot.title = element_blank(),
              axis.title = element_blank())

plot(p_cor_solo)

vars_hum_gee <- c(vars_humicas, "CO₂", "CH₄")

df_hum_cor <- df %>%
        select(
                matriz,
                all_of(vars_humicas) ) %>%
        mutate(across(
                c(all_of(vars_humicas)),
                as.numeric
        ))

plot_cor_hum <- function(data, nome){
        
        mat <- data %>%
                select(all_of(vars_humicas)) %>%
                mutate(across(everything(), as.numeric)) %>%
                drop_na()
        
        cor_mat <- cor(mat, method = "spearman")
        
        p_mat <- ggcorrplot::cor_pmat(mat, method = "spearman")
        
        p <- ggcorrplot(
                cor_mat,
                type = "lower",
                lab = TRUE,
                p.mat = p_mat,
                sig.level = 0.05,
                insig = "blank",      # esconde não-significativos
                colors = c("#2166AC", "white", "#B2182B")
        ) +
                theme_bw(base_family = "Helvetica") +
                theme(
                        plot.title = element_blank(),
                        axis.title = element_blank(),
                        axis.text.x = element_text(
                                angle = 90,
                                hjust = 1,
                                vjust = 0.5
                        )
                )
        
        return(p)
}

df_AF <- df_hum_cor %>% filter(matriz == "FA")

p_cor_AF <- plot_cor_hum(df_AF, "FA")
p_cor_AF

df_SHT <- df_hum_cor %>% dplyr::filter(matriz == "THS")

p_cor_SHT <- plot_cor_hum(df_SHT, "THS")
p_cor_SHT

cor_plot <- df %>%
        ggplot(aes(x =`Functional aromaticity`, y =`ΔABS(nat-red)`)) +
        geom_point(aes(color = tipo, shape = estacao, group = replicata), size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE, fill = "grey80") +
        ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
        
        labs(   x = "Aromaticidade Funcional",
                y = "Δ ABS (Capacidade de aceitação de elétrons)",
                color = "Matriz",
                shape = "Estação"
        ) +
        
        theme_bw(base_family = "Helvetica") +
        theme(
                panel.grid.minor = element_blank(),
                legend.position = "right"
        )+
        facet_grid(~matriz)

print(cor_plot)

# PCA

library(factoextra)

dados_AF <- df[1:72,-c(4)]
dados_AH <- df[73:144,-c(4)]

dados_pca_AF <- dados_AF[, sapply(dados_AF, is.numeric)]
dados_pca_AH <- dados_AH[, sapply(dados_AH, is.numeric)]

dados_pca_AF <- dados_pca_AF[, apply(dados_pca_AF, 2, var, na.rm=TRUE) > 0]
dados_pca_AH <- dados_pca_AH[, apply(dados_pca_AH, 2, var, na.rm=TRUE) > 0]

res.pca_AF <- prcomp(dados_pca_AF, scale. = TRUE)
res.pca_AH <- prcomp(dados_pca_AH, scale. = TRUE)

contrib_AF<-fviz_contrib(res.pca_AF, choice = "var", axes = 1:2)+
        labs(title = NULL) +
        theme_light(base_family = "Helvetica")+
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_text(
                      angle = 90,
                      hjust = 1,
                      vjust = 0.5
              ))

contrib_AH<-fviz_contrib(res.pca_AH, choice = "var", axes = 1:2)+
        labs(title = NULL) +
        theme_light(base_family = "Helvetica")+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(
                      angle = 90,
                      hjust = 1,
                      vjust = 0.5
              ))

plot(contrib_AF)
plot(contrib_AH)

dados_AF <- df[1:72,-c(4,8:10,12:14,16:18,20:22,24,27,26,28)]
dados_AH <- df[73:144,-c(4,8,10:12, 13, 14,16:18,20:22,24,27,26,28)]

dados_pca_AF <- dados_AF[, sapply(dados_AF, is.numeric)]
dados_pca_AH <- dados_AH[, sapply(dados_AH, is.numeric)]

dados_pca_AF <- dados_pca_AF[, apply(dados_pca_AF, 2, var, na.rm=TRUE) > 0]
dados_pca_AH <- dados_pca_AH[, apply(dados_pca_AH, 2, var, na.rm=TRUE) > 0]

res.pca_AF <- prcomp(dados_pca_AF, scale. = TRUE)
res.pca_AH <- prcomp(dados_pca_AH, scale. = TRUE)

pca_AF<-fviz_pca_biplot(res.pca_AF,
                        habillage = dados_AF$tipo, 
                        addEllipses = F,
                        col.var = "grey40",
                        label = "var",
                        repel = T)+
        labs(title = NULL) +
        theme_light(base_family = "Helvetica")+
        theme(
                legend.position = "bottom"
        )

#fviz_contrib(res.pca, choice = "var", axes = 1:2)
pca_AH<-fviz_pca_biplot(res.pca_AH,
                        habillage = dados_AH$tipo, 
                        addEllipses = F,
                        col.var = "grey40",
                        label = "var",
                        repel = T)+
        labs(title = NULL) +
        theme_light(base_family = "Helvetica")+
        theme(
                legend.position = "bottom"
        )

pca_AH
pca_AF

library(patchwork)

pca_AF_final <- (contrib_AF|pca_AF) +
        plot_layout(widths = c(1, 1.5)) +
        plot_annotation(tag_levels = "A")

pca_AH_final <-(contrib_AH|pca_AH) +
        plot_layout(widths = c(1, 1.5)) +
        plot_annotation(tag_levels = "A") 

ggsave("pca_AF_final.png", width = 10, height = 6)
ggsave("pca_AH_final.png", width = 10, height = 6)

# Variação temporal

dt<-df %>% rename(
        "TCsolo" = CT,
        "Umidade do solo" = umidade,
        "pH" = pH,
        "Eh" = Eh,
        "NO₃⁻" = Nitrato,
        "SO₄²⁻" = Sulfato,
        "Fe²⁺" = `Fe(II)`,
        "Fe³⁺" = `Fe(III)`,
        "CO₂" = CO2,
        "CH₄" = CH4,
        "SUVA₂₅₄" = SUVA,
        "E₂/E₃" = E2_E3,
        "E₄/E₅" = E4_E5,
        "E₄/E₆" = E4_E6,
        "SR" = SR,
        "Ácidos carboxílicos" = carboxilicos,
        "Grupos fenólicos" = fenolicos,
        "Acidez total" = total,
        "Aromaticidade funcional" = aromaticidade_func,
        "Aromaticidade espectroscópica" = aromaticidade_esp,
        "ΔABS(nat-red)" = EAC
)

vars_quim_t <- dt %>%
        dplyr::select(
                TCsolo,
                `Umidade do solo`,
                pH,
                Eh,
                `NO₃⁻`,
                `SO₄²⁻`,
                `Fe²⁺`,
                `Fe³⁺`,
                `CO₂`,
                `CH₄`,
                `SUVA₂₅₄`,
                `E₂/E₃`,
                `E₄/E₅`,
                `E₄/E₆`,
                SR,
                `Ácidos carboxílicos`,
                `Grupos fenólicos`,
                `Acidez total`,
                `Aromaticidade funcional`,
                `Aromaticidade espectroscópica`,
                `ΔABS(nat-red)`
        )

vars_quim_t <- names(vars_quim_t)

temp_plot <- dt %>%
        pivot_longer(
                cols = all_of(vars_quim_t),
                names_to = "variable",
                values_to = "value"
        ) %>%
        mutate(value = as.numeric(value)) %>%
        ggplot(aes(as.factor(coleta), value,
                   shape=estacao, color=matriz, group = matriz)) +
        stat_summary(fun=mean, geom="line") +
        stat_summary(fun=mean, geom="point") +
        facet_wrap(
                tipo~variable,
                scales="free_y") +
        theme_bw()

ggsave("variacao_temporal_cientifico.png",
       temp_plot,
       width=16,height=12,dpi=300)



library(tidyverse)

ordem_vars_t <- c(
        TCsolo,
        `Umidade do solo`,
        pH,
        Eh,
        `NO₃⁻`,
        `SO₄²⁻`,
        `Fe²⁺`,
        `Fe³⁺`,
        `CO₂`,
        `CH₄`,
        `SUVA₂₅₄`,
        `E₂/E₃`,
        `E₄/E₅`,
        `E₄/E₆`,
        SR,
        `Ácidos carboxílicos`,
        `Grupos fenólicos`,
        `Acidez total`,
        `Aromaticidade funcional`,
        `Aromaticidade espectroscópica`,
        `ΔABS(nat-red)`
)

# Ordem
ordem_vars <- c(
        "TCsoil",
        "Soil moisture",
        "pH",
        "Eh",
        "NO₃⁻",
        "SO₄²⁻",
        "Fe²⁺",
        "Fe³⁺",
        "CO₂",
        "CH₄",
        "SUVA₂₅₄",
        "E₂/E₃",
        "E₄/E₅",
        "E₄/E₆",
        "SR",
        "Carboxylic acids",
        "Phenolic hydroxyls",
        "Total acidity",
        "Functional aromaticity",
        "Spectroscopic aromaticity",
        "ΔABS(nat-red)"
)

temp_plot <- dt %>%
        
        pivot_longer(
                cols = all_of(vars_quim_t),
                names_to = "variable",
                values_to = "value"
        ) %>%
        
        mutate(
                value = as.numeric(value),
                variable = factor(variable, levels = ordem_vars_t),
                coleta = as.factor(coleta)
        ) %>%
        
        ggplot(
                aes(
                        x = coleta,
                        y = value,
                        color = matriz,
                        shape = estacao,
                        group = matriz
                )
        ) +
        
        stat_summary(
                fun = mean,
                geom = "line",
                linewidth = 0.9
        ) +
        
        stat_summary(
                fun = mean,
                geom = "point",
                size = 2.2
        ) +
        
        facet_wrap(
                tipo~variable,
                scales = "free_y",
                nrow = 4
        ) +
        labs(
                x = "Coleta",
                y = NULL,
                color = "Matrix",
                shape = "Season"
        ) +
        
        theme_bw(base_family = "Helvetica") +
        
        theme(
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(size = 9),
                
                axis.text.x = element_text(
                        angle = 90,
                        hjust = 1,
                        vjust = 0.5
                ),
                
                legend.position = "bottom",
                legend.box = "horizontal"
        )

temp_plot

ggsave(
        "variacao_temporal_cientifico.png",
        temp_plot,
        width = 18,
        height = 14,
        dpi = 600
)