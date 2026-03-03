setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data")

dados_quimicos <- read.csv2("Grupos.csv", stringsAsFactors = FALSE)

library(tidyverse)
library(rstatix)

dados_quimicos <- dados_quimicos %>%
        mutate(
                matriz = factor(matriz),
                tipo = factor(tipo),
                estacao = factor(estacao)
        )

dados_quimicos_long <- dados_quimicos %>%
        pivot_longer(
                cols = c(carboxilicos, fenolicos, total, Aromaticidade),
                names_to = "variavel",
                values_to = "valor"
        )

normalidade <- dados_quimicos_long %>%
        group_by(variavel) %>%
        shapiro_test(valor)

normalidade

dados_quimicos_long %>%
        group_by(variavel, estacao) %>%
        shapiro_test(valor)

#####################################
dados_quimicos_long %>%
        group_by(variavel) %>%
        kruskal_test(valor ~ estacao)

dados_quimicos_long %>%
        group_by(variavel) %>%
        kruskal_test(valor ~ tipo)

dados_quimicos_long %>%
        group_by(variavel) %>%
        kruskal_test(valor ~ matriz)

dados_quimicos_long %>%
        group_by(variavel) %>%
        dunn_test(valor ~ estacao, p.adjust.method = "bonferroni")

dados_quimicos_long %>%
        group_by(variavel, estacao) %>%
        summarise(
                media = mean(valor),
                sd = sd(valor),
                mediana = median(valor),
                .groups = "drop"
        )

library(ggpubr)

ggplot(dados_quimicos_long, aes(x = estacao, y = valor)) +
        geom_boxplot() +
        facet_wrap(~ variavel, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Estação",
                y = "Valor",
                title = "Comparação por estação"
        )

ggplot(dados_quimicos_long, aes(x = tipo, y = valor)) +
        geom_boxplot() +
        facet_wrap(~ variavel+matriz, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Tipo",
                y = "Valor",
                title = "Comparação por tipo de fração húmica"
        )
ggplot(dados_quimicos_long, aes(x = matriz, y = valor)) +
        geom_boxplot() +
        facet_wrap(~ variavel, scales = "free_y") +
        theme_bw() +
        labs(
                x = "Matriz",
                y = "Valor",
                title = "Comparação por matriz"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(dados_quimicos_long,
       aes(x = coleta, y = valor, group = interaction(matriz, tipo), color = variavel)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        facet_wrap(tipo~ matriz) +
        theme_bw() +
        labs(
                x = "Tempo (coletas)",
                y = "Valor",
                color = "Fração",
                title = "Evolução temporal dos grupos funcionais"
        )

#Build table with all indexes

library(dplyr)

funcionais_med <- dados_quimicos %>%
        group_by(amostra, coleta, matriz) %>% 
        summarise(
                carboxilicos_media = mean(carboxilicos, na.rm = TRUE),
                carboxilicos_sd    = sd(carboxilicos, na.rm = TRUE),
                fenolicos_media    = mean(fenolicos, na.rm = TRUE),
                fenolicos_sd       = sd(fenolicos, na.rm = TRUE),
                total_media        = mean(total, na.rm = TRUE),
                total_sd           = sd(total, na.rm = TRUE),
                aromaticidade_media = mean(Aromaticidade, na.rm = TRUE),
                aromaticidade_sd    = sd(Aromaticidade, na.rm = TRUE),
                .groups = "drop"
        )

funcionais_med

funcionais_med <- funcionais_med %>%
        rename(
                Amostra = amostra,
                Coleta  = coleta,
                Matriz  = matriz
        )

tabela_final_indices <- tabela_final_indices %>%
        mutate(Coleta = as.character(Coleta))

funcionais_med <- funcionais_med %>%
        mutate(Coleta = as.character(Coleta))


tabela_final_completa <- tabela_final_indices %>%
        left_join(
                funcionais_med,
                by = c("Amostra", "Coleta", "Matriz")
        )
glimpse(tabela_final_completa)

library(ggplot2)

ggplot(tabela_final_completa, aes(x = Matriz, y = SR_media, colour = Amostra)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.6) +
        labs(y = "Slope Ratio (SR)", x = "Matriz") +
        theme_minimal()

ggplot(tabela_final_completa, aes(x = Matriz, y = SUVA_media, colour = Amostra)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.6) +
        labs(y = "SUVA", x = "Matriz") +
        theme_minimal()

ggplot(tabela_final_completa, aes(x = Matriz, y = E2_E3, colour = Amostra)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.6) +
        labs(y = "E2/E3", x = "Matriz") +
        theme_minimal()

ggplot(tabela_final_completa, aes(x = Matriz, y = aromaticidade_media, colour = Amostra)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.6) +
        labs(y = "Aromaticidade (fenólicos/carboxílicos)", x = "Matriz") +
        theme_minimal()


ggplot(tabela_final_completa, aes(x = Matriz, y = carboxilicos_media, colour = Amostra)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.6) +
        labs(y = " Teor de carboxílicos", x = "Matriz") +
        theme_minimal()


ggplot(tabela_final_completa, aes(x = Matriz, y = fenolicos_media, colour = Amostra)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.6) +
        labs(y = "Teor de fenólicos", x = "Matriz") +
        theme_minimal()


ggplot(tabela_final_completa, aes(x = Matriz, y = total_media, colour = Amostra)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.6) +
        labs(y = "Acidez total", x = "Matriz") +
        theme_minimal()

ggplot(tabela_final_completa,
       aes(x = SR_media, y = E2_E3, color = Matriz)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        theme_minimal() +
        labs(x = "Slope Ratio (SR)", y = "E2/E3")

ggplot(tabela_final_completa,
       aes(x = SR_media, y = aromaticidade_media, color = Matriz)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        theme_minimal() +
        labs(x = "Slope Ratio (SR)", y = "Aromaticidade")

ggplot(tabela_final_completa,
       aes(x = SUVA_media, y = carboxilicos_media, color = Matriz)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        theme_minimal()

ggplot(tabela_final_completa,
       aes(x = SUVA_media, y = fenolicos_media, color = Matriz)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        theme_minimal()

dados_pca <- tabela_final_completa %>%
        select(
                SUVA_media, E2_E3, E4_E5, SR_media,
                carboxilicos_media, fenolicos_media, aromaticidade_media
        ) %>%
        scale()

pca <- prcomp(dados_pca, center = TRUE, scale. = TRUE)

biplot(pca, scale = 0)

ggplot(tabela_final_completa,
       aes(x = Coleta, y = SR_media, color = Matriz, group = Matriz)) +
        geom_point(size = 3) +
        geom_line() +
        theme_minimal()

library(dplyr)
library(tidyr)

indices <- tabela_final_completa %>%
        mutate(Estacao = case_when(
                Coleta %in% c(1, 4, 6) ~ "Chuvosa",
                Coleta %in% c(2, 3, 5) ~ "Seca"
        ))

rodar_estatistica <- function(df, variavel) {
        
        norm_seca <- shapiro.test(df[[variavel]][df$Estacao == "Seca"])$p.value
        norm_chuv <- shapiro.test(df[[variavel]][df$Estacao == "Chuvosa"])$p.value
        
        if(norm_seca > 0.05 & norm_chuv > 0.05) {
                teste <- t.test(df[[variavel]] ~ df$Estacao)
                metodo <- "Teste t (Paramétrico)"
        } else {
                teste <- wilcox.test(df[[variavel]] ~ df$Estacao)
                metodo <- "Wilcoxon (Não-paramétrico)"
        }
        
        return(list(Metodo = metodo, P_Valor = teste$p.value))
}

library(purrr)

vars_interesse<- c("SUVA_media","E2_E3", "E4_E5","E4_E6","SR_media","carboxilicos_media", "fenolicos_media", "total_media", "aromaticidade_media")

tabela_estatistica <- indices %>%
        group_by(Matriz) %>%
        group_modify(~ {
                map_df(vars_interesse, function(v) {
                        res <- rodar_estatistica(.x, v)
                        tibble(Variavel = v, Metodo = res$Metodo, P_Valor = res$P_Valor)
                })
        })

tabela_estatistica %>% 
        arrange(P_Valor) %>% 
        print(n = Inf)

library(ggplot2)
library(ggpubr)

variavel_plot <- c("E4_E5")

ggplot(indices, aes(x = Estacao, y = .data[[variavel_plot]], fill = Estacao)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, alpha = 0.5) +
        facet_wrap(~Matriz, scales = "free_y") +
        stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
        scale_fill_manual(values = c("Chuvosa" = "#3B9AB2", "Seca" = "#EBCC2A")) +
        labs(
                title = paste("Variação Sazonal de", variavel_plot),
                subtitle = "Comparação não-paramétrica (Wilcoxon)",
                x = "Estação",
                y = variavel_plot
        ) +
        theme_minimal() +
        theme(legend.position = "none", strip.text = element_text(face = "bold", size = 12))
###############################################

indices <- indices %>%
        mutate(Tipo = if_else(str_detect(Amostra, "^P"), "Permanente", "Sazonal"))

tabela_estatistica_cruzada <- indices %>%
        group_by(Matriz, Tipo) %>%
        group_modify(~ {
                map_df(vars_interesse, function(v) {
                        res <- rodar_estatistica(.x, v)
                        tibble(Variavel = v, P_Valor = res$P_Valor)
                })
        }) %>%
        arrange(P_Valor)

print(tabela_estatistica_cruzada, n = Inf)

ggplot(indices, aes(x = Estacao, y = .data[[variavel_plot]], fill = Estacao)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, alpha = 0.5) +
        # Facet_grid cria uma matriz: Matriz nas linhas, Tipo nas colunas
        facet_grid(Matriz ~ Tipo, scales = "free_y") + 
        stat_compare_means(method = "wilcox.test", label = "p.format") +
        scale_fill_manual(values = c("Chuvosa" = "#3B9AB2", "Seca" = "#EBCC2A")) +
        labs(
                title = paste("Variação de", variavel_plot, "por Tipo de Wetland e Estação"),
                x = "Estação",
                y = variavel_plot
        ) +
        theme_bw()

ggplot(indices, aes(x = Estacao, y = SR_media, fill = Estacao)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.15, alpha = 0.5) +
        facet_grid(Matriz ~ Tipo, scales = "free_y") + 
        stat_compare_means(method = "wilcox.test", label = "p.format") +
        scale_fill_manual(values = c("Chuvosa" = "#3B9AB2", "Seca" = "#EBCC2A")) +
        labs(
                title = paste("Variação de", variavel_plot, "por Tipo de Wetland e Estação"),
                x = "Estação",
                y = variavel_plot
        ) +
        theme_bw()
