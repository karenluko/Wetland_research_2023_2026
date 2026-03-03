
# ANALISE REDOX das amostras de AH e AF

library(tidyverse)
library(dplyr)
library(ggplot2)
library(purrr)
library(ggpubr)

# 1) Leitura dos arquivos

pasta <- "~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data/reducao"
arquivos <- list.files(pasta, pattern = "\\.csv$", full.names = TRUE)

extrair_meta <- function(nome) {
        
        meta <- tibble(
                arquivo = nome,
                
                amostra = str_extract(nome, "PB|PM|PC|SB|SM|SC"),
                
                coleta = str_extract(nome, "(?<= )[1-6](?= )"),
                
                matriz = str_extract(nome, "AH|AF"),
                
                tratamento = str_extract(nome, "NTL|RDZ"),
                
                diluicao = str_extract(nome, "\\d+X|X\\d+"),
                
                replicata = str_extract(nome, "[ab](?=\\.csv$)")
        )
        
        meta <- meta %>%
                mutate(
                        diluicao = ifelse(is.na(diluicao), "ND", diluicao),
                        diluicao = str_replace_all(diluicao, "X", ""),
                        
                        tipo = case_when(
                                amostra %in% c("PB","PM","PC") ~ "Permanente",
                                amostra %in% c("SB","SM","SC") ~ "Sazonal",
                        )
                )
        
        return(meta)
}


ler_espectro <- function(file) {
        nome <- basename(file)
        dados <- read.csv2(file, skip = 10)
        if(ncol(dados) == 1) dados <- read.csv(file, skip = 10)
        
        dados <- dados[,1:2]
        colnames(dados) <- c("WL","Abs")
        
        meta <- extrair_meta(nome)
        
        dados %>%
                mutate(across(everything(), as.numeric)) %>%
                dplyr::filter(!is.na(WL), WL >= 300, WL <= 450) %>%
                bind_cols(meta %>% select(-arquivo))
}


dados_completos <- map_dfr(arquivos, ler_espectro) %>%
        mutate(across(where(is.character), str_trim)) %>%
        group_by(tipo, matriz, tratamento, amostra, coleta)


dados_filtro <- dados_completos %>%
        dplyr::filter(WL >= 335, WL <= 440) %>%
        arrange(tipo, matriz, tratamento, amostra, coleta, WL)

dados_limpos <- dados_filtro %>%
        group_by(amostra, coleta, tipo, matriz, tratamento, diluicao, WL) %>%
        summarise(
                Abs = mean(Abs, na.rm = TRUE), 
                .groups = "drop"
        )

# selecionar mesma diluição para NTL e RDZ
amostras_comparaveis <- dados_limpos %>%
        group_by(amostra, coleta, matriz) %>%
        dplyr::filter(n_distinct(diluicao) == 1) %>% 
        select(amostra, coleta, matriz, tipo) %>%
        distinct()

analise_final <- dados_limpos %>%
        inner_join(amostras_comparaveis, by = c("amostra", "coleta", "matriz", "tipo"))

#library(writexl)

#write_xlsx(
        #analise_final,
        #"analise_final.xlsx"
#)


cv_por_amostra <- analise_final %>%
        select(amostra, coleta, matriz, tipo, WL, tratamento, Abs) %>%
        pivot_wider(
                names_from = tratamento,
                values_from = Abs
        ) %>%
        dplyr::filter(!is.na(NTL), !is.na(RDZ)) %>% 
        mutate(
                media_par = (NTL + RDZ) / 2,
                desvio_par = abs(NTL - RDZ) / sqrt(2),
                cv_amostra = (desvio_par / media_par) * 100
        )

cv_resumo <- cv_por_amostra %>%
        group_by(matriz, WL) %>%
        summarise(cv_mean = mean(cv_amostra, na.rm = TRUE),
                  .groups = "drop")

ggplot(cv_resumo, aes(x = WL, y = cv_mean)) +
        geom_line(size = 0.9) +
        facet_wrap(~matriz, scales = "free_y") +
        labs(
                title = "Comparação entre NTL e RDZ dentro de cada par",
                x = "Wavelength (nm)",
                y = "CV (%)"
        ) +
        theme_bw()

# gráfico
ggplot(cv_por_amostra, aes(x = WL, y = cv_amostra, colour=amostra)) +
        geom_line(size = 0.7, alpha = 0.6) + 
        facet_wrap(~matriz, scales = "free_y") + 
        labs(
                title = "Comparação entre NTL e RDZ dentro de cada par",
                x = "Wavelength (nm)",
                y = "CV (%)",
        ) +
        theme_bw() +
        theme(legend.position = "bottom")

razoes_redox <- dados_completos %>%
        select(tipo, matriz, amostra, coleta, tratamento, WL, Abs) %>%
        group_by(tipo, matriz, amostra, coleta, tratamento) %>%
        reframe(
                soma_Ox = case_when(
                        tipo == "Sazonal"  & matriz == "AH"~ sum(Abs[WL >= 325 & WL <= 375], na.rm = TRUE),
                        tipo == "Sazonal"  & matriz == "AF"~ sum(Abs[WL >= 325 & WL <= 375], na.rm = TRUE),
                        tipo == "Permanente" & matriz == "AH" ~ sum(Abs[WL >= 325 & WL <= 375], na.rm = TRUE),
                        tipo == "Permanente" & matriz == "AF" ~ sum(Abs[WL >= 325 & WL <= 375], na.rm = TRUE)
                ),
                soma_An = case_when(
                        tipo == "Sazonal"  & matriz == "AH" ~ sum(Abs[WL >= 427 & WL <= 440], na.rm = TRUE),
                        tipo == "Sazonal"  & matriz == "AF" ~ sum(Abs[WL >= 427 & WL <= 440], na.rm = TRUE),
                        tipo == "Permanente"  & matriz == "AH" ~ sum(Abs[WL >= 427 & WL <= 440], na.rm = TRUE),
                        tipo == "Permanente"  & matriz == "AF" ~ sum(Abs[WL >= 427 & WL <= 440], na.rm = TRUE)
                )
        ) %>%
        distinct() 
        

ggplot(razoes_redox[razoes_redox$tipo == "Sazonal", ], aes(x = soma_An, y = soma_Ox, color = tratamento)) +
        #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(aes(shape = amostra), size = 3.5, alpha = 0.7) +
        #geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, linetype = "dashed", alpha = 0.3) +
        #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
        #label.x.npc = "left", label.y.npc = "top") +
        labs(  x = "Σ325-375",
               y = "Σ437-440",
        ) +
        theme_bw() +
        facet_wrap(tipo~matriz+coleta, scales = "free") 

ggplot(razoes_redox[razoes_redox$tipo == "Permanente", ], aes(x = soma_An, y = soma_Ox, color = tratamento)) +
        #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(aes(shape = amostra), size = 3.5, alpha = 0.7) +
        #geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, linetype = "dashed", alpha = 0.3) +
        #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                 #label.x.npc = "left", label.y.npc = "top") +
        labs(  x = "Σ325-375",
               y = "Σ437-440",
        ) +
        theme_bw() +
        facet_wrap(tipo~matriz+coleta, scales = "free")

dados_norm <- dados_completos %>%
        group_by(amostra, coleta, matriz, tratamento,
                 diluicao, replicata, tipo) %>%
        summarise(
                abs_325_375 = sum(Abs[WL >= 325 & WL <= 375], na.rm = TRUE),
                abs_437_440 = sum(Abs[WL >= 437 & WL <= 440], na.rm = TRUE),
                abs_norm = abs_325_375 / abs_437_440,
                .groups = "drop"
        )

dados_delta <- dados_norm %>%
        select(amostra, coleta, matriz, tipo, replicata,
               tratamento, abs_norm) %>%
        tidyr::pivot_wider(
                names_from = tratamento,
                values_from = abs_norm
        ) %>%
        mutate(delta_abs = NTL-RDZ)

delta_resumo <- dados_delta %>%
        group_by(amostra, coleta, tipo, matriz) %>%
        summarise(
                delta_media = mean(delta_abs, na.rm = TRUE),
                delta_sd = sd(delta_abs, na.rm = TRUE),
                n = n(),
                .groups = "drop"
        )

write_xlsx(delta_resumo, "delta_resumo.xlsx")

ggplot(delta_resumo,
        aes(coleta, delta_media,
            color = amostra,
            group = amostra)) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        theme_bw() +
        facet_wrap(~matriz)

ggplot(delta_resumo,
       aes(x = coleta,
           y = delta_media,
           color = matriz,
           group = interaction(amostra, matriz))) +
        
        geom_hline(yintercept = 0, linetype = 2) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        
        facet_wrap(~tipo) +
        
        labs(
                x = "Coleta",
                y = expression(Delta*"ABS (NTL-RDZ)"),
                color = "Fração"
        ) +
        theme_bw()

ggplot(delta_resumo,
       aes(x = matriz, y = delta_media, fill = matriz)) +
        geom_boxplot(alpha = 0.7) +
        facet_wrap(~tipo) +
        theme_bw() +
        labs(
                x = "Fração",
                y = expression(Delta*"ABS")
        )

ggplot(delta_resumo,
       aes(x = coleta,
           y = amostra,
           fill = delta_media)) +
        geom_tile() +
        facet_wrap(~matriz) +
        scale_fill_gradient2(
                low = "blue",
                mid = "white",
                high = "red",
                midpoint = 0
        ) +
        theme_bw()

delta_plot <- delta_resumo %>%
        dplyr::mutate(coleta = as.integer(coleta)) %>%
        left_join(
                funcionais_med %>%
                        select(amostra, coleta, matriz, fenolicos_media, carboxilicos_media, total_media, aromaticidade_media),
                by = c("amostra", "coleta", "matriz")
        )

library(ggpubr)
library(dplyr)
library(broom)


ggplot(delta_plot,
       aes(x = fenolicos_media,
           y = delta_media,
           color = tipo)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        stat_regline_equation(
                aes(label = paste(..rr.label.., sep = "  "),
                    color = matriz),
                label.x.npc = "left",
                label.y.npc = "top"
        ) +
        facet_wrap(tipo~matriz) +
        theme_bw() +
        labs(
                x = expression("Teor de fenólicos (mmol L"^{-1}*")"),
                y = expression(Delta*"ABS"))

ggplot(delta_plot,
       aes(x = coleta,
           y = delta_media,
           group = amostra,
           color = aromaticidade_media)) +
        
        geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
        
        geom_line(size = 1) +
        geom_point(size = 3) +
        
        facet_wrap(~matriz) +
        
        scale_color_viridis_c(name = "Aromaticidade") +
        
        theme_bw() +
        
        labs(
                x = "Coleta",
                y = expression(Delta*"ABS (RDZ - NTL)")
        )

library(dplyr)

dados_pca <- delta_plot %>%
        select(
                amostra,
                coleta,
                tipo,
                matriz,
                delta_media,
                fenolicos_media,
                carboxilicos_media,
                total_media,
                aromaticidade_media
        ) %>%
        na.omit()

pca <- prcomp(
        dados_pca %>%
                select(delta_media,
                       fenolicos_media,
                       carboxilicos_media,
                       total_media,
                       aromaticidade_media),
        scale. = TRUE
)
summary(pca)

library(factoextra)

fviz_pca_biplot(
        pca,
        geom.ind = "point",
        habillage = dados_pca$matriz,
        addEllipses = TRUE,
        label = "var",
        repel = TRUE
)

fviz_pca_biplot(
        pca,
        geom.ind = "point",
        habillage = dados_pca$tipo,
        addEllipses = TRUE,
        label = "var",
        repel = TRUE
)

scores <- as.data.frame(pca$x)

scores <- cbind(
        dados_pca %>%
                select(amostra, coleta, tipo, matriz),
        scores
)

ggplot(scores,
       aes(x = PC1, y = PC2,
           color = matriz)) +
        
        geom_point(size = 3) +
        geom_text(aes(label = coleta), vjust = -1)+
        geom_path(
                aes(group = interaction(amostra, matriz)),
                arrow = arrow(type = "closed", length = unit(0.15,"cm"))
        )+
        facet_wrap(~tipo) +
        
        theme_bw() +
        
        labs(
                x = "PC1",
                y = "PC2",
                color = "Fração"
        )

library(dplyr)
library(rstatix)

dados_delta <- dados_delta %>%
        mutate(estacao = case_when(
                coleta %in% c(1, 4, 6) ~ "chuvosa",
                coleta %in% c(2, 3, 5) ~ "seca",
                TRUE ~ NA_character_
        ))

kruskal_eac <- dados_delta %>%
        group_by(matriz, estacao) %>%
        kruskal_test(delta_abs ~ tipo)

kruskal_eac


fenolicos_medio <- dados_quimicos_long %>%
        dplyr::filter(variavel == "fenolicos") %>%
        group_by(matriz, coleta, tipo, amostra, estacao) %>%
        summarise(valor_medio = mean(valor, na.rm = TRUE), .groups = "drop")

delta_medio <- dados_delta %>%
        group_by(matriz, coleta, tipo, amostra, estacao) %>%
        summarise(delta_medio = mean(delta_abs, na.rm = TRUE), .groups = "drop")

dados_correlacao <- merge(
        fenolicos_medio, 
        delta_medio, 
        by = c("matriz", "coleta", "tipo", "amostra", "estacao")
)

library(ggpubr)
library(showtext)

font_add("Helvetica", regular = "Helvetica.ttc")
showtext_auto()


ggplot(dados_correlacao, aes(x = valor_medio, y = delta_medio)) +
        geom_point(aes(color = as.factor(coleta)), size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black") +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                 label.x.npc = "left", label.y.npc = "top", family = "Helvetica") +
        facet_grid(matriz ~ tipo, scales = "free") +
        scale_color_brewer(palette = "RdYlBu", name = "Coleta") +
        theme_bw(base_family = "Helvetica") +
        labs(
                x = "Teor de sítios fenólicos (mmol g⁻¹ C)",
                y = "EAC (ΔAbs)"
        )

ggplot(tabela_final_completa, aes(x = fenolicos_media, y = E4_E5)) +
        geom_point(aes(color = as.factor(Coleta)), size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black") +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                 label.x.npc = "left", label.y.npc = "top", family = "Helvetica") +
        facet_grid(Matriz ~ ., scales = "free") +
        scale_color_brewer(palette = "RdYlBu", name = "Coleta") +
        theme_bw(base_family = "Helvetica") +
        labs(
                x = "Teor de sítios fenólicos (mmol g⁻¹ C)",
                y = "E4/E6"
        )

ggplot(tabela_final_completa, aes(x = fenolicos_media, y = SUVA_media)) +
        geom_point(aes(color = as.factor(Coleta)), size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black") +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                 label.x.npc = "left", label.y.npc = "top", family = "Helvetica") +
        facet_grid(Matriz ~ ., scales = "free") +
        scale_color_brewer(palette = "RdYlBu", name = "Coleta") +
        theme_bw(base_family = "Helvetica") +
        labs(
                x = "Teor de sítios fenólicos (mmol g⁻¹ C)",
                y = "SUVA"
        )

ggplot(tabela_final_completa, aes(x = carboxilicos_media, y = E4_E5)) +
        geom_point(aes(color = as.factor(Coleta)), size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black") +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                 label.x.npc = "left", label.y.npc = "top", family = "Helvetica") +
        facet_grid(Matriz ~ ., scales = "free") +
        scale_color_brewer(palette = "RdYlBu", name = "Coleta") +
        theme_bw(base_family = "Helvetica") +
        labs(
                x = "Carboxílicos (mmol g⁻¹ C)",
                y = "E4/E6"
        )

ggplot(tabela_final_completa, aes(x = carboxilicos_media, y = SUVA_media )) +
        geom_point(aes(color = as.factor(Coleta)), size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black") +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                 label.x.npc = "left", label.y.npc = "top", family = "Helvetica") +
        facet_grid(Matriz ~ ., scales = "free") +
        scale_color_brewer(palette = "RdYlBu", name = "Coleta") +
        theme_bw(base_family = "Helvetica") +
        labs(
                x = "Carboxílicos (mmol g⁻¹ C)",
                y = "SUVA"
        )

ggplot(tabela_final_completa, aes(x = E4_E6, y = E2_E3 )) +
        geom_point(aes(color = as.factor(Coleta)), size = 3, alpha = 0.8) +
        geom_smooth(method = "lm", color = "black") +
        stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
                 label.x.npc = "left", label.y.npc = "top", family = "Helvetica") +
        facet_grid(Matriz ~ ., scales = "free") +
        scale_color_brewer(palette = "RdYlBu", name = "Coleta") +
        theme_bw(base_family = "Helvetica") +
        labs(
                x = "E4/E5",
                y = "E2_E3"
        )

