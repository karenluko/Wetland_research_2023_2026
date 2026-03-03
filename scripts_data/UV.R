library(tidyverse)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data")

pasta <- "~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data/UV-vis"

arquivos <- list.files(pasta, pattern = "\\.csv$", full.names = TRUE)

ler_espectro <- function(arq) {
        
        linhas <- readLines(arq)
        inicio <- which(grepl("^WL", linhas))
        
        dados <- read.csv(arq, skip = inicio - 1, check.names = TRUE)
        
        dados <- dados[, 1:2]          # WL e Abs
        colnames(dados) <- c("WL", "Abs")
        
        nome_amostra <- tools::file_path_sans_ext(basename(arq))
        colnames(dados)[2] <- nome_amostra
        
        return(dados)
}


tabela_final <- arquivos %>%
        purrr::map(ler_espectro) %>%
        purrr::reduce(full_join, by = "WL")

write.csv(tabela_final, "Espectros_AH_AF_unificados.csv", row.names = FALSE)

#################################################################################

dados_long <- tabela_final %>%
        pivot_longer(
                cols = -WL,
                names_to = "nome_bruto",
                values_to = "Abs"
        )

extrair_info <- function(nome) {
        
        partes <- str_split(nome, " ", simplify = TRUE)
        
        amostra  <- partes[1]
        coleta   <- partes[2]
        matriz   <- partes[3]
        
        ultimo <- partes[4]
        
        diluicao <- str_extract(ultimo, "\\d+X")
        replicata <- str_extract(ultimo, "[a-zA-Z]$")
        
        tibble(
                Amostra = amostra,
                Coleta = coleta,
                Matriz = matriz,
                Diluicao = diluicao,
                Replicata = replicata
        )
}

metadados <- map_dfr(dados_long$nome_bruto, extrair_info)

dados_long <- bind_cols(dados_long, metadados)

dados_long <- dados_long %>%
        mutate(
                Matriz = if_else(Matriz %in% c("AH", "AF"), Matriz, NA_character_)
        )

dados_long <- dados_long %>%
        mutate(
                Replicata = if_else(Replicata == "X", "a", Replicata)
        )

dados_long <- dados_long %>%
        mutate(
                Diluicao = if_else(is.na(Diluicao), "1X", Diluicao)
        )



glimpse(dados_long)
        
distinct(dados_long, Matriz)
distinct(dados_long, Diluicao)
count(dados_long, Replicata)


###########################################################################
#Calculo das razões

dados_long <- dados_long %>%
        mutate(
                fator_diluicao = as.numeric(str_extract(Diluicao, "\\d+"))
        )

abs_254 <- dados_long %>%
        filter(WL == 254) %>%
        filter(Abs >= 0.1, Abs <= 1.05)

melhor_diluicao <- abs_254 %>%
        group_by(Amostra, Coleta, Matriz) %>%
        slice_min(fator_diluicao, n = 1) %>%
        ungroup()

dados_best <- dados_long %>%
        semi_join(
                melhor_diluicao,
                by = c("Amostra", "Coleta", "Matriz", "Diluicao")
        )

#SUVA

suva <- dados_best %>%
        filter(WL == 254) %>%
        mutate(
                SUVA = (Abs * fator_diluicao*100) /30
        ) %>%
        select(Amostra, Coleta, Matriz, Diluicao, Replicata, SUVA)

suva_med <- suva %>%
        group_by(Amostra, Coleta, Matriz) %>%
        summarise(
                SUVA_media = mean(SUVA, na.rm = TRUE),
                SUVA_sd = sd(SUVA, na.rm = TRUE),
                .groups = "drop"
        )

#Razões espectrais

dados_ratio <- dados_best %>%
        group_by(Amostra, Coleta, Matriz, Diluicao, WL) %>%
        summarise(
                Abs = mean(Abs, na.rm = TRUE),
                .groups = "drop"
        )

calc_razao <- function(df, wl1, wl2, nome) {
        
        df %>%
                filter(WL %in% c(wl1, wl2)) %>%
                select(Amostra, Coleta, Matriz, Diluicao, WL, Abs) %>%
                pivot_wider(
                        names_from = WL,
                        values_from = Abs
                ) %>%
                mutate(
                        !!nome := .[[as.character(wl1)]] / .[[as.character(wl2)]]
                )
}

E2E3 <- calc_razao(dados_ratio, 250, 365, "E2_E3")
E4E5 <- calc_razao(dados_ratio, 465, 520, "E4_E5")
E4E6 <- calc_razao(dados_ratio, 465, 665, "E4_E6")

summary(E2E3$E2_E3)
summary(E4E5$E4_E5)
summary(E4E6$E4_E6)

indices <- suva_med %>%
        left_join(E2E3, by = c("Amostra", "Coleta", "Matriz")) %>%
        left_join(E4E5, by = c("Amostra", "Coleta", "Matriz")) %>%
        left_join(E4E6, by = c("Amostra", "Coleta", "Matriz"))

#Slopes S275–295 e S350–400 e Slope ratio (SR) = S275–295 / S350–400

dados_nest <- dados_best %>%
        group_by(Amostra, Coleta, Matriz, Diluicao, Replicata) %>%
        nest()

calc_slope <- function(df, wl_min, wl_max) {
        
        df_int <- df %>%
                filter(WL >= wl_min, WL <= wl_max)
        
        if (nrow(df_int) < 5) return(NA_real_)
        
        coef(lm(Abs ~ WL, data = df_int))[2]
}

slopes <- dados_nest %>%
        mutate(
                S275_295 = map_dbl(data, calc_slope, 275, 295),
                S350_400 = map_dbl(data, calc_slope, 350, 400),
                SR = S275_295 / S350_400
        ) %>%
        select(-data)

summary(slopes$S275_295)
summary(slopes$S350_400)
summary(slopes$SR)


plot_slope_diag <- function(df, wl_min, wl_max, titulo = NULL) {
        
        df_int <- df %>%
                filter(WL >= wl_min, WL <= wl_max)
        
        mod <- lm(Abs ~ WL, data = df_int)
        r2 <- summary(mod)$r.squared
        
        ggplot(df_int, aes(x = WL, y = Abs)) +
                geom_point(size = 2) +
                geom_smooth(method = "lm", se = FALSE) +
                labs(
                        title = titulo,
                        subtitle = paste0(
                                wl_min, "-", wl_max, " nm | ",
                                "Slope = ", signif(coef(mod)[2], 3),
                                " | R² = ", round(r2, 3)
                        ),
                        x = "WL (nm)",
                        y = "Absorbância"
                ) +
                theme_minimal()
}

exemplo <- dados_best %>%
        filter(
                Amostra == "PB",
                Coleta == "1",
                Matriz == "AF",
                Diluicao == "6X",
                Replicata == "a"
        )

plot_slope_diag(
        exemplo,
        wl_min = 275,
        wl_max = 295,
        titulo = "PB 1 AF 6Xa – S275–295"
)

plot_slope_diag(
        exemplo,
        wl_min = 350,
        wl_max = 400,
        titulo = "PB 1 AF 6Xa – S350–400"
)

calc_r2 <- function(df, wl_min, wl_max) {
        
        df_int <- df %>%
                filter(WL >= wl_min, WL <= wl_max)
        
        if (nrow(df_int) < 5) return(NA_real_)
        
        summary(lm(Abs ~ WL, data = df_int))$r.squared
}

slopes_qc <- dados_nest %>%
        mutate(
                R2_275_295 = map_dbl(data, calc_r2, 275, 295),
                R2_350_400 = map_dbl(data, calc_r2, 350, 400)
        )

summary(slopes_qc$R2_275_295)
summary(slopes_qc$R2_350_400)

#> 0.95

casos_nao_lineares <- slopes_qc %>%
        filter(R2_350_400 < 0.95) %>%
        arrange(R2_350_400)

casos_nao_lineares %>%
        count(Matriz)
casos_nao_lineares %>%
        count(Diluicao)

caso_exemplo <- casos_nao_lineares %>%
        slice(1)

dados_exemplo <- dados_best %>%
        semi_join(
                caso_exemplo,
                by = c("Amostra", "Coleta", "Matriz", "Diluicao", "Replicata")
        )

plot_slope_diag(
        dados_exemplo,
        wl_min = 350,
        wl_max = 400,
        titulo = paste(
                caso_exemplo$Amostra,
                caso_exemplo$Coleta,
                caso_exemplo$Matriz,
                caso_exemplo$Diluicao,
                caso_exemplo$Replicata
        )
)

slopes_qc %>%
        mutate(flag = if_else(R2_350_400 < 0.95, "baixo_R2", "ok")) %>%
        count(flag)


########################################################

slopes_med <- slopes %>%
        group_by(Amostra, Coleta, Matriz, Diluicao) %>%
        summarise(
                S275_295_media = mean(S275_295, na.rm = TRUE),
                S275_295_sd    = sd(S275_295, na.rm = TRUE),
                S350_400_media = mean(S350_400, na.rm = TRUE),
                S350_400_sd    = sd(S350_400, na.rm = TRUE),
                SR_media       = mean(SR, na.rm = TRUE),
                SR_sd          = sd(SR, na.rm = TRUE),
                .groups = "drop"
        )

head(slopes_med)

indices_opticos <- indices %>%
        select(
                Amostra, Coleta, Matriz, Diluicao,
                SUVA_media, SUVA_sd,
                E2_E3, E4_E5, E4_E6
        )

tabela_final_indices <- indices_opticos %>%
        left_join(
                slopes_med,
                by = c("Amostra", "Coleta", "Matriz", "Diluicao")
        )

glimpse(tabela_final_indices)

tabela_final_indices <- tabela_final_indices %>%
        arrange(Amostra, Coleta, Matriz)

write.csv(
        tabela_final_indices,
        "Tabela_indices_opticos_e_slopes.csv",
        row.names = FALSE
)


library(writexl)

write_xlsx(
        tabela_final_indices,
        "Tabela_indices_opticos_e_slopes.xlsx"
)

