setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data")
dados_raw <- read.csv2("EGA_dados.csv", header = FALSE, stringsAsFactors = FALSE)

library(ggplot2)
library(dplyr)
library(tidyverse)
library(janitor)
library(purrr)
library(ggh4x)
library(ggforce)
library(pracma)
library(writexl)
library(ggpubr)

# metadados
metadados <- tibble(
        Amostra = c(rep("P01", 4), rep("P03", 4), rep("P02", 4),  
                    rep("P01", 4), rep("P03", 4), rep("P02", 4)), 
        Matriz = rep(c("HA", "Soil", "HA", "Soil", "HA", "Soil", "HA", "Soil", "HA", "Soil", "HA", "Soil"), 2),
        Tipo_Area = c(rep("Permanent", 12), rep("Seasonal", 12)),
        Estacao = rep(c("Rainy", "Rainy", "Dry", "Dry", "Rainy", "Rainy", "Dry", "Dry", "Rainy", "Rainy", "Dry", "Dry"), 2)
)

dados_raw<-dados_raw[-c(1:4),]
dados_raw <- dados_raw %>%
        row_to_names(row_number = 1)
colnames(dados_raw)[-1] <- paste0("Counts_", seq_along(colnames(dados_raw)[-1]))

dados_long <- dados_raw %>%
        pivot_longer(
                cols = -Time,
                names_to = "Amostra", 
                values_to = "Counts"
        )

dados_long <- dados_long %>%
        mutate(
                Tipo_Area = rep(metadados$Tipo_Area, length.out = nrow(dados_long)),
                Estacao = rep(metadados$Estacao, length.out = nrow(dados_long)),
                Matriz = rep(metadados$Matriz, length.out = nrow(dados_long)),
                Amostra = rep(metadados$Amostra, length.out = nrow(dados_long))
        )

head(dados_long)
table(dados_long$Matriz, dados_long$Amostra)

dados <- dados_long %>% 
        mutate(Time = as.numeric(Time),
               Temperatura = 50 + (20 * Time),
               Counts = as.numeric(Counts))

table(dados$Matriz, dados$Amostra)

dados <- dados %>% 
mutate(Estacao = factor(Estacao, levels = c("Rainy","Dry")),
       Tipo_Area = factor(Tipo_Area, levels = c("Permanent","Seasonal")),
       Matriz = factor(Matriz, levels = c("Soil","HA")),
       Amostra = factor(Amostra, levels = c("P01","P02", "P03"))) 


# REMOÇÃO DO PICO INTERFERENTE
dados_filtrados <- dados %>% dplyr::mutate()

ultimo_valor <- dados_filtrados %>%
        filter(Temperatura < 600, !is.na(Counts)) %>%
        summarize(ultimo_counts = last(Counts, order_by = Temperatura)) %>%
        pull(ultimo_counts)

dados_filtrados <- dados_filtrados %>%
        arrange(Temperatura) %>%
        mutate(Counts = case_when(
                (Tipo_Area == "Permanent" & Matriz == "HA" & Estacao == "Dry" & Temperatura >= 590 & Temperatura <= 750) ~ 
                        lag(Counts, default = ultimo_valor) * (1 - (Temperatura - 590) / (750 - 590)),
                
                (Tipo_Area == "Seasonal" & Matriz == "HA" & Estacao == "Dry" & Temperatura >= 485 & Temperatura <= 750) ~ 
                        lag(Counts, default = ultimo_valor) * (1 - (Temperatura - 485) / (750 - 485)),
                
                TRUE ~ Counts
        ))


y_limite <- range(dados_filtrados$Counts, na.rm = TRUE)

# Gráfico
dados_filtrados <- dados_filtrados %>%
        mutate(Counts_Ajustado = ifelse(Matriz == "Soil", Counts * 4, Counts))#aumentar escala para solos em 4x

linhas_verticais <- read.csv2("Linhas_verticais_temperatura.csv", header = TRUE, stringsAsFactors = FALSE)
linhas_verticais <- linhas_verticais %>%
        mutate(Estacao = factor(Estacao, levels = c("Rainy", "Dry")))%>%
        mutate(Matriz = factor(Matriz, levels = c("Soil", "HA")))

p <- ggplot(dados_filtrados, aes(x = Temperatura, y = Counts_Ajustado)) +
        stat_smooth(method = "loess", formula = y ~ x, span = 0.1, color = "black", size = 0.5, se = FALSE) +
        geom_vline(data = linhas_verticais, aes(xintercept = Temperatura), linetype = "dashed", color = "red", linewidth = 0.3) +
        geom_text(data = linhas_verticais, aes(x = Temperatura, 
                                               y = max(dados_filtrados$Counts_Ajustado, na.rm = TRUE) * 0.85, 
                                               label = paste0(Temperatura, "°C")), 
                  angle = 90, vjust = -0.5, size = 2.0, color = "gray40")+
        facet_nested(fct_relevel(Tipo_Area, "Permanent", "Seasonal") + fct_relevel(Estacao, "Rainy", "Dry") ~ Amostra + Matriz, space = "free_x") +
        labs(x = "Temperature (°C)", y = "") +
        scale_x_continuous(limits = c(50, 750), breaks = seq(50, 750, by = 50)) +
        scale_y_continuous(limits = y_limite) +
        theme_light() +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
              legend.position = "none")

print(p)

ggsave("grafico_picos.tiff", plot = p, width = 10, height = 8, dpi = 300, units = "in", device = "tiff")
                                
                                ########################################

#Cálculos dos pools de carbono

linhas_verticais <- read.csv2("Linhas_verticais_temperatura.csv", stringsAsFactors = FALSE)
pools_ega <- read.csv2("Pools_EGA.csv", stringsAsFactors = FALSE)
dados_raw <- read.csv2("EGA_dados.csv", header = FALSE, stringsAsFactors = FALSE)

metadados <- tibble(
        Amostra = c(rep("P01", 4), rep("P03", 4), rep("P02", 4),  
                    rep("P01", 4), rep("P03", 4), rep("P02", 4)), 
        Matriz = rep(c("HA", "Soil", "HA", "Soil", "HA", "Soil", "HA", "Soil", "HA", "Soil", "HA", "Soil"), 2),
        Tipo_Area = c(rep("Permanent", 12), rep("Seasonal", 12)),
        Estacao = rep(c("Rainy", "Rainy", "Dry", "Dry", "Rainy", "Rainy", "Dry", "Dry", "Rainy", "Rainy", "Dry", "Dry"), 2)
)

dados_raw<-dados_raw[-c(1:4),]
dados_raw <- dados_raw %>%
        row_to_names(row_number = 1)
colnames(dados_raw)[-1] <- paste0("Counts_", seq_along(colnames(dados_raw)[-1]))

dados_long <- dados_raw %>%
        pivot_longer(
                cols = -Time,
                names_to = "Amostra", 
                values_to = "Counts"
        )

dados_long <- dados_long %>%
        mutate(
                Tipo_Area = rep(metadados$Tipo_Area, length.out = nrow(dados_long)),
                Estacao = rep(metadados$Estacao, length.out = nrow(dados_long)),
                Matriz = rep(metadados$Matriz, length.out = nrow(dados_long)),
                Amostra = rep(metadados$Amostra, length.out = nrow(dados_long))
        )

ega_dados <- dados_long %>%
        mutate(Time = as.numeric(Time),
               Temperatura = 50 + (20 * Time),  
               Counts = as.numeric(Counts))

linhas_verticais <- linhas_verticais %>%
        arrange(Amostra, Tipo_Area, Estacao, Matriz, Temperatura) %>%
        group_by(Amostra, Tipo_Area, Estacao, Matriz) %>%
        mutate(Temp_inicio = Temperatura,
               Temp_fim = lead(Temperatura)) %>%
        filter(!is.na(Temp_fim)) %>%  
        ungroup()


ega_dados <- ega_dados %>%
        rename(Temp_medida = Temperatura) 

intervalo_valido <- linhas_verticais %>%
        group_by(Tipo_Area, Matriz, Estacao, Amostra) %>%
        summarise(
                temp_min = min(Temp_inicio),  
                temp_max = max(Temp_fim),    
                .groups = "drop"
        )

area_total_df <- ega_dados %>%
        inner_join(intervalo_valido, by = c("Tipo_Area", "Matriz", "Estacao", "Amostra")) %>%
        filter(Temp_medida >= temp_min & Temp_medida <= temp_max) %>%  # Limita a faixa de integração
        group_by(Tipo_Area, Matriz, Estacao, Amostra) %>%
        summarise(area_total = trapz(Temp_medida, Counts), .groups = "drop")


resultados <- ega_dados %>%
        inner_join(linhas_verticais, by = c("Amostra", "Tipo_Area", "Estacao", "Matriz"), relationship = "many-to-many") %>%
        filter(Temp_medida >= Temp_inicio & Temp_medida < Temp_fim) %>%
        group_by(Tipo_Area, Matriz, Estacao, Amostra, Temp_inicio, Temp_fim) %>%
        summarise(area_intervalo = trapz(Temp_medida, Counts), .groups = "drop") %>%
        left_join(area_total_df, by = c("Tipo_Area", "Matriz", "Estacao", "Amostra")) %>%
        left_join(pools_ega, by = c("Tipo_Area", "Matriz", "Estacao", "Amostra")) %>%
        mutate(
                massa_carbono = (area_intervalo / area_total) * C
        ) %>%
        select(Tipo_Area, Matriz, Estacao, Amostra, Temp_inicio, Temp_fim, area_intervalo, massa_carbono)
            
                    ################################################################


write_xlsx(resultados, "Pools_Carbono.xlsx")

shapiro.test(resultados$massa_carbono[resultados$Tipo_Area == "Permanent"& resultados$Matriz=="HA"])#non-normal
shapiro.test(resultados$massa_carbono[resultados$Tipo_Area == "Permanent"& resultados$Matriz=="Soil"])#non-normal
shapiro.test(resultados$massa_carbono[resultados$Tipo_Area == "Seasonal"& resultados$Matriz=="HA"])#normal
shapiro.test(resultados$massa_carbono[resultados$Tipo_Area == "Seasonal"& resultados$Matriz=="Soil"])#non-normal

# Kruskal-Wallis test - Tipo_Area e Matriz
kruskal.test(massa_carbono ~ Tipo_Area, data = resultados[resultados$Matriz=="HA", ]) #p>0.05
kruskal.test(massa_carbono ~ Tipo_Area, data = resultados[resultados$Matriz=="Soil", ]) #p>0.05

kruskal.test(massa_carbono ~ Estacao, data = resultados[resultados$Tipo_Area=="Permanent" & resultados$Matriz=="HA", ]) #p>0.05
kruskal.test(massa_carbono ~ Estacao, data = resultados[resultados$Tipo_Area=="Permanent" & resultados$Matriz=="Soil", ]) #p>0.05
kruskal.test(massa_carbono ~ Estacao, data = resultados[resultados$Tipo_Area=="Seasonal" & resultados$Matriz=="HA", ]) #p>0.05
kruskal.test(massa_carbono ~ Estacao, data = resultados[resultados$Tipo_Area=="Seasonal" & resultados$Matriz=="Soil", ]) #p>0.05

kruskal.test(massa_carbono ~ Matriz, data = resultados[resultados$Tipo_Area == "Permanent", ]) #p<0.05 HA>Soil
kruskal.test(massa_carbono ~ Matriz, data = resultados[resultados$Tipo_Area == "Seasonal", ]) #p>0.05
