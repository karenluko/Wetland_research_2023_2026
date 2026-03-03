setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Pós doc/Areas umidas/Wetland_research_GitHub/scripts_data")
dados.na<-read.csv2("Biomarkers.csv", stringsAsFactors = T, na.strings = "NA")

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggh4x)
library(ggpattern)

dados_long <- dados.na %>%
        pivot_longer(cols = -c(Matrix, Sample, Type, Season), 
                     names_to = "Alkane_number", 
                     values_to = "Abundance") %>%
        mutate(Season = factor(Season, levels = c("Rainy","Dry")),
               Type = factor(Type, levels = c("Permanent","Seasonal")),
               Matrix = factor(Matrix, levels = c("Soil","HA")),
               Sample = factor(Sample, levels = c("P01","P02", "P03")),
               Alkane_number = as.numeric(gsub("\\D", "", Alkane_number)))  # Extrai número do alcano

names(dados_long)[names(dados_long) == "Type"] <- "Wetland type"
names(dados_long)[names(dados_long) == "Sample"] <- "Sampling point"

cores_personalizadas <- c("Rainy" = "plum4", "Dry" = "peachpuff3")

plot <- ggplot(dados_long, aes(x = Alkane_number, y = Abundance, fill = Season, pattern = Matrix, colour = Matrix )) +
        geom_bar_pattern(stat = "identity", 
                         position = position_dodge2(preserve = "single", padding = 0.2), 
                         width = 0.8,
                         pattern_density = 0.02,
                         pattern_spacing = 0.05,
                         pattern_fill = "black",
                         pattern_alpha = 0.8,
                         linewidth=0.5)+
        facet_nested(`Wetland type`+ Season ~ `Sampling point` + Matrix, scales = "free", space = "free_x") +
        scale_fill_manual(values = cores_personalizadas) +
        scale_pattern_manual(values = c("Soil" = "circle", "HA" = "none")) +
        scale_color_manual(values = c("Soil" = "black", "HA" = "black")) + 
        labs(x = expression(italic(n)*"-Alkane number (C"[n]*")"),
        y = expression("Abundance of each " * italic(n) * "-alkane relative to the sum of " * italic(n) * "-alkanes per sample (%)"),
             fill = "Season",
             pattern = "Matrix") +
        scale_x_continuous(breaks = seq(min(dados_long$Alkane_number), max(dados_long$Alkane_number), by = 2)) +
        theme_light() +
        theme(legend.position = "bottom",
              strip.text = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(angle = 90, hjust = 1),
              panel.spacing.x = unit(1, "mm"))+
        guides(fill = "none", pattern = "none", color = "none")

print(plot)
ggsave("histogram_alkanes_revised.tiff", plot, width = 12, height = 8, dpi = 300)

###################################################

# Filtra alcanos de C24 a C35
dados_long_filtrado <- subset(dados_long, Alkane_number >= 24 & Alkane_number <= 35)

plot <- ggplot(dados_long_filtrado, aes(x = Alkane_number, y = Abundance, fill = Season, pattern = Matrix, colour = Matrix)) +
        geom_bar_pattern(stat = "identity", 
                         position = position_dodge2(preserve = "single", padding = 0.2), 
                         width = 1,
                         pattern_density = 0.1,
                         pattern_spacing = 0.03,
                         pattern_fill = "black",
                         pattern_alpha = 0.8,
                         linewidth = 3) +
        facet_nested(`Wetland type` + Season ~ `Sampling point` + Matrix, scales = "free", space = "free_x") +
        scale_fill_manual(values = cores_personalizadas) +
        scale_pattern_manual(values = c("Soil" = "circle", "HA" = "none")) + 
        scale_color_manual(values = c("Soil" = "black", "HA" = "black")) + 
        labs(x = "Alkane number (Cn)", 
             y = "Abundance of each alkane relative to the sum of alkanes per sample (%)",
             fill = "Season",
             pattern = "Matrix") +
        scale_x_continuous(breaks = seq(24, 35, by = 2)) +
        theme_light() +
        theme(legend.position = "bottom",
              strip.text = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(angle = 90, hjust = 1),
              panel.spacing.x = unit(1, "mm")) +
        guides(fill = "none", pattern = "none", color = "none")  

plot

