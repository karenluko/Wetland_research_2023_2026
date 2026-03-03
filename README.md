# Wetland_research_2022_2026
Raw and processed data and scripts related to research of two tropical wetlands

Project Title:
Composition of humic substances and their role in electron transfer in wetland areas of the São Paulo Peripheral Depression

Researcher:
Dr. Karen Luko-Sulato

Supervisor: 
Prof Dr Vania Silvia Rosolen

Institution:
São Paulo State University (UNESP)

Funding Agency:
São Paulo Research Foundation(FAPESP) - Processo no 2021/06332-1 and 2023/15396-9

Project Period:
01/11/2022 to 19/02/2026 (DD/MM/YY)

----------------------------------------------------------------------------------------------------------
1. GENERAL DESCRIPTION
----------------------------------------------------------------------------------------------------------

This dataset refers to the environmental, chemical, microbiological, and statistical analyses carried out within the scope of the project involving studies of soils and greenhouse gases in wetlands under different hydrological regimes.
The data include raw and processed information used in the preparation of scientific articles and technical reports.

----------------------------------------------------------------------------------------------------------
2. DATA STRUCTURE
----------------------------------------------------------------------------------------------------------

The repository is organized as follows:

/raw_data/
/EGA/ → .csv files from evolved gas cromatography
/Pyrolisys/ → raw data from the software for analytical pyrolysis (humic acid and soil)
/UV_Vis/ → Spectral scans (absorbance × wavelength)
 	TGs_Karen_2024.pptx → .pptx file with termogravimetric curve from thermogravimetry analysis

/processed/
/GC / → .xlsx files from gas chromatography (FID-GC) (CO2 and CH4)
Data_TG_Karen_2024.xlsx → .xlsx file with dataset from TGA

Metabolism_Wetlands.xlsx → spreadsheet with putative metabolism information from 16S sequencing - Functional Annotation of Prokaryotic Taxa (FAPROTAX), program version 1.2.7 (Louca et al., 2016).

Resultados_Kruskal_Dunn_formatada.xlsx → spreadsheet with data from statistical anlysis

Taxonomy_Wetlands_Arq.xlsx → spreadsheet with data of the relative abundancies of Archaea found in the soil samples from the wetlands by Phylum, Class, Order, Family and genus*
Taxonomy_Wetlands_Bact.xlsx → spreadsheet with data of the relative abundancies of Bacteria in the soil samples from the wetlands by Phylum, Class, Order, Family and genus*

* Raw data processing: QIIME2 version 2023.9 (Bolyen et al., 2019). 
Generation of table of ASVs:DADA2 pipeline (Callahan et al., 2016).
Taxonomic identification: SILVA database (97% cutoff threshold)

/scripts_data/ → folder with R scripts and their respective required dataset.

#Analytical pyrolisys:
Histograms.R
Biomarkers.csv

#Evolved gas analysis
EGA.R
EGA_dados.csv
Linhas_verticais_temperatura.csv

#Characterization of the microorganisms
16S_revised_paper.R
Dados_16S_depth.csv
Dados_16S

Rarefacao_alpha_diversity.R
Feature-table_Wetland 1.tsv

#Characterization of the humic fractions

Functional_groups.R
Grupos.csv

Quinonas.R
/Reducao/

tabela_mestre.csv → experimental dataset used to run the R scripts

UV.R
/UV-vis/

#Integrated data from other analysis. Includes: characterization of humic fraction, GHG fluxes, quantification of inorganic electron acceptors and soil physical-chemical parameters.

Integrated_data.R
Tabela_mestre.csv

In some cases, measurements were performed manually (UV at selected absorbance, pH, EH, soil moisture measures, total carbon) and later digitized and thereby raw data from the analyzer is not available.

----------------------------------------------------------------------------------------------------------
3. ANALYTICAL METHODOLOGY
----------------------------------------------------------------------------------------------------------

Main techniques employed:

- Gas Chromatography (FID-GC)
- Pyrolysis-GC/MS
- Evolved gas analysis
- Thermogravimetric analysis
- UV-Vis Spectroscopy
- 16S rRNA gene sequencing

Full methodological details are described in the articles associated
with the project.

----------------------------------------------------------------------------------------------------------
4. DATA PROCESSING AND QUALITY CONTROL
----------------------------------------------------------------------------------------------------------

The data underwent:

- Consistency checking
- Removal of invalid values
- Calculation of means and standard deviations

----------------------------------------------------------------------------------------------------------
5. 16S SEQUENCING DATA
----------------------------------------------------------------------------------------------------------

Sequencing data are available in the NCBI repository under BioProject ID PRJNA1254781

----------------------------------------------------------------------------------------------------------
6. SOFTWARE AND REPRODUCIBILITY
----------------------------------------------------------------------------------------------------------

Statistical analyses were performed using:

R (version 2026.01.1+403)

Main packages: 

AICcmodavg
agricolae
broom
corrplot
cowplot
dplyr
FactoMineR
factoextra
FSA
ggh4x
ggforce
ggfortify
ggpattern
ggplot2
ggpubr
ggthemes
janitor
pracma
purrr
readr
reshape2
rstatix
tidyr
tidyverse
wesanderson
writexlR 

----------------------------------------------------------------------------------------------------------
7. DATA ACCESS
----------------------------------------------------------------------------------------------------------

All materials are provided for scientific research purposes.

For further information or questions, please contact:

karen.luko@unesp.br

----------------------------------------------------------------------------------------------------------
8. CITATION
----------------------------------------------------------------------------------------------------------

When using these data, please cite:

