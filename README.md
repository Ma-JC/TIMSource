# TIMSource  
Complement inhibition reshapes the tumor microenvironment and enhances clinical efficacy in immunotherapy  
![image](overall.png)

This repository contains all the data analysis scripts utilized in our study entitled **"Complement inhibition reshapes the tumor microenvironment and enhances clinical efficacy in immunotherapy."** The project is organized into three core components: **Pan-cancer analysis**, **single-cell RNA sequencing (scRNA-seq) analysis**, and **TIMSource web server**. Below, we provide a detailed overview of each component along with guidance on how to navigate the repository.

## Repository Structure

### 1. Pan-cancer Analysis

The **Pan-cancer analysis** directory comprises scripts for the preprocessing and analysis of multiple datasets, including **ICB cohort datasets**, **TCGA datasets**, and **CPTAC datasets**. In this analysis, we developed ImmSNVscore to explore potential immunotherapy biomarkers and therapeutic targets, and further explore and verify the role of complement pathway mutations in immunotherapy through different datasets and methods:

- **Data_Preprocess**: Contains scripts for preprocessing the data from the aforementioned sources.
- **Data_Analysis**:
  - **Overview.ipynb**: Provides an summary of the data from the three sources.
  - **ImmuneSNVscore_SingleGene_OS+PFS.ipynb**: Assesses the prognostic value of individual gene mutations using the ImmSNVscore, specifically in the context of survival outcomes (Overall Survival and Progression-Free Survival) for patients receiving immune checkpoint blockade (ICB) therapies.
  - **ImmuneSNVscore_Complement_OS+PFS.ipynb**: Evaluates mutations within pathways/biological processes, with a specific focus on complement-related pathways, and their association with patient survival under ICB therapy. Particular emphasis is placed on the analysis and validation of the GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION pathway.
  - **Pan-Cancer Analysis.ipynb**: Investigates expression profile alterations and functional pathway changes in tumors harboring complement mutations across different cancer types.

### 2. scRNA-seq Analysis

This directory contains scripts essential for the analysis of single-cell RNA sequencing (scRNA-seq) data, which is a critical component of our study. These analyses enable us to dissect the tumor microenvironment at a single-cell resolution:

- **Mouse single-cell analysis**: Conducts comprehensive analyses of single-cell transcriptomic data derived from AKR and CT26 murine models. This includes cell type identification, differential expression analysis, and pathway enrichment.
- **Human single-cell analysis**: Involves whole-exome sequencing and single-cell transcriptomic analyses of human esophageal cancer samples, focusing on cellular heterogeneity and the impact of complement inhibition on the tumor microenvironment.

### 3. TIMSource Web Server

The **TIMSource** directory contains the source code for the TIMSource database, an interactive web server built on the Shiny framework, designed to facilitate exploration and visualization of our findings:

- **app.R**: Defines the user interface (UI) layout and the core server-side logic.
- **script**: Contains the detailed server logic and data processing scripts required for the various modules of the web application. Each module is designed to provide users with interactive access to key datasets and analyses performed in the study.
