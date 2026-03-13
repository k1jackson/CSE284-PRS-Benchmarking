# Comparative Analysis of PRSice2 and LDPred to Compute Polygenic Risk Scores for Human Height   
**Authors:** Kate Jackson (A16781414) and Fatemeh Heydari (A69027046)

## Respository Structure
- **results_graphing.R:** RScript used to produce all output plots 
- **PGP_results.png:** Plot of reported vs. predicted heights for Harvard PGP data with LDPred-2 and PRSice-2
- **1000G_LDPred2_results.png:** Histogram of predicted heights for 1000 Genomes data with LDPred-2
- **1000G_PRSice2_results.png:** Histogram of predicted heights for 1000 Genomes data with PRSice-2
- **LDPred2:** Data and scripts to run LDPred-2
  - **1000Genomes:** Selected (<100MB) 1000 Genomes data
    - **1000G_phenotypes.tsv:** Height summary data by subpopulation from Yengo et al. (2022)
    - **igsr_samples.tsv** Supopulation data from 1000 Genomes
  - **HarvardPGP:** Selected (<100MB) Harvard PGP data
    - **PGP_data_processing.sh:** Bash script with commands used to preprocess PGP data
    - **PGP_phenotypes.csv:** Height data from Harvard PGP survey
  - **LDPred2_running.R**: RScript used to run LDPred-2
  - **ldpred2_1000G_results.tsv:** PRS for 1000 Genomes data with LDPred-2
  - **ldpred2_PGP_results.tsv:** PRS for Harvard PGP data with LDPred-2
- **PRSice2:** Data and scripts to run PRSice-2
  - **PRSice2_running.R**: Bash script used to run PRSice-2
  - **1000Genomes:** Selected (<100MB) 1000 Genomes data and outputs
    - **PRSice_1000G.score:** PRS for 1000 Genomes with PRSice-2
  - **HarvardPGP:** Selected (<100MB) Harvard PGP data and outputs
    - **PRSice_PGP.best:** PRS for Harvard PGP with PRSice-2
