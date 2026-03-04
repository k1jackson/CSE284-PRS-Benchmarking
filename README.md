# CSE284-PRS-Benchmarking



Title: Comparative Analysis of PRSice2 and LDPred to Compute PRS for Human Height   
Partners: Kate Jackson (A16781414) and Fatemeh Heydari (A69027046)





Introduction

A polygenic risk score (PRS) is a numerical metric to estimate an individual's genetic predisposition for a complex trait or disease. PRSs are computed by summing the effects of thousands of small genetic variants (SNPs) across the entire genome, to quantify a personalized risk score that is independent of environmental risk. Various methods have been developed to calculate PRS. In this project, we compare PRSice-2 and LDpred2, which utilize fundamentally different algorithmic approaches, to predict height, a classic polygenic trait.

PRSice-2 (Choi & O’Reilly, 2019) is the gold standard for the Clumping and Thresholding (C+T) method. It simplifies the genetic architecture by removing correlated SNPs (clumping) and including only those that meet a specific p-value threshold. This method is computationally efficient and interpretable, making it ideal for a preliminary analysis. In contrast, LDpred2 (Privé et al., 2020) utilizes a Bayesian approach to model the LD structure across the entire genome. Rather than discarding SNPs, it re-weights them based on a prior distribution of effect sizes and an LD reference panel. This typically results in higher predictive accuracy, especially for highly polygenic traits where many small-effect variants are masked by strict p-value thresholds in C+T. In this project we will assess and compare the performance of these two methods.

Methods

Data
As most PRS methods require GWAS summary statistics, we obtained these from Yengo et al. (2022), a comprehensive study conducted across diverse populations to characterize genetic variants associated with human height. From this dataset, we identified relevant SNPs and extracted critical summary statistics, such as effect sizes and p-values. To obtain individual-level genotype and phenotype (height) data, we utilized the Harvard Personal Genome Project (PGP), a public repository where participants share their genomic data for scientific research. From this database, we extracted 23andMe data for 155 individuals.


Calculating Polygenic Risk Scores
To prepare the 23andMe data as input for PRS methods, we used PLINK to generate binary files (BED/BIM/FAM) for each individual. During the initial merging process, PLINK generated a .missnp file, identifying variants with allelic inconsistencies (such as multiallelic sites) that prevented a successful merge. We resolved this by excluding these problematic SNPs and re-running the merge to create a unified dataset.
Height data was standardized to inches and because sex is a major determinant of height, it was included as a covariate. Finally, after calculating the scores using PRSice-2 and LDpred2, we fitted a linear regression model using both the PRS and sex as predictors to evaluate the model's accuracy.

Results


Next step
We will conduct this analysis using the 1000 Genomes dataset. Although height data is unavailable for this cohort, we will see what insights can be gained.
