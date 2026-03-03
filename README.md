# CSE284-PRS-Benchmarking



Title: Comparative Analysis of PRSice2 and LDPred to Compute PRS for Human Height   
Partners: Kate Jackson (A16781414) and Fatemeh Heydari (A69027046)





Introduction

A polygenic risk score (PRS) is a numerical calculation that estimates an individual's genetic predisposition to a specific disease or trait, such as heart disease or diabetes. By summing the effects of thousands of small genetic variants (SNPs) across the entire genome, it provides a personalized, lifetime risk score that is independent of environmental risk. Various methods have been developed to calculate PRS; in this project, we compare PRSice-2 and LDpred2, which utilize fundamentally different algorithmic approaches, to predict height, a classic polygenic trait.

PRSice-2 (Choi & O’Reilly, 2019) is the gold standard for the Clumping and Thresholding (C+T) method. It simplifies the genetic architecture by removing correlated SNPs (clumping) and including only those that meet a specific p-value threshold. This method is computationally efficient and highly transparent, making it ideal for a preliminary analysis. In contrast, LDpred2 (Privé et al., 2020) utilizes a Bayesian approach to model the LD structure across the entire genome. Rather than discarding SNPs, it re-weights them based on a prior distribution of effect sizes and an LD reference panel. This typically results in higher predictive accuracy, especially for highly polygenic traits where many small-effect variants are masked by strict p-value thresholds in C+T. In this project we will assess and compare the performance of these two methods.


