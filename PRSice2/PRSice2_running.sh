#!/bin/bash

Rscript PRSice.R \
    --prsice ./PRSice_linux \
    --base HarvardPGP/PRSice_PGP_sumstats.txt \
    --target HarvardPGP/PRSice_PGP_data \
    --snp rsid \
    --chr chr \
    --bp pos \
    --A1 a1 \
    --A2 a0 \
    --stat beta \
    --pvalue p \
    --pheno HarvardPGP/PRSice_PGP_phenotypes.txt \
    --pheno-col height \
    --cov HarvardPGP/PRSice_PGP_phenotypes.txt \
    --cov-col sex_num \
    --out PRSice_PRS

Rscript PRSice.R \
    --prsice ./PRSice_linux \
    --base 1000Genomes/PRSice_1000G_sumstats.txt \
    --target 1000Genomes/PRSice_1000G_data \
    --snp rsid \
    --chr chr \
    --bp pos \
    --A1 a1 \
    --A2 a0 \
    --stat beta \
    --pvalue p \
    --no-regress \
    --bar-levels 1 \
    --out 1000Genomes/PRSice_1000G