#!/bin/bash

Rscript PRSice.R \
    --prsice ./PRSice_linux \
    --base PRSice_sumstats.txt \
    --target PRSice_PGP_data \
    --snp rsid \
    --chr chr \
    --bp pos \
    --A1 a1 \
    --A2 a0 \
    --stat beta \
    --pvalue p \
    --pheno PRSice_phenotypes.txt \
    --pheno-col height \
    --cov PRSice_phenotypes.txt \
    --cov-col sex_num \
    --out PRSice_PRS