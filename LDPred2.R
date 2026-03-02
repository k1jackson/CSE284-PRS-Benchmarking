library(bigsnpr)
library(bigreadr)
library(tidyverse)
library(data.table)

# load Harvard PGP data
rds_path <- snp_readBed("PGP_filtered_data.bed") #"clean_plink_data.bed")
obj.bigSNP <- snp_attach(rds_path)
genotypes <- snp_fastImputeSimple(obj.bigSNP$genotypes, method = "mean0", ncores = nb_cores())
SNPs <- obj.bigSNP$map
names(SNPs) <- c("chr", "rsid", "dist", "pos", "a1", "a0")
fam.order <- as.data.table(obj.bigSNP$fam)
setnames(fam.order, c("family.ID", "sample.ID"), c("FID", "IID"))

# load summary statistics
sum_stats <- fread2("GIANT_GWAS_summary_stats.tsv")
names(sum_stats) <- c("id", "rsid", "chr", "pos", "a1", "a0", "EAF", "beta", "beta_se", "p", "n_eff")
sum_stats <- sum_stats %>% filter(chr %in% 1:22)

# filter SNPs
hapmap3 <- readRDS("map_hm3_plus.rds")
info_SNPs <- snp_match(sum_stats, hapmap3, join_by_pos = TRUE)
setnames(info_SNPs, old = "_NUM_ID_", new = "hapmap3_index")
info_SNPs <- snp_match(info_SNPs, SNPs, join_by_pos = FALSE, match.min.prop = 0)
setnames(info_SNPs, old = "_NUM_ID_", new = "plink_index")

# construct LD matrix
SNPs_pos <- snp_asGeneticPos(info_SNPs$chr, info_SNPs$pos, dir = ".")
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
corr <- NULL
ld <- NULL
for (chr in 1:22) {
    ind.chr <- which(info_SNPs$chr == chr)
    offset <- which(hapmap3$chr == chr)
    ind.chr2 <- info_SNPs$hapmap3_index[ind.chr] - min(offset) + 1
    file_ld <- paste0("ldref_hm3_plus/LD_with_blocks_chr", chr, ".rds")
    corr0 <- readRDS(file_ld)[ind.chr2, ind.chr2]
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

# LD score regression
df_beta <- info_SNPs[,c("beta", "beta_se", "n_eff", "hapmap3_index")]
ldsc <- snp_ldsc(ld, length(ld), chi2 = (df_beta$beta / df_beta$beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL)
h2_est <- ldsc[["h2"]]

# compute PRS
final_beta <- snp_ldpred2_inf(corr, info_SNPs, h2 = h2_est)
prs <- big_prodVec(genotypes, final_beta, ind.row = 1:nrow(genotypes), ind.col = info_SNPs$plink_index)
results <- data.frame(sample = sub("_genome", "", obj.bigSNP$fam$family.ID), PRS = prs)
write.table(results, "ldpred2_results.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

# compare heights
heights <- fread("PGP_phenotypes.csv")#"height_data.csv")
names(heights) <- c("sample", "height", "sex")
results <- merge(results, heights, by = "sample")
setorder(results, height)

results <- results %>% group_by(sex) %>% mutate(pred_height = scale(PRS) * sd(height) * 2 + mean(height)) %>% ungroup()
r2_values <- results %>% group_by(sex) %>% summarize(r2 = round(summary(lm(pred_height ~ height)) $ r.squared, 3))
r2_label <- paste0("Female R²: ", r2_values$r2[r2_values$sex == "Female"], ", Male R²: ", r2_values$r2[r2_values$sex == "Male"])

plot <- ggplot(results, aes(x = height, y = pred_height, color = sex, group = sex)) + 
    geom_point(size = 3) + geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
    labs(title = r2_label, x = "Reported Height (in)", y = "Predicted Height from PRS (in)")

ggsave("ldpred2_results.png", plot = plot, width = 6, height = 4, units = "in", dpi = 300)