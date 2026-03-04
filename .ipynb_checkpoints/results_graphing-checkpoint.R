library(tidyverse)
library(data.table)
library(patchwork)

# HARVARD PGP RESULTS

# Load PRS and phenotype data
heights <- fread("PRSice2/PRSice_phenotypes.txt")
ldpred2_results <- fread("LDPred2/ldpred2_PGP_results.tsv")
names(ldpred2_results) <- c("IID", "PRS_LDPred2")
prsice2_results <- fread("PRSice2/PRSice_PRS.best")
setnames(prsice2_results, old = "PRS", new = "PRS_PRSice2")

results <- merge(ldpred2_results, prsice2_results, by = "IID")
results <- merge(results, heights, by = "IID")
results <- results %>% select("IID", "PRS_PRSice2", "PRS_LDPred2", "height", "sex", "sex_num")
setorder(results, height)

# Perform regression to predict heights from PRS
ldpred2_model <- lm(height ~ PRS_LDPred2 + sex_num, data = results)
prsice2_model <- lm(height ~ PRS_PRSice2 + sex_num, data = results)
results$height_LDPred2 <- predict(ldpred2_model, newdata = results)
results$height_PRSice2 <- predict(prsice2_model, newdata = results)

# Compute coefficient of determination (R^2)
ldpred2_r2 <- results %>% group_by(sex_num) %>% summarize(r2 = round(summary(lm(height_LDPred2 ~ height)) $ r.squared, 3))
prsice2_r2 <- results %>% group_by(sex_num) %>% summarize(r2 = round(summary(lm(height_PRSice2 ~ height)) $ r.squared, 3))

# Plot reported height vs. predicted height
r2_label <- paste0("LDPred-2 (Female R²: ", ldpred2_r2$r2[1], ", Male R²: ", ldpred2_r2$r2[2], ")")
plot1 <- ggplot(results, aes(x = height, y = height_LDPred2, color = sex, group = sex)) + 
    geom_point(size = 3) + geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) + # xlim(60,77) + ylim(60,77) +
    labs(title = r2_label, x = "Reported Height (in)", y = "Predicted Height from PRS (in)")

r2_label <- paste0("PRSice-2 (Female R²: ", prsice2_r2$r2[1], ", Male R²: ", prsice2_r2$r2[2], ")")
plot2 <- ggplot(results, aes(x = height, y = height_PRSice2, color = sex, group = sex)) + 
    geom_point(size = 3) + geom_smooth(method = "lm", se = TRUE, fullrange = TRUE) +
    labs(title = r2_label, x = "Reported Height (in)", y = "Predicted Height from PRS (in)") +
    theme(legend.position = "none")

plot <- plot2 | plot1
ggsave("PGP_results.png", plot = plot, width = 12, height = 5, units = "in", dpi = 300)


# 1000 GENOMES RESULTS

# Load phenotype data
heights_1000G <- fread("LDPred2/1000Genomes/1000G_phenotypes.csv")
heights_1000G_male <- heights_1000G %>% select("superpopulation", "Mean_Male", "SD_Male")
names(heights_1000G_male) <- c("superpopulation", "mean", "sd")
heights_1000G_female <- heights_1000G %>% select("superpopulation", "Mean_Female", "SD_Female")
names(heights_1000G_female) <- c("superpopulation", "mean", "sd")

pop_1000G <- fread("LDPred2/1000Genomes/igsr_samples.tsv")
pop_1000G <- pop_1000G %>% select("Sample name", "Sex", "Superpopulation code")
names(pop_1000G) <- c("IID", "sex", "superpopulation")

# Load PDS data for LDPred-2
ldpred2_1000G_results <- fread("LDPred2/ldpred2_1000G_results.tsv")
names(ldpred2_1000G_results) <- c("IID", "PRS_LDPred2")
results <- merge(ldpred2_1000G_results, pop_1000G, by = "IID")
results <- results %>% filter(superpopulation != "EUR,AFR")

# Perform regression to predict heights from PRS 
# (split by superpopulation and sex)
results_comb <- NULL
for (pop in heights_1000G$superpopulation) {
    results_pop_female <- results %>% filter(superpopulation == pop) %>% filter(sex == "Female")
    results_pop_female$ZScore_LDPred2 <- (results_pop_female$PRS_LDPred2 - mean(results_pop_female$PRS_LDPred2)) / sd(results_pop_female$PRS_LDPred2)
    results_pop_female$height_LDPred2 <- heights_1000G_female[superpopulation == pop, "mean"][[1]] + results_pop_female$ZScore_LDPred2 * heights_1000G_female[superpopulation == pop, "sd"][[1]]
    results_pop_male <- results %>% filter(superpopulation == pop) %>% filter(sex == "Male")
    results_pop_male$ZScore_LDPred2 <- (results_pop_male$PRS_LDPred2 - mean(results_pop_male$PRS_LDPred2)) / sd(results_pop_male$PRS_LDPred2)
    results_pop_male$height_LDPred2 <- heights_1000G_male[superpopulation == pop, "mean"][[1]] + results_pop_male$ZScore_LDPred2 * heights_1000G_male[superpopulation == pop, "sd"][[1]]
    results_comb <- rbindlist(list(results_comb, results_pop_female, results_pop_male))
}

# Plot distributions of predicted heights
plot <- ggplot(results_comb, aes(x = height_LDPred2, fill = sex)) + xlim(50, 85) +
    geom_density(aes(color = sex), alpha = 0.5) + facet_wrap(~superpopulation, nrow = 1) +
    labs(title = "1000Genomes (LDPred-2)", x = "Predicted Height from PRS (in)") + theme(plot.title = element_text(hjust = 0.5))
ggsave("1000G_LDPred2_results.png", plot = plot, width = 20, height = 3, units = "in", dpi = 300)

# Load PDS data for PRSice-2
prsice2_1000G_results <- fread("PRSice2/PRSice_1000G.score")
names(prsice2_1000G_results) <- c("IID", "FID", "PRS_PRSice2")
results <- merge(prsice2_1000G_results, pop_1000G, by = "IID")
results <- results %>% filter(superpopulation != "EUR,AFR")

# Perform regression to predict heights from PRS 
# (split by superpopulation and sex)
results_comb <- NULL
for (pop in heights_1000G$superpopulation) {
    results_pop_female <- results %>% filter(superpopulation == pop) %>% filter(sex == "Female")
    results_pop_female$ZScore_PRSice2 <- (results_pop_female$PRS_PRSice2 - mean(results_pop_female$PRS_PRSice2)) / sd(results_pop_female$PRS_PRSice2)
    results_pop_female$height_PRSice2 <- heights_1000G_female[superpopulation == pop, "mean"][[1]] + results_pop_female$ZScore_PRSice2 * heights_1000G_female[superpopulation == pop, "sd"][[1]]
    results_pop_male <- results %>% filter(superpopulation == pop) %>% filter(sex == "Male")
    results_pop_male$ZScore_PRSice2 <- (results_pop_male$PRS_PRSice2 - mean(results_pop_male$PRS_PRSice2)) / sd(results_pop_male$PRS_PRSice2)
    results_pop_male$height_PRSice2 <- heights_1000G_male[superpopulation == pop, "mean"][[1]] + results_pop_male$ZScore_PRSice2 * heights_1000G_male[superpopulation == pop, "sd"][[1]]
    results_comb <- rbindlist(list(results_comb, results_pop_female, results_pop_male))
}

# Plot distributions of predicted heights
plot <- ggplot(results_comb, aes(x = height_PRSice2, fill = sex)) + xlim(50, 85) +
    geom_density(aes(color = sex), alpha = 0.5) + facet_wrap(~superpopulation, nrow = 1) +
    labs(title = "1000Genomes (PRSice-2)", x = "Predicted Height from PRS (in)") + theme(plot.title = element_text(hjust = 0.5))
ggsave("1000G_PRSice2_results.png", plot = plot, width = 20, height = 3, units = "in", dpi = 300)

