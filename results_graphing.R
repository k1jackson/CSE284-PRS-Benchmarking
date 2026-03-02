library(tidyverse)
library(data.table)
library(patchwork)

heights <- fread("PRSice2/PRSice_phenotypes.txt")
ldpred2_results <- fread("LDPred2/ldpred_results.tsv")
names(ldpred2_results) <- c("IID", "PRS_LDPred2")
prsice2_results <- fread("PRSice2/PRSice_PRS.best")
setnames(prsice2_results, old = "PRS", new = "PRS_PRSice2")

results <- merge(ldpred2_results, prsice2_results, by = "IID")
results <- merge(results, heights, by = "IID")
results <- results %>% select("IID", "PRS_PRSice2", "PRS_LDPred2", "height", "sex", "sex_num")
setorder(results, height)

ldpred2_model <- lm(height ~ PRS_LDPred2 + sex_num, data = results)
prsice2_model <- lm(height ~ PRS_PRSice2 + sex_num, data = results)
results$height_LDPred2 <- predict(ldpred2_model, newdata = results)
results$height_PRSice2 <- predict(prsice2_model, newdata = results)

# results <- results %>% group_by(sex) %>% mutate(pred_height = scale(PRS) * sd(height) * 2 + mean(height)) %>% ungroup()
ldpred2_r2 <- results %>% group_by(sex_num) %>% summarize(r2 = round(summary(lm(height_LDPred2 ~ height)) $ r.squared, 3))
prsice2_r2 <- results %>% group_by(sex_num) %>% summarize(r2 = round(summary(lm(height_PRSice2 ~ height)) $ r.squared, 3))

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
ggsave("results.png", plot = plot, width = 12, height = 5, units = "in", dpi = 300)