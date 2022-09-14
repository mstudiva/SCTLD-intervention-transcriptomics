#### packages ####

library(ggplot2)
library(ggpubr)
library(rcompanion)
library(MASS)
library(stringr)
library(survival)
library(survminer)


#### data import ####

transmission <- read.csv("sctld omics transmission.csv", head=T)
str(transmission)
head(transmission)

# subsets the data frame to remove healthy sediment (since no disease observed)
trans <- subset(trans, treatment!="control")
trans

#### data normality/transformation ####

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test(trans$lesion)
# not normal

# BoxCox transformation finds best exponent value for data transformation
trans_bc_days <- boxcox(trans$lesion ~ trans$species_treatment)
trans_lambda_days <- trans_bc_days$x[which.max(trans_bc_days$y)]

# Shapiro test, p-values below 0.05 indicate violations of normality assumptions
shapiro.test((trans$lesion^trans_lambda_days-1)/trans_lambda_days)
# now normal

pdf("sctld omics normality.pdf")
par(mfrow=c(2,2))
hist(trans$lesion)
hist((trans$lesion^trans_lambda_days-1)/trans_lambda_days)
qqnorm(trans$lesion)
qqline(trans$lesion)
qqnorm((trans$lesion^trans_lambda_days-1)/trans_lambda_days)
qqline((trans$lesion^trans_lambda_days-1)/trans_lambda_days)
dev.off()
# continuing with analysis on transformed data


#### statistical tests ####

# ANOVA on Box-Cox transformed data
trans_anova_days_bc <- aov(((lesion^trans_lambda_days-1)/trans_lambda_days) ~ species, data=trans)
summary(trans_anova_days_bc)

# Q-Q plot for transformed data
pdf("sctld omics normality ANOVA.pdf")
par(mfrow=c(2,2))
qqnorm(trans_anova_days_bc$residuals, main="days transformed"); qqline(trans_anova_days_bc$residuals)
dev.off()

# Tukey post hoc tests
trans_tukey_days <- TukeyHSD(trans_anova_days_bc)
trans_tukey_days

# creating dataframes of the pairwise comparisons needed for plots and doing a bit of table reformatting
trans_letters_days <- data.frame(trans_tukey_days$`species`)
trans_letters_days$Var <- rownames(trans_letters_days)
names(trans_letters_days)[5] <- "comparison"
trans_letters_days$comparison = str_replace_all(trans_letters_days$comparison,":","_")
trans_letters_days$p.adj[is.na(trans_letters_days$p.adj)] <- 1
trans_letters_days

# creates compact letter display of significant pairwise differences for figure
trans_cld_days <- cldList(p.adj ~ comparison, data = trans_letters_days, threshold = 0.05)
trans_cld_days
trans_cld_days$Group <- c("Mc","Of")


#### time to transmission figures ####

# keeps treatment order as imported
trans$species=factor(trans$species, levels=unique(trans$species)) 

# boxplots comparing time to transmission among species/treatments
transmission_days <-
  ggboxplot(
    trans,
    x = "species",
    y = "lesion",
    color = "grey30",
    fill = "species",
    add = "jitter",
    add.params = list(size = 1, jitter = 0.25),
    width = 0.7,
    size = 0.5
  ) + labs(x = "Species",
           y = "Time to Transmission (d)",
           fill = 'Species') + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),  legend.position = "right")+
  geom_text(data=trans_cld_days, aes(x = Group, y=0, label=Letter)) 
transmission_days

ggsave("sctld omics transmission.pdf", plot= transmission_days, width=3, height=4, units="in", dpi=300)


#### trans transmission rate figure ####

trans_rate <- read.csv("sctld omics rate.csv", head = T)
trans_rate
# trans_rate$treatment=factor(trans_rate$treatment, levels=unique(trans_rate$treatment)) 
omics_transrate <- ggplot(trans_rate, aes(fill=forcats::fct_rev(condition), y=rate, x=species)) + 
  # scale_fill_manual(values=c("green", "red")) +
  geom_col(width = 0.5) +
  theme_bw() 
omics_transrate

ggsave("sctld omics rate.pdf", plot= omics_transrate, width=4.25, height=2, units="in", dpi=300)


#### trans survivorship ####

# create survival object (using the Kaplan-Meier method)
omics_surv <- Surv(time = trans$days, event = trans$status)
omics_surv

# run survival model
omics_fit <- survfit(omics_surv ~ species, data = trans)
summary(omics_fit)

# Kaplan-Meier plot 
omics_survival<-ggsurvplot(omics_fit, data = trans, pval = TRUE, xlab="Days", ylab="Health probability",
                             conf.int = T, risk.table=T, 
                             break.time.by=2, xlim=c(0,14), risk.table.y.text = FALSE) + ggtitle("Species") 
omics_survival

# hazard ratio by species
omics_haz <- coxph(omics_surv ~ species, data = trans)
summary(omics_haz)

# plot
omics_hazard <- ggforest(omics_haz, data = trans)


#### survivorship figure ####

omics_survival_multiplot<-ggarrange(omics_survival$plot,
                                 omics_survival$table, 
                                 omics_hazard,
                                 heights = c(2, 0.6, 0.75),
                                 ncol = 1, nrow = 3)
omics_survival_multiplot

ggsave("sctld omics survivorship.pdf", omics_survival_multiplot, width=8, height=8,dpi = 300)

