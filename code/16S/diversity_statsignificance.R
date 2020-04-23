alpha_diversity <- read.delim("data/16S/mothur_output/alpha_diversity/groups.summary", stringsAsFactors = F)

alpha_diversity$group[1:30] <- "Adenoma"
alpha_diversity$group[31:60] <- "Cancer"
alpha_diversity$group[61:90] <- "Healthy"

library(vegan)

#anosim(x=alpha_diversity[,13],grouping = alpha_diversity$group)
#using anova over anosim because couldn't specify distance as shannon and default is bray curtis?
anova_result <- aov(shannon ~ group, data = alpha_diversity)
summary(anova_result)

beta_diversity <- read.delim("data/16S/mothur_output/beta_diversity/braycurtis.0.03.pcoa.axes", stringsAsFactors = F)

beta_diversity$group[1:30] <- "Adenoma"
beta_diversity$group[31:60] <- "Cancer"
beta_diversity$group[61:90] <- "Healthy"

#anosim(x=beta_diversity[,2:90], grouping = beta_diversity$group)
#doesn't work with negative values
anova_result2a <- aov(axis1~group, data = beta_diversity)
summary(anova_result2a)
anova_result2b <- aov(axis2~group, data = beta_diversity)
summary(anova_result2b)
#only looks at statistical significance for two axes that represent large fraction of total variance