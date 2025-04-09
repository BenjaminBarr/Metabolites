library(here) # here makes a project transportable
library(janitor) # clean_names
library(readxl) # read excel, duh!
library(data.table) # magical data frames
library(magrittr) # pipes
library(stringr) # string functions
library(forcats) # factor functions

# analysis packages
library(emmeans) # the workhorse for inference
library(nlme) # gls and some lmm
library(lme4) # linear mixed models
library(lmerTest) # linear mixed model inference
library(afex) # ANOVA linear models
library(glmmTMB) # generalized linear models
library(MASS) # negative binomial and some other functions
library(car) # model checking and ANOVA
library(DHARMa) # model checking

# graphing packages
library(ggsci) # color palettes
library(ggpubr) # publication quality plots
library(ggforce) # better jitter
library(cowplot) # combine plots
library(knitr) # kable tables
library(kableExtra) # kable_styling tables

library(insight)
library(lazyWeave)

# Data Read In

library(ggplot2)
library(reshape2)
library(matrixStats)
library(dplyr)
library(AICcmodavg)
library(tidyverse)
library(rstatix)
library(ggdark)
library(multcompView)
library(RColorBrewer)

library(cowplot)
library(gridExtra)

getwd()

gluc <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Neg.xlsx", sheet = "Glucose")
gluc <- data.frame(gluc)

lactate <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Neg.xlsx", sheet = "Lactate")
lactate <- data.frame(lactate)

ala <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Pos.xlsx", sheet = "Alanine Pos")
ala <- data.frame(ala)

gln <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Pos.xlsx", sheet = "Glutamine Pos")
gln <- data.frame(gln)

carn <- read_xlsx("./Metabolites/Metabolites Data/Metabolites Pos.xlsx", sheet = "Carnitine Pos")
carn <- data.frame(carn)

gluc <- melt(gluc)
lactate <- melt(lactate)
ala <- melt(ala)
gln <- melt(gln)
carn <- melt(carn)


met <- cbind(gluc, lactate[3], ala[3], gln[3], carn[3])
names(met) <- c("label", "variable", "glucose", "lactate", 
                "alanine", "glutamine", "carnitine")

met.labs <- str_split_fixed(met$label, " ", 3)
met.labs <- data.frame(met.labs)
names(met.labs) <- c("age", "group", "sex")
met.labs

met <- cbind(met, met.labs)
met <- na.omit(met)
met <- met[, -2]
met

met$fat <- "HF"
met$fat[which(met$group %in% c("G1", "G2"))] <- "C"
met

met$protein <- "C"
met$protein[which(met$group %in% c("G3", "G4"))] <- "Beef"
met

met$trt <- "yes"
met$trt[which(met$group %in% c("G2", "G4", "G6"))] <- "no"
met

f.met <- met[which(met$sex == "F"),] # Females
hff.met <- f.met[which(f.met$group %in% c("G3", "G4", "G5", "G6")),] #HF only

hff.met
## Females ##
#Glucose First#

ggqqplot(hff.met$glucose)

hf.glucose.f.lm <- lm(glucose ~ trt * protein * age, data = hff.met)
anova(hf.glucose.f.lm)
summary(hf.glucose.f.lm)


f.glucose.em <- emmeans(hf.glucose.f.lm, specs = c("trt", "protein", "age"))
f.glucose.em

f.glucose.dt <- data.table(summary(f.glucose.em))
f.glucose.dt


# Contrast "pairwise comparisions"
f.glucose.pairs <- contrast(f.glucose.em,
                            method = "revpairwise",
                            #simple = "each",
                            combine = T,
                            adjust = "sidak") %>%
  summary(infer = T, )

f.glucose.pairs.dt <- data.table(f.glucose.pairs)
f.glucose.pairs.dt[which(f.glucose.pairs.dt$p.value <= 0.05)]

# Coordinates for significance indicators
comp <- f.glucose.pairs.dt$contrast

g <- gsub("no Beef", "group 4", as.character(comp))
g <- gsub("yes Beef", "group 3", as.character(g))
g <- gsub("no C", "group 6", as.character(g))
g <- gsub("yes C", "group 5", as.character(g))

g <- as.factor(g)
age <- gsub("group [3-6] ", "", as.character(g))
age1 <- gsub(" - H[1-3]", "", as.character(age))
age2 <- gsub("H[1-3] - ", "", as.character(age))
g1 <- gsub(" H[1-3] - group [3-6] H[1-3]", "", as.character(g))
g2 <- gsub("group [3-6] H[1-3] - ", "", as.character(g))
g2 <- gsub(" H[1-3]", "", as.character(g2))

f.glucose.pairs.dt$group1 <- g1
f.glucose.pairs.dt$group2 <- g2
f.glucose.pairs.dt$age1 <- age1
f.glucose.pairs.dt$age2 <- age2

# Nice p-values for graphs
f.glucose.pairs.dt[, p_rounded := p_round(p.value,
                                          digits = 2)]
f.glucose.pairs.dt[, p_pretty := p_format(p_rounded,
                                          digits = 2,
                                          accuracy = 1e-04,
                                          add.p = TRUE)]
f.glucose.pairs.dt

# Graphing labels ()
f.glucose.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5")


f.glucose.pairs.sig <- f.glucose.pairs.dt[which(f.glucose.pairs.dt$p.value <= 0.05),]
names(f.glucose.pairs.sig)

# Identify Significance #
f.glucose.pairs.sig

# Within Group Indications #
f.glucose.pairs.sig[which(f.glucose.pairs.sig$group1 == f.glucose.pairs.sig$group2),]

f.glucose.pairs.sig <- f.glucose.pairs.sig[,-c(1,5:7)]
f.glucose.pairs.sig

### Customize Here ###
fg.lab <- paste("group", rep(c(3:6), each = 3))
fg.lab <- data.frame(fg.lab)
fg.lab$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
fg.lab$lab <- c("A", "A", "A",
                "A", "A", "A",
                "A", "A", "A",
                "A", "A", "A")
names(fg.lab)[1] <- "group1"
fg.lab$group2 <- NA
fg.lab$group <- fg.lab$group1

# Customize for each set of comp
alab <- c(`H1` = "6 Months", `H2` = "12 Months", 
          `H3` = "18 Months")

glab <- c(`group 3` = "HFBN",
          `group 4` = "HFB",
          `group 5` = "HFCN",
          `group 6` = "HFC")

names(fg.lab)[1:2] <- c("group1", "age")

#Plot
f.glu.sum <- hff.met%>% group_by(group, age) %>%
  summarise(emmean = mean(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

fg.lab$y <- c(0, f.glu.sum$emmean + f.glu.sum$SE + 1)
g3f.h1.sum <- data.frame(group = "G3", 
                         age = "H1", 
                         emmean = NA,
                         SE = NA,
                         n = NA)

f.glu.sum <- rbind(g3f.h1.sum, f.glu.sum)

f.glu.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
f.glu.sum$age <- as.factor(f.glu.sum$age)
f.glu.sum <- f.glu.sum[order(f.glu.sum$age),]
f.glucose.dt <- f.glucose.dt[order(f.glucose.dt$age,f.glucose.dt$group),]
f.glu.sum <- f.glu.sum[order(f.glu.sum$group),]

names(f.glu.sum)[c(1,6)] <- c("lab", "group")

f.glu.plot <- ggplot(data = f.glucose.dt, aes(x = age, y = emmean, fill = group))+
  geom_col(fill = rep(c("purple","chartreuse4",
                        "deepskyblue3","brown1"), each = 3), 
           color = "black") +
  facet_wrap(~group, labeller = as_labeller(glab),
             nrow = 1) +
  scale_x_discrete(labels = rep(c("6", "12", "18"), times = 3)) +
  theme_pubr() +
  theme(legend.position = "none",
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "black", linewidth = 1.5)) +
  geom_errorbar(data = f.glu.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = fg.lab,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Females Glucose", y ="Concentration (mg/mL)",
       x = "Age (months)")
f.glu.plot

#Lactate
ggqqplot(hff.met$lactate)

hf.lactate.f.lm <- lm(lactate ~ trt * protein * age, data = hff.met)
anova(hf.lactate.f.lm)
summary(hf.lactate.f.lm)

f.lactate.em <- emmeans(hf.lactate.f.lm, specs = c("trt", "protein", "age"))
f.lactate.em

f.lactate.dt <- data.table(summary(f.lactate.em))
f.lactate.dt

# Contrast "pairwise comparisions"
f.lactate.pairs <- contrast(f.lactate.em,
                            method = "revpairwise",
                            #simple = "each",
                            combine = T,
                            adjust = "sidak") %>%
  summary(infer = T, )

f.lactate.pairs.dt <- data.table(f.lactate.pairs)

# Coordinates for significance indicators

g <- gsub("no Beef", "group 4", as.character(comp))
g <- gsub("yes Beef", "group 3", as.character(g))
g <- gsub("no C", "group 6", as.character(g))
g <- gsub("yes C", "group 5", as.character(g))

g <- as.factor(g)
age <- gsub("group [3-6] ", "", as.character(g))
age1 <- gsub(" - H[1-3]", "", as.character(age))
age2 <- gsub("H[1-3] - ", "", as.character(age))
g1 <- gsub(" H[1-3] - group [3-6] H[1-3]", "", as.character(g))
g2 <- gsub("group [3-6] H[1-3] - ", "", as.character(g))
g2 <- gsub(" H[1-3]", "", as.character(g2))

f.lactate.pairs.dt$group1 <- g1
f.lactate.pairs.dt$group2 <- g2
f.lactate.pairs.dt$age1 <- age1
f.lactate.pairs.dt$age2 <- age2

# Nice p-values for graphs
f.lactate.pairs.dt[, p_rounded := p_round(p.value,
                                          digits = 2)]
f.lactate.pairs.dt[, p_pretty := p_format(p_rounded,
                                          digits = 2,
                                          accuracy = 1e-04,
                                          add.p = TRUE)]
f.lactate.pairs.dt

# Graphing labels ()
f.lactate.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5")

f.lactate.pairs.sig <- f.lactate.pairs.dt[which(f.lactate.pairs.dt$p.value <= 0.05),]

f.lactate.pairs.sig <- f.lactate.pairs.sig[,-c(5:10)]


### Customize Here ###
f.lac <- paste("group", rep(c(3:6), each = 3))
f.lac <- data.frame(f.lac)
f.lac$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
f.lac$lab <- c("A", "A", "A",
               "A", "A", "A",
               "A", "B", "AB",
               "A", "A", "A")
names(f.lac)[1] <- "group1"
f.lac$group2 <- NA
f.lac$group <- f.lac$group1


f.lac.sum <- hff.met%>% group_by(group, age) %>%
  summarise(emmean = mean(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))

f.lac.sum <- rbind(g3f.h1.sum, f.lac.sum)

f.lac.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
f.lac.sum$age <- as.factor(f.lac.sum$age)
f.lac.sum <- f.lac.sum[order(f.lac.sum$age),]
f.lactate.dt <- f.lactate.dt[order(f.lactate.dt$age,f.lactate.dt$group),]
f.lac.sum <- f.lac.sum[order(f.lac.sum$group),]

names(f.lac.sum)[c(1,6)] <- c("lab", "group")

f.lac$y <- f.lac.sum$emmean + f.lac.sum$SE + .05

f.lac

# Customize for each set of comp

#Plot

f.lac.plot <- ggplot(data = f.lactate.dt, aes(x = age, y = emmean, fill = group))+
  geom_col(fill = rep(c("purple","chartreuse4",
                        "deepskyblue3","brown1"), each = 3), 
           color = "black") +
  facet_wrap(~group, labeller = as_labeller(glab),
             nrow = 1) +
  scale_x_discrete(labels = rep(c("6", "12", "18"), times = 3)) +
  theme_pubr() +
  theme(legend.position = "none",
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "black", linewidth = 1.5)) +
  geom_errorbar(data = f.lac.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = f.lac,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Females Lactate", y ="Concentration (mg/mL)",
       x = "Age (months)")
f.lac.plot

## Alanine ##
ggqqplot(hff.met$alanine)

hf.ala.f.lm <- lm(alanine ~ trt * protein * age, data = hff.met)
anova(hf.ala.f.lm)
summary(hf.ala.f.lm)

f.ala.em <- emmeans(hf.ala.f.lm, specs = c("trt", "protein", "age"))
f.ala.em


f.ala.dt <- data.table(summary(f.ala.em))
f.ala.dt

# Contrast "pairwise comparisions"
f.ala.pairs <- contrast(f.ala.em,
                        method = "revpairwise",
                        #simple = "each",
                        combine = T,
                        adjust = "sidak") %>%
  summary(infer = T, )

f.ala.pairs.dt <- data.table(f.ala.pairs)
f.ala.pairs.dt

# Coordinates for significance indicators
comp <- f.ala.pairs.dt$contrast

g <- gsub("no Beef", "group 4", as.character(comp))
g <- gsub("yes Beef", "group 3", as.character(g))
g <- gsub("no C", "group 6", as.character(g))
g <- gsub("yes C", "group 5", as.character(g))

g <- as.factor(g)
age <- gsub("group [3-6] ", "", as.character(g))
age1 <- gsub(" - H[1-3]", "", as.character(age))
age2 <- gsub("H[1-3] - ", "", as.character(age))
g1 <- gsub(" H[1-3] - group [3-6] H[1-3]", "", as.character(g))
g2 <- gsub("group [3-6] H[1-3] - ", "", as.character(g))
g2 <- gsub(" H[1-3]", "", as.character(g2))

f.ala.pairs.dt$group1 <- g1
f.ala.pairs.dt$group2 <- g2
f.ala.pairs.dt$age1 <- age1
f.ala.pairs.dt$age2 <- age2

# Nice p-values for graphs
f.ala.pairs.dt[, p_rounded := p_round(p.value,
                                      digits = 2)]
f.ala.pairs.dt[, p_pretty := p_format(p_rounded,
                                      digits = 2,
                                      accuracy = 1e-04,
                                      add.p = TRUE)]
f.ala.pairs.dt

# Graphing labels ()
f.ala.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5")
f.ala.pairs.sig <- f.ala.pairs.dt[which(f.ala.pairs.dt$p.value <= 0.05),]
f.ala.pairs.sig[which(f.ala.pairs.sig$group1 == f.ala.pairs.sig$group2),]


f.ala.pairs.sig

#Custom Comparisons

f.ala <- paste("group", rep(c(3:6), each = 3))
f.ala <- data.frame(f.ala)
f.ala$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
f.ala$lab <- c("A", "A", "A",
               "A", "A", "A",
               "A", "A", "A",
               "A", "A", "A")
names(f.ala)[1] <- "group1"
f.ala$group2 <- NA
f.ala$group <- f.ala$group1


f.ala.sum <- hff.met%>% group_by(group, age) %>%
  summarise(emmean = mean(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))
f.ala.sum <- rbind(g3f.h1.sum, f.ala.sum)

f.ala.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
f.ala.sum$age <- as.factor(f.ala.sum$age)
f.ala.sum <- f.ala.sum[order(f.ala.sum$age),]
f.ala.dt <- f.ala.dt[order(f.ala.dt$age,f.ala.dt$group),]
f.ala.sum <- f.ala.sum[order(f.ala.sum$group),]

names(f.ala.sum)[c(1,6)] <- c("lab", "group")

f.ala$y <- f.ala.sum$emmean + f.ala.sum$SE + 10

f.ala
#Plot
f.ala.plot <- ggplot(data = f.ala.dt, aes(x = age, y = emmean, fill = group))+
  geom_col(fill = rep(c("purple","chartreuse4",
                        "deepskyblue3","brown1"), each = 3), 
           color = "black") +
  facet_wrap(~group, labeller = as_labeller(glab),
             nrow = 1) +
  scale_x_discrete(labels = rep(c("6", "12", "18"), times = 3)) +
  theme_pubr() +
  theme(legend.position = "none",
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "black", linewidth = 1.5)) +
  geom_errorbar(data = f.ala.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = f.ala,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Females Alanine", y ="Concentration (\U00B5g/mL)",
       x = "Age (months)")
f.ala.plot

## Glutamine ##
ggqqplot(hff.met$glutamine)

hf.gln.f.lm <- lm(glutamine ~ trt * protein * age, data = hff.met)
anova(hf.gln.f.lm)
summary(hf.gln.f.lm)

f.gln.em <- emmeans(hf.gln.f.lm, specs = c("trt", "protein", "age"))
f.gln.em


f.gln.dt <- data.table(summary(f.gln.em))
f.gln.dt

# Contrast "pairwise comparisions"
f.gln.pairs <- contrast(f.gln.em,
                        method = "revpairwise",
                        #simple = "each",
                        combine = T,
                        adjust = "sidak") %>%
  summary(infer = T, )

f.gln.pairs.dt <- data.table(f.gln.pairs)
f.gln.pairs.dt

# Coordinates for significance indicators
comp <- f.gln.pairs.dt$contrast

g <- gsub("no Beef", "group 4", as.character(comp))
g <- gsub("yes Beef", "group 3", as.character(g))
g <- gsub("no C", "group 6", as.character(g))
g <- gsub("yes C", "group 5", as.character(g))

g <- as.factor(g)
age <- gsub("group [3-6] ", "", as.character(g))
age1 <- gsub(" - H[1-3]", "", as.character(age))
age2 <- gsub("H[1-3] - ", "", as.character(age))
g1 <- gsub(" H[1-3] - group [3-6] H[1-3]", "", as.character(g))
g2 <- gsub("group [3-6] H[1-3] - ", "", as.character(g))
g2 <- gsub(" H[1-3]", "", as.character(g2))

f.gln.pairs.dt$group1 <- g1
f.gln.pairs.dt$group2 <- g2
f.gln.pairs.dt$age1 <- age1
f.gln.pairs.dt$age2 <- age2

# Nice p-values for graphs
f.gln.pairs.dt[, p_rounded := p_round(p.value,
                                      digits = 2)]
f.gln.pairs.dt[, p_pretty := p_format(p_rounded,
                                      digits = 2,
                                      accuracy = 1e-04,
                                      add.p = TRUE)]
f.gln.pairs.dt

# Graphing labels ()
f.gln.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5")
f.gln.pairs.sig <- f.gln.pairs.dt[which(f.gln.pairs.dt$p.value <= 0.05),]
f.gln.pairs.sig[which(f.gln.pairs.sig$group1 == f.gln.pairs.sig$group2),]


f.gln.pairs.sig

#Custom Comparisons

f.gln <- paste("group", rep(c(3:6), each = 3))
f.gln <- data.frame(f.gln)
f.gln$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
f.gln$lab <- c("A", "A", "A",
               "A", "A", "A",
               "A", "A", "B",
               "A", "A", "A")
names(f.gln)[1] <- "group1"
f.gln$group2 <- NA
f.gln$group <- f.gln$group1


f.gln.sum <- hff.met%>% group_by(group, age) %>%
  summarise(emmean = mean(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

f.gln.sum <- rbind(g3f.h1.sum, f.gln.sum)

f.gln.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
f.gln.sum$age <- as.factor(f.gln.sum$age)
f.gln.sum <- f.gln.sum[order(f.gln.sum$age),]
f.gln.dt <- f.gln.dt[order(f.gln.dt$age,f.gln.dt$group),]
f.gln.sum <- f.gln.sum[order(f.gln.sum$group),]

names(f.gln.sum)[c(1,6)] <- c("lab", "group")

f.gln$y <- f.gln.sum$emmean + f.gln.sum$SE + 10

f.gln
#Plot
f.gln.plot <- ggplot(data = f.gln.dt, aes(x = age, y = emmean, fill = group))+
  geom_col(fill = rep(c("purple","chartreuse4",
                        "deepskyblue3","brown1"), each = 3), 
           color = "black") +
  facet_wrap(~group, labeller = as_labeller(glab),
             nrow = 1) +
  scale_x_discrete(labels = rep(c("6", "12", "18"), times = 3)) +
  theme_pubr() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "black", linewidth = 1.5)) +
  geom_errorbar(data = f.gln.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = f.gln,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Females Glutamine", y ="Concentration (\U00B5g/mL)",
       x = "Age (months)")
f.gln.plot


# Carnitine #
ggqqplot(hff.met$carnitine)

hf.crn.f.lm <- lm(carnitine ~ trt * protein * age, data = hff.met)
anova(hf.crn.f.lm)
summary(hf.crn.f.lm)

f.crn.em <- emmeans(hf.crn.f.lm, specs = c("trt", "protein", "age"))
f.crn.em


f.crn.dt <- data.table(summary(f.crn.em))
f.crn.dt

# Contrast "pairwise comparisions"
f.crn.pairs <- contrast(f.crn.em,
                        method = "revpairwise",
                        #simple = "each",
                        combine = T,
                        adjust = "sidak") %>%
  summary(infer = T, )

f.crn.pairs.dt <- data.table(f.crn.pairs)
f.crn.pairs.dt

# Coordinates for significance indicators
comp <- f.crn.pairs.dt$contrast

g <- gsub("no Beef", "group 4", as.character(comp))
g <- gsub("yes Beef", "group 3", as.character(g))
g <- gsub("no C", "group 6", as.character(g))
g <- gsub("yes C", "group 5", as.character(g))

g <- as.factor(g)
age <- gsub("group [3-6] ", "", as.character(g))
age1 <- gsub(" - H[1-3]", "", as.character(age))
age2 <- gsub("H[1-3] - ", "", as.character(age))
g1 <- gsub(" H[1-3] - group [3-6] H[1-3]", "", as.character(g))
g2 <- gsub("group [3-6] H[1-3] - ", "", as.character(g))
g2 <- gsub(" H[1-3]", "", as.character(g2))

f.crn.pairs.dt$group1 <- g1
f.crn.pairs.dt$group2 <- g2
f.crn.pairs.dt$age1 <- age1
f.crn.pairs.dt$age2 <- age2

# Nice p-values for graphs
f.crn.pairs.dt[, p_rounded := p_round(p.value,
                                      digits = 2)]
f.crn.pairs.dt[, p_pretty := p_format(p_rounded,
                                      digits = 2,
                                      accuracy = 1e-04,
                                      add.p = TRUE)]
f.crn.pairs.dt

# Graphing labels ()
f.crn.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5")
f.crn.pairs.sig <- f.crn.pairs.dt[which(f.crn.pairs.dt$p.value <= 0.05),]
f.crn.pairs.sig[which(f.crn.pairs.sig$group1 == f.crn.pairs.sig$group2),]

f.crn.pairs.sig

#Custom Comparisons

f.crn <- paste("group", rep(c(3:6), each = 3))
f.crn <- data.frame(f.crn)
f.crn$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
f.crn$lab <- c("A", "A", "A",
               "A", "A", "A",
               "A", "A", "A",
               "A", "A", "A")
names(f.crn)[1] <- "group1"
f.crn$group2 <- NA
f.crn$group <- f.crn$group1


f.crn.sum <- hff.met%>% group_by(group, age) %>%
  summarise(emmean = mean(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

f.crn.sum <- rbind(g3f.h1.sum, f.crn.sum)

f.crn.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
f.crn.sum$age <- as.factor(f.crn.sum$age)
f.crn.sum <- f.crn.sum[order(f.crn.sum$age),]
f.crn.dt <- f.crn.dt[order(f.crn.dt$age,f.crn.dt$group),]
f.crn.sum <- f.crn.sum[order(f.crn.sum$group),]

names(f.crn.sum)[c(1,6)] <- c("lab", "group")

f.crn$y <- f.crn.sum$emmean + f.crn.sum$SE + 2

f.crn
#Plot
f.crn.plot <- ggplot(data = f.crn.dt, aes(x = age, y = emmean, fill = group))+
  geom_col(fill = rep(c("purple","chartreuse4",
                        "deepskyblue3","brown1"), each = 3), 
           color = "black") +
  facet_wrap(~group, labeller = as_labeller(glab),
             nrow = 1) +
  scale_x_discrete(labels = rep(c("6", "12", "18"), times = 3)) +
  theme_pubr() +
  theme(legend.position = "none",
        panel.border = element_rect(fill = "transparent", # Needed to add the border
                                    color = "black", linewidth = 1.5)) +
  geom_errorbar(data = f.crn.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = f.crn,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Females Carnitine", y ="Concentration (\U00B5g/mL)",
       x = "Age (months)")
f.crn.plot

ggarrange(f.glu.plot, f.lac.plot)
ggarrange(f.ala.plot, f.gln.plot)
ggarrange(NULL,f.crn.plot, NULL, ncol = 3, widths = c(.5, 1, .5))

write.csv(anova(hf.glucose.f.lm), "./Data/f.glucose.anova.csv")
write.csv(anova(hf.lactate.f.lm), "./Data/f.lactate.anova.csv")
write.csv(anova(hf.ala.f.lm), "./Data/f.ala.anova.csv")
write.csv(anova(hf.gln.f.lm), "./Data/f.gln.anova.csv")
write.csv(anova(hf.crn.f.lm), "./Data/f.crn.anova.csv")
