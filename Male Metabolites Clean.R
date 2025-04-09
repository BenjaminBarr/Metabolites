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

m.met <- met[which(met$sex == "M"),] # Males
hfm.met <- m.met[which(m.met$group %in% c("G3", "G4", "G5", "G6")),] #HF only

hfm.met
## Males ##
#Glucose First#

ggqqplot(hfm.met$glucose)

hf.glucose.m.lm <- lm(glucose ~ trt * protein * age, data = hfm.met)
anova(hf.glucose.m.lm)
summary(hf.glucose.m.lm)


m.glucose.em <- emmeans(hf.glucose.m.lm, specs = c("trt", "protein", "age"))
m.glucose.em

m.glucose.dt <- data.table(summary(m.glucose.em))
m.glucose.dt


# Contrast "pairwise comparisions"
m.glucose.pairs <- contrast(m.glucose.em,
                            method = "revpairwise",
                            #simple = "each",
                            combine = T,
                            adjust = "sidak") %>%
  summary(infer = T, )

m.glucose.pairs.dt <- data.table(m.glucose.pairs)
m.glucose.pairs.dt[which(m.glucose.pairs.dt$p.value <= 0.05)]

# Coordinates for significance indicators
comp <- m.glucose.pairs.dt$contrast

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

m.glucose.pairs.dt$group1 <- g1
m.glucose.pairs.dt$group2 <- g2
m.glucose.pairs.dt$age1 <- age1
m.glucose.pairs.dt$age2 <- age2

# Nice p-values for graphs
m.glucose.pairs.dt[, p_rounded := p_round(p.value,
                                          digits = 2)]
m.glucose.pairs.dt[, p_pretty := p_format(p_rounded,
                                          digits = 2,
                                          accuracy = 1e-04,
                                          add.p = TRUE)]
m.glucose.pairs.dt

# Graphing labels ()
m.glucose.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5")


m.glucose.pairs.sig <- m.glucose.pairs.dt[which(m.glucose.pairs.dt$p.value <= 0.05),]
names(m.glucose.pairs.sig)

# Identify Significance #
m.glucose.pairs.sig

# Within Group Indications #
m.glucose.pairs.sig[which(m.glucose.pairs.sig$group1 == m.glucose.pairs.sig$group2),]

m.glucose.pairs.sig <- m.glucose.pairs.sig[,-c(1,5:7)]
m.glucose.pairs.sig

### Customize Here ###
mg.lab <- paste("group", rep(c(3:6), each = 3))
mg.lab <- data.frame(mg.lab)
mg.lab$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
mg.lab$lab <- c("A", "B", "B",
                "A", "A", "A",
                "A", "AB", "B",
                "A", "A", "A")
names(mg.lab)[1] <- "group1"
mg.lab$group2 <- NA
mg.lab$group <- mg.lab$group1


# Customize for each set of comp
alab <- c(`H1` = "6 Months", `H2` = "12 Months", 
          `H3` = "18 Months")

glab <- c(`group 3` = "HFBN",
          `group 4` = "HFB",
          `group 5` = "HFCN",
          `group 6` = "HFC")

names(mg.lab)[1:2] <- c("group1", "age")

#Plot
g.glu.sum <- hfm.met%>% group_by(group, age) %>%
  summarise(emmean = mean(glucose),
            SE = sd(glucose)/sqrt(length(glucose)),
            n = length(glucose))

mg.lab$y <- g.glu.sum$emmean + g.glu.sum$SE + 1

g.glu.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
g.glu.sum$age <- as.factor(g.glu.sum$age)
g.glu.sum <- g.glu.sum[order(g.glu.sum$age),]
m.glucose.dt <- m.glucose.dt[order(m.glucose.dt$age,m.glucose.dt$group),]
g.glu.sum <- g.glu.sum[order(g.glu.sum$group),]

names(g.glu.sum)[c(1,6)] <- c("lab", "group")

m.glu.plot <- ggplot(data = m.glucose.dt, aes(x = age, y = emmean, fill = group))+
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
  geom_errorbar(data = g.glu.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = mg.lab,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Males Glucose", y ="Concentration (mg/mL)",
       x = "Age (months)")
m.glu.plot

#Lactate
ggqqplot(hfm.met$lactate)

hf.lactate.m.lm <- lm(lactate ~ trt * protein * age, data = hfm.met)
anova(hf.lactate.m.lm)
summary(hf.lactate.m.lm)

m.lactate.em <- emmeans(hf.lactate.m.lm, specs = c("trt", "protein", "age"))
m.lactate.em

m.lactate.dt <- data.table(summary(m.lactate.em))
m.lactate.dt

# Contrast "pairwise comparisions"
m.lactate.pairs <- contrast(m.lactate.em,
                            method = "revpairwise",
                            #simple = "each",
                            combine = T,
                            adjust = "sidak") %>%
  summary(infer = T, )

m.lactate.pairs.dt <- data.table(m.lactate.pairs)

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

m.lactate.pairs.dt$group1 <- g1
m.lactate.pairs.dt$group2 <- g2
m.lactate.pairs.dt$age1 <- age1
m.lactate.pairs.dt$age2 <- age2

# Nice p-values for graphs
m.lactate.pairs.dt[, p_rounded := p_round(p.value,
                                          digits = 2)]
m.lactate.pairs.dt[, p_pretty := p_format(p_rounded,
                                          digits = 2,
                                          accuracy = 1e-04,
                                          add.p = TRUE)]
m.lactate.pairs.dt

# Graphing labels ()
m.lactate.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5",
                        "group 4", "group 3", "group 6", "group 5")

m.lactate.pairs.sig <- m.lactate.pairs.dt[which(m.lactate.pairs.dt$p.value <= 0.05),]

m.lactate.pairs.sig <- m.lactate.pairs.sig[,-c(5:10)]

### Customize Here ###
m.lac <- paste("group", rep(c(3:6), each = 3))
m.lac <- data.frame(m.lac)
m.lac$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
m.lac$lab <- c("A", "AB", "B",
               "A", "A", "A",
               "A", "A", "A",
               "A", "A", "A")
names(m.lac)[1] <- "group1"
m.lac$group2 <- NA
m.lac$group <- m.lac$group1


m.lac.sum <- hfm.met%>% group_by(group, age) %>%
  summarise(emmean = mean(lactate),
            SE = sd(lactate)/sqrt(length(lactate)),
            n = length(lactate))

m.lac.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
m.lac.sum$age <- as.factor(m.lac.sum$age)
m.lac.sum <- m.lac.sum[order(m.lac.sum$age),]
m.lactate.dt <- m.lactate.dt[order(m.lactate.dt$age,m.lactate.dt$group),]
m.lac.sum <- m.lac.sum[order(m.lac.sum$group),]

names(m.lac.sum)[c(1,6)] <- c("lab", "group")

m.lac$y <- m.lac.sum$emmean + m.lac.sum$SE + .05

m.lac

# Customize for each set of comp

#Plot

m.lac.plot <- ggplot(data = m.lactate.dt, aes(x = age, y = emmean, fill = group))+
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
  geom_errorbar(data = m.lac.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m.lac,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Males Lactate", y ="Concentration (mg/mL)",
       x = "Age (months)")
m.lac.plot

## Alanine ##
ggqqplot(hfm.met$alanine)

hf.ala.m.lm <- lm(alanine ~ trt * protein * age, data = hfm.met)
anova(hf.ala.m.lm)
summary(hf.ala.m.lm)

m.ala.em <- emmeans(hf.ala.m.lm, specs = c("trt", "protein", "age"))
m.ala.em


m.ala.dt <- data.table(summary(m.ala.em))
m.ala.dt

# Contrast "pairwise comparisions"
m.ala.pairs <- contrast(m.ala.em,
                        method = "revpairwise",
                        #simple = "each",
                        combine = T,
                        adjust = "sidak") %>%
  summary(infer = T, )

m.ala.pairs.dt <- data.table(m.ala.pairs)
m.ala.pairs.dt

# Coordinates for significance indicators
comp <- m.ala.pairs.dt$contrast

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

m.ala.pairs.dt$group1 <- g1
m.ala.pairs.dt$group2 <- g2
m.ala.pairs.dt$age1 <- age1
m.ala.pairs.dt$age2 <- age2

# Nice p-values for graphs
m.ala.pairs.dt[, p_rounded := p_round(p.value,
                                      digits = 2)]
m.ala.pairs.dt[, p_pretty := p_format(p_rounded,
                                      digits = 2,
                                      accuracy = 1e-04,
                                      add.p = TRUE)]
m.ala.pairs.dt

# Graphing labels ()
m.ala.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5")
m.ala.pairs.sig <- m.ala.pairs.dt[which(m.ala.pairs.dt$p.value <= 0.05),]
m.ala.pairs.sig[which(m.ala.pairs.sig$group1 == m.ala.pairs.sig$group2),]


m.ala.pairs.sig

#Custom Comparisons

m.ala <- paste("group", rep(c(3:6), each = 3))
m.ala <- data.frame(m.ala)
m.ala$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
m.ala$lab <- c("A", "AB", "B",
               "AB", "A", "B",
               "A", "A", "A",
               "A", "A", "A")
names(m.ala)[1] <- "group1"
m.ala$group2 <- NA
m.ala$group <- m.ala$group1


m.ala.sum <- hfm.met%>% group_by(group, age) %>%
  summarise(emmean = mean(alanine),
            SE = sd(alanine)/sqrt(length(alanine)),
            n = length(alanine))

m.ala.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
m.ala.sum$age <- as.factor(m.ala.sum$age)
m.ala.sum <- m.ala.sum[order(m.ala.sum$age),]
m.ala.dt <- m.ala.dt[order(m.ala.dt$age,m.ala.dt$group),]
m.ala.sum <- m.ala.sum[order(m.ala.sum$group),]

names(m.ala.sum)[c(1,6)] <- c("lab", "group")

m.ala$y <- m.ala.sum$emmean + m.ala.sum$SE + 10

m.ala
#Plot
m.ala.plot <- ggplot(data = m.ala.dt, aes(x = age, y = emmean, fill = group))+
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
  geom_errorbar(data = m.ala.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m.ala,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Males Alanine", y ="Concentration (\U00B5g/mL)",
       x = "Age (months)")
m.ala.plot

## Glutamine ##
ggqqplot(hfm.met$glutamine)

hf.gln.m.lm <- lm(glutamine ~ trt * protein * age, data = hfm.met)
anova(hf.gln.m.lm)
summary(hf.gln.m.lm)

m.gln.em <- emmeans(hf.gln.m.lm, specs = c("trt", "protein", "age"))
m.gln.em


m.gln.dt <- data.table(summary(m.gln.em))
m.gln.dt

# Contrast "pairwise comparisions"
m.gln.pairs <- contrast(m.gln.em,
                        method = "revpairwise",
                        #simple = "each",
                        combine = T,
                        adjust = "sidak") %>%
  summary(infer = T, )

m.gln.pairs.dt <- data.table(m.gln.pairs)
m.gln.pairs.dt

# Coordinates for significance indicators
comp <- m.gln.pairs.dt$contrast

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

m.gln.pairs.dt$group1 <- g1
m.gln.pairs.dt$group2 <- g2
m.gln.pairs.dt$age1 <- age1
m.gln.pairs.dt$age2 <- age2

# Nice p-values for graphs
m.gln.pairs.dt[, p_rounded := p_round(p.value,
                                      digits = 2)]
m.gln.pairs.dt[, p_pretty := p_format(p_rounded,
                                      digits = 2,
                                      accuracy = 1e-04,
                                      add.p = TRUE)]
m.gln.pairs.dt

# Graphing labels ()
m.gln.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5")
m.gln.pairs.sig <- m.gln.pairs.dt[which(m.gln.pairs.dt$p.value <= 0.05),]
m.gln.pairs.sig[which(m.gln.pairs.sig$group1 == m.gln.pairs.sig$group2),]


m.gln.pairs.sig

#Custom Comparisons

m.gln <- paste("group", rep(c(3:6), each = 3))
m.gln <- data.frame(m.gln)
m.gln$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
m.gln$lab <- c("A", "A", "A",
               "A", "A", "A",
               "A", "AB", "B",
               "A", "A", "A")
names(m.gln)[1] <- "group1"
m.gln$group2 <- NA
m.gln$group <- m.gln$group1


m.gln.sum <- hfm.met%>% group_by(group, age) %>%
  summarise(emmean = mean(glutamine),
            SE = sd(glutamine)/sqrt(length(glutamine)),
            n = length(glutamine))

m.gln.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
m.gln.sum$age <- as.factor(m.gln.sum$age)
m.gln.sum <- m.gln.sum[order(m.gln.sum$age),]
m.gln.dt <- m.gln.dt[order(m.gln.dt$age,m.gln.dt$group),]
m.gln.sum <- m.gln.sum[order(m.gln.sum$group),]

names(m.gln.sum)[c(1,6)] <- c("lab", "group")

m.gln$y <- m.gln.sum$emmean + m.gln.sum$SE + 10

m.gln
#Plot
m.gln.plot <- ggplot(data = m.gln.dt, aes(x = age, y = emmean, fill = group))+
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
  geom_errorbar(data = m.gln.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m.gln,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Males Glutamine", y ="Concentration (\U00B5g/mL)",
       x = "Age (months)")
m.gln.plot


# Carnitine #
ggqqplot(hfm.met$carnitine)

hf.crn.m.lm <- lm(carnitine ~ trt * protein * age, data = hfm.met)
anova(hf.crn.m.lm)
summary(hf.crn.m.lm)

m.crn.em <- emmeans(hf.crn.m.lm, specs = c("trt", "protein", "age"))
m.crn.em


m.crn.dt <- data.table(summary(m.crn.em))
m.crn.dt

# Contrast "pairwise comparisions"
m.crn.pairs <- contrast(m.crn.em,
                        method = "revpairwise",
                        #simple = "each",
                        combine = T,
                        adjust = "sidak") %>%
  summary(infer = T, )

m.crn.pairs.dt <- data.table(m.crn.pairs)
m.crn.pairs.dt

# Coordinates for significance indicators
comp <- m.crn.pairs.dt$contrast

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

m.crn.pairs.dt$group1 <- g1
m.crn.pairs.dt$group2 <- g2
m.crn.pairs.dt$age1 <- age1
m.crn.pairs.dt$age2 <- age2

# Nice p-values for graphs
m.crn.pairs.dt[, p_rounded := p_round(p.value,
                                      digits = 2)]
m.crn.pairs.dt[, p_pretty := p_format(p_rounded,
                                      digits = 2,
                                      accuracy = 1e-04,
                                      add.p = TRUE)]
m.crn.pairs.dt

# Graphing labels ()
m.crn.dt$group <- c("group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5",
                    "group 4", "group 3", "group 6", "group 5")
m.crn.pairs.sig <- m.crn.pairs.dt[which(m.crn.pairs.dt$p.value <= 0.05),]
m.crn.pairs.sig[which(m.crn.pairs.sig$group1 == m.crn.pairs.sig$group2),]


m.crn.pairs.sig

#Custom Comparisons

m.crn <- paste("group", rep(c(3:6), each = 3))
m.crn <- data.frame(m.crn)
m.crn$age <- paste(rep(c("H1", "H2", "H3"), times = 4))
m.crn$lab <- c("A", "A", "A",
               "A", "B", "A",
               "A", "A", "A",
               "A", "A", "A")
names(m.crn)[1] <- "group1"
m.crn$group2 <- NA
m.crn$group <- m.crn$group1


m.crn.sum <- hfm.met%>% group_by(group, age) %>%
  summarise(emmean = mean(carnitine),
            SE = sd(carnitine)/sqrt(length(carnitine)),
            n = length(carnitine))

m.crn.sum$lab <- rep(c("group 3", "group 4", "group 5", "group 6"), each = 3)
m.crn.sum$age <- as.factor(m.crn.sum$age)
m.crn.sum <- m.crn.sum[order(m.crn.sum$age),]
m.crn.dt <- m.crn.dt[order(m.crn.dt$age,m.crn.dt$group),]
m.crn.sum <- m.crn.sum[order(m.crn.sum$group),]

names(m.crn.sum)[c(1,6)] <- c("lab", "group")

m.crn$y <- m.crn.sum$emmean + m.crn.sum$SE + 2

m.crn
#Plot
m.crn.plot <- ggplot(data = m.crn.dt, aes(x = age, y = emmean, fill = group))+
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
  geom_errorbar(data = m.crn.sum, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = age), width = 0.1, color = "black") +
  stat_pvalue_manual(data = m.crn,
                     label = "lab",
                     x = "age",
                     y.position = "y") +
  labs(title = "HF Males Carnitine", y ="Concentration (\U00B5g/mL)",
       x = "Age (months)")
m.crn.plot

ggarrange(m.glu.plot, m.lac.plot)
ggarrange(m.ala.plot, m.gln.plot)
ggarrange(NULL,m.crn.plot, NULL, ncol = 3, widths = c(.5, 1, .5))

write.csv(anova(hf.glucose.m.lm), "./Data/m.glucose.anova.csv")
write.csv(anova(hf.lactate.m.lm), "./Data/m.lactate.anova.csv")
write.csv(anova(hf.ala.m.lm), "./Data/m.ala.anova.csv")
write.csv(anova(hf.gln.m.lm), "./Data/m.gln.anova.csv")
write.csv(anova(hf.crn.m.lm), "./Data/m.crn.anova.csv")
