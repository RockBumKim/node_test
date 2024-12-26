library(haven)
library(epiDisplay)
library(tidyverse)
library(gtsummary)
library(knitr)
library(kableExtra)
library(rmarkdown)
library(broom)
library(xtable)
library(PSweight)
library(cobalt)
library(WeightIt)
library(sjstats)
library(questionr)
library(weights)
library(nnet)
library(mlogit)


aa <- read.csv("D:/csv/jsh.csv")

table(aa$Sex)
table(aa$group4)

tabpct(aa$Sex, aa$group4)
# Total complication 1151 pts
aa$total <- if_else(aa$Complication1 == 0, 0, 1)


# severe complication 362 pts
aa$severe <- if_else(aa$CD31 == 31, 1, 0)

table(aa$severe)
aa$severe <- factor(aa$severe)
# Age
# 479 missing
cc <- aa %>% filter (!is.na(Age))
cc <- cc %>% filter (!is.na(Approach))
cc <- cc %>% filter(!is.na(Resectionextent))
cc<- cc %>% filter(!is.na(Resectionextent))

# complete cases: cc

cc$group4 <- if_else(cc$group4 ==1, 4,
                     if_else(cc$group4 ==2, 3,
                             if_else(cc$group4 == 3, 2, 1)) )
cc$group4 <- factor(cc$group4)
tab1(cc$group4)

cc <- cc %>% mutate(group3 = if_else(group4 == 3, 2, 
                                     if_else(group4 == 4, 3, 1)))
tab1(cc$group3)

# group3 : 1 = more than 200 beds,  2 = 100-199 beds,  3 = less than 100 beds



str(cc$Age)
cc$Age <- as.numeric(cc$Age)


# Sex
cc$Sex <- factor(cc$Sex, levels = c(1,2), labels = c("Male", "Female"))


# BMI
str(cc$BMI)
cc$BMI <- as.numeric(cc$BMI)


# Tumor size
str(cc$Tumorsizecm)
cc$Tumorsizecm <- as.numeric(cc$Tumorsizecm)


# TNM stage
str(cc$stage4)
cc$stage4 <- factor(cc$stage4, levels = c(1,2,3), labels = c("I", "II", "III"))


# histologic type
# 3: moderately differentiated tubular adenocarcinoma
# 4: poorly differentiated adenocarcinoma
# 5: poorly cohesive carcinoma 
str(cc$Histologictype)
table(cc$Histologictype)

cc <- cc %>% mutate(Histologictype2 = if_else(is.na(Histologictype), 4, 
                                                    if_else(Histologictype == 3, 1, 
                                                            if_else(Histologictype == 4, 2,
                                                                    if_else(Histologictype ==5, 3, 4)))))

table(cc$Histologictype2)
cc$Histologictype2 <- factor(cc$Histologictype2, levels = c(1:4), labels = c("moderately differentiated ADC","poorly differentiated ADC", "poorly cohesive carcinoma", "others or unknown"))


# Pre OP chemotherapy
tab1(cc$Neoadjuvantchemotherapy)
cc$Neoadjuvantchemotherapy <- if_else(cc$Neoadjuvantchemotherapy  %in% c(2,3), 1,0)


# comorbidity
str(cc$Comorbidity2)
table(cc$Comorbidity2)

cc$Comorbidity2 <- if_else(is.na(cc$Comorbidity2), 2, cc$Comorbidity2)
cc$Comorbidity2 <- factor(cc$Comorbidity2, levels = c(0,1,2), labels = c("No", "Yes", "unknown"))


# ASA score
str(cc$ASAscore)
table(cc$ASAscore)
cc$ASAscore <- as.numeric(cc$ASAscore)

cc$ASAscore <- if_else(is.na(cc$ASAscore), 6, cc$ASAscore)
cc$ASAscore <- factor(cc$ASAscore, levels = c(1:6), labels = c("1","2","3","4","5","unknown"))


# ECOG
table(cc$ECOGscore)

cc$ECOGscore <- as.numeric(cc$ECOGscore)
cc$ECOGscore <- if_else(is.na(cc$ECOGscore) , 5, cc$ECOGscore)
cc$ECOGscore <- factor(cc$ECOGscore, levels = c(0:5), labels = c("0","1","2","3","4","unknown"))


# Hospital stay
str(cc$hostpialstayday)
cc$hostpialstayday <- as.numeric(cc$hostpialstayday)


# Approach methods.
# 1. totally laparoscopic , 2: laparoscopy-assisted, 3: open gastrectomy , 4: robotic gastrectomy
str(cc$Approach)
tab1(cc$Approach)

cc$Approach <- factor(cc$Approach, levels = c(1:4), labels = c("totally laparoscopic", "laparoscopy-assisted", "open gastrectomy","robotic gastrectomy"))


str(cc$group3)
str(cc$group4)

table(cc$group3)
table(cc$group4)

cc$group3 <- factor(cc$group3)
cc$group4 <- factor(cc$group4)
#============================================================
#PSW : 3 group
#============================================================
ps_model1 <- group3 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent
ps_model1_w <- SumStat(ps.formula = ps_model1,  data = cc, weight = c("overlap"))
ps_model1_w$ess
summary(ps_model1_w, weighted.var = TRUE, metric = "ASD")
plot(ps_model1_w, type = "balance")

library(svglite)
svglite("D:/csv/love_plot_3group.svg", width = 18, height = 9)
plot(ps_model1_w, type = "balance")
dev.off()


tiff("D:/csv/love_plot_3group.tiff", units="in", width=18, height=9, res=600, compression = 'lzw')
plot(ps_model1_w, type = "balance")
dev.off()

# Saving weights
ps_weights <- ps_model1_w$ps.weights

ps_weights <- ps_weights %>%  select(-group3)

cc <- cc %>%  cbind(ps_weights)
names(cc)


print(tab1m, showAllLevels = FALSE, smd = FALSE, test = TRUE) 

#=======================
# comparison of 3 groups
# ======================
# Hospital day
model.2way <- lm(hostpialstayday ~ group3, cc, weights = overlap)

summary(model.2way)

cc1 <- cc %>% filter (group3 ==1)
cc1 <- cc1 %>% filter(!is.na(hostpialstayday))

cc2 <- cc %>% filter (group3 ==2)
cc2 <- cc2 %>% filter(!is.na(hostpialstayday))

cc3 <- cc %>% filter (group3 ==3)
cc3 <- cc3 %>% filter(!is.na(hostpialstayday))

library(Hmisc)
c3 <- weighted.mean(cc3$hostpialstayday, cc3$overlap)
mean3 <- weighted.mean(cc3$hostpialstayday, cc3$overlap)
s3 <- sqrt(sum(cc3$overlap * (cc3$hostpialstayday - mean3)^2))

c2 <- weighted.mean(cc2$hostpialstayday, cc2$overlap)
mean2 <- weighted.mean(cc2$hostpialstayday, cc2$overlap)
s2 <- sqrt(sum(cc2$overlap * (cc2$hostpialstayday - mean2)^2))

c1 <- weighted.mean(cc1$hostpialstayday, cc1$overlap)
mean1 <- weighted.mean(cc1$hostpialstayday, cc1$overlap)
s1 <- sqrt(sum(cc1$overlap * (cc1$hostpialstayday - mean1)^2))

result <- c(c3, s3, c2 ,s2, c1, s1)
result

# Hospital day p-value of 3 groups 
summary(aov(hostpialstayday ~ group3, data = cc, weight = overlap))

# Tumorsizecm
model.2way <- lm(Tumorsizecm ~ group3, cc, weights = overlap)

summary(model.2way)

cc1 <- cc %>% filter (group3 ==1)
cc1 <- cc1 %>% filter(!is.na(Tumorsizecm))

cc2 <- cc %>% filter (group3 ==2)
cc2 <- cc2 %>% filter(!is.na(Tumorsizecm))

cc3 <- cc %>% filter (group3 ==3)
cc3 <- cc3 %>% filter(!is.na(Tumorsizecm))

c3 <- weighted.mean(cc3$Tumorsizecm, cc3$overlap)
mean3 <- weighted.mean(cc3$Tumorsizecm, cc3$overlap)
s3 <- sqrt(sum(cc3$overlap * (cc3$Tumorsizecm - mean3)^2))

c2 <- weighted.mean(cc2$Tumorsizecm, cc2$overlap)
mean2 <- weighted.mean(cc2$Tumorsizecm, cc2$overlap)
s2 <- sqrt(sum(cc2$overlap * (cc2$Tumorsizecm - mean2)^2))

c1 <- weighted.mean(cc1$Tumorsizecm, cc1$overlap)
mean1 <- weighted.mean(cc1$Tumorsizecm, cc1$overlap)
s1 <- sqrt(sum(cc1$overlap * (cc1$Tumorsizecm - mean1)^2))

result <- c(c3, s3, c2 ,s2, c1, s1)
result

#Tumorsizecm p-value of 3 groups
summary(aov(Tumorsizecm ~ group3, data = cc, weight = overlap))


# No.ofharvestedLN
model.2way <- lm(No.ofharvestedLN ~ group3, cc, weights = overlap)

summary(model.2way)

cc1 <- cc %>% filter (group3 ==1)
cc1 <- cc1 %>% filter(!is.na(No.ofharvestedLN))

cc2 <- cc %>% filter (group3 ==2)
cc2 <- cc2 %>% filter(!is.na(No.ofharvestedLN))

cc3 <- cc %>% filter (group3 ==3)
cc3 <- cc3 %>% filter(!is.na(No.ofharvestedLN))

c3 <- weighted.mean(cc3$No.ofharvestedLN, cc3$overlap)
mean3 <- weighted.mean(cc3$No.ofharvestedLN, cc3$overlap)
s3 <- sqrt(sum(cc3$overlap * (cc3$No.ofharvestedLN - mean3)^2))

c2 <- weighted.mean(cc2$No.ofharvestedLN, cc2$overlap)
mean2 <- weighted.mean(cc2$No.ofharvestedLN, cc2$overlap)
s2 <- sqrt(sum(cc2$overlap * (cc2$No.ofharvestedLN - mean2)^2))

c1 <- weighted.mean(cc1$No.ofharvestedLN, cc1$overlap)
mean1 <- weighted.mean(cc1$No.ofharvestedLN, cc1$overlap)
s1 <- sqrt(sum(cc1$overlap * (cc1$No.ofharvestedLN - mean1)^2))

result <- c(c3, s3, c2 ,s2, c1, s1)
result

#No.ofharvestedLN p-value of 3 groups
summary(aov(No.ofharvestedLN ~ group3, data = cc, weight = overlap))














# Weighted t test
cc$weights <- cc$overlap

weighted_ttest(hostpialstayday ~ group3 + weights, cc)

cc1 <- cc %>% filter(group2 == "Distal")
cc2 <- cc %>% filter(group2 == "PPG")
weighted_sd(cc1$Age, cc1$overlap)
weighted_sd(cc2$No.ofharvestedLN, cc2$overlap)

# weighted chi-square test

weighted_chisqtest(Sex ~ group2 + weights, aa)
wtd.table(aa$Sex, aa$group2, weights=aa$w2)

cc_severe <- cc %>% filter (Approach != "open gastrectomy")
cc_open <- cc %>% filter (Approach == "open gastrectomy")



wtd.chi.sq(aa$Sex, aa$group2,weight= aa$overlap)

xtabs(overlap ~  Sex+ group2, data=cc) 

cc$o2 <- cc$overlap*100

table(aa$total)















library(weights)
wtd.chi.sq(bb$Sex, bb$group4)
wtd.chi.sq(bb$Resectionextent , bb$group4, weight=bb$overlap)
weighted.table(bb$Resectionextent , bb$group4, stat = "cprop",weights =bb$overlap)

library(survey)
bbs=svydesign(ids = ~ bb$No,data =bb,weights =bb$overlap)
svychisq(~group4+Sex,bbs)


library(descriptio)
weighted.table(cc$Sex, cc$group3, stat = "cprop",weights =cc$overlap)

install.packages("matrixStats")                    # Install matrixStats package
library("matrixStats") 

weightedMean(cc$hostpialstayday, cc$overlap)                               # Apply weightedMean function

library(jtools)
wtd.mean(cc1$hostpialstayday, cc1$overlap)



bb$overlap2 <- bb$overlap *100

wpct (bb$group4, bb$overlap)



fit <- glm(severe ~ group4+Age+ Sex+ stage4+ Comorbidity2+ ASAscore+ Approach, data=bb, family = "quasibinomial", weights = overlap)
summary(fit)
exp(fit$coefficients)
exp(confint(fit))







#============================================================
#PSW : 3 group
#============================================================
ps_model1 <- group3 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent
ps_model1_w <- SumStat(ps.formula = ps_model1,  data = bb, weight = c("overlap"))
ps_model1_w$ess
summary(ps_model1_w, weighted.var = TRUE, metric = "ASD")
plot(ps_model1_w, type = "balance")

library(svglite)
svglite("D:/csv/love_plot_3group.svg", width = 18, height = 9)
plot(ps_model1_w, type = "balance")
dev.off()


tiff("D:/csv/love_plot_3group.tiff", units="in", width=18, height=9, res=600, compression = 'lzw')
plot(ps_model1_w, type = "balance")
dev.off()
#============================================================
#PSW : 4 group
#============================================================
table(bb$ECOGscore)

# Make model
# ECOG score : too many missing (unknown = 5168)
ps_model1 <- group4 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent

# PS_weight
ps_model1_w <- SumStat(ps.formula = ps_model1,  data = bb, weight = c("overlap"))

ps_model1_w$ess

summary(ps_model1_w, weighted.var = TRUE, metric = "ASD")


# Plotting PS_weight
plot(ps_model1_w, type = "balance")



# Saving weights
ps_weights <- ps_model1_w$ps.weights

ps_weights <- ps_weights %>%  select(-group4)

bb <- bb %>%  cbind(ps_weights)
names(bb)


model.2way <- lm(Tumorsizecm ~ group4, bb, weights = overlap)

summary(model.2way)

bb1 <- bb %>% filter (group4 ==1)
bb2 <- bb %>% filter (group4 ==2)
bb3 <- bb %>% filter (group4 ==3)
bb4 <- bb %>% filter (group4 ==4)

str(bb$Tumorsizecm)

library(Hmisc)

c4 <- weighted.mean(bb4$No.ofharvestedLN, bb4$overlap)

mean4 <- weighted.mean(bb4$Tumorsizecm, bb4$overlap)

s4 <- sqrt(sum(bb4$overlap * (bb4$Tumorsizecm - mean4)^2))

c3 <- weighted.mean(bb3$Tumorsizecm, bb3$overlap)

mean2 <- weighted.mean(bb3$Tumorsizecm, bb3$overlap)

s3 <- sqrt(sum(bb3$overlap * (bb3$Tumorsizecm - mean2)^2))

c2 <- weighted.mean(bb2$Tumorsizecm, bb2$overlap)

mean2 <- weighted.mean(bb2$Tumorsizecm, bb2$overlap)

s2 <- sqrt(sum(bb2$overlap * (bb2$Tumorsizecm - mean2)^2))

c1 <- weighted.mean(bb1$Tumorsizecm, bb1$overlap)

mean1 <- weighted.mean(bb1$Tumorsizecm, bb1$overlap)

s1<- sqrt(sum(bb4$overlap2 * (bb4$Tumorsizecm - 3.643512)^2))

result <- c(c4, s4, c3, s3, c2 ,s2, c1, s1)
result

library(weights)
wtd.chi.sq(bb$Sex, bb$group4)
wtd.chi.sq(bb$Resectionextent , bb$group4, weight=bb$overlap)
weighted.table(bb$Resectionextent , bb$group4, stat = "cprop",weights =bb$overlap)

library(survey)
bbs=svydesign(ids = ~ bb$No,data =bb,weights =bb$overlap)
svychisq(~group4+Sex,bbs)


library(descriptio)
weighted.table(bb$Sex, bb$group4, stat = "cprop",weights =bb$overlap2)

install.packages("matrixStats")                    # Install matrixStats package
library("matrixStats") 

weightedMean(bb4$BMI, bb4$overlap)                               # Apply weightedMean function

library(jtools)
wtd.mean(bb1$No.ofharvestedLN, bb1$overlap)



bb$overlap2 <- bb$overlap *100

wpct (bb$group4, bb$overlap)



fit <- glm(severe ~ group4+Age+ Sex+ stage4+ Comorbidity2+ ASAscore+ Approach, data=bb, family = "quasibinomial", weights = overlap)
summary(fit)
exp(fit$coefficients)
exp(confint(fit))










bb$sex2 <- as.integer(bb$sex2)
bb$group42 <- as.integer(bb$group42)



fit0 <- glm(CD31 ~ group4, data=aa, family = "binomial")
summary(fit0)
exp(fit0$coefficients)
logistic.display(fit0)

fit1 <- glm(CD31 ~ group4 + Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent, data=aa_ps, family = "binomial", weights = overlap)
summary(fit1)
exp(fit1$coefficients)
logistic.display(fit1)




fit2 <- glm(CD31 ~ group4 , data=aa_ps, weights = overlap, family = "binomial")
summary(fit2)
exp(fit1$coefficients)
logistic.display(fit2)





cc <- read.csv("D:/csv/jsh2.csv")

#============================================================
#PSW: 2group
#============================================================
str(cc$group2)
tab1(cc$group2)

# group2 : 0=Distal, 1=PPG
cc$group2 <- if_else(cc$group2 == 1, 0, 1)
  

#missing Age, Sex, TMN stage, Comorbitiies, ASA score, and Approach method
sum(is.na(cc$Age))

cc <- cc %>% filter (!is.na(cc$Age))

table(cc$group2)
names(cc)
table(cc$stage4)
#PSW
model <-group2 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach

ps_model1_w <- SumStat(ps.formula = model,  data = cc, weight = c("overlap"))

ps_model1_w$ess

summary(ps_model1_w, weighted.var = TRUE, metric = "ASD")
# Plotting PS_weight
plot(ps_model1_w, type = "balance")


library(svglite)
svglite("D:/csv/love_plot2.svg", width = 18, height = 9)
plot(ps_model1_w, type = "balance")
dev.off()






#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#================================================================

bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))

# Love plot
tiff("D:/csv/love_plot4.tiff", width = 1800, height = 1800, res = 300, compression = "lzw")

love.plot(W.out,stars = "raw", threshold = .1, abs=T,sample.names = c("Original", "Weighted"), position = "top")

dev.off()



# save the weight value to dataset
aa$weights <- W.out$weights



#===============
#PSW second
#==============
model2 <- group2 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach
ps_model1_w <- SumStat(ps.formula = model2, data = cc, weight = c("overlap"))

ps_model1_w$ess

summary(ps_model1_w, weighted.var = TRUE, metric = "ASD")


# Plotting PS_weight
plot(ps_model1_w, type = "balance")


# Saving weights
ps_weights <- ps_model1_w$ps.weights

ps_weights <- ps_weights %>%  select(-group2)

aa <- aa %>%  cbind(ps_weights)

names(aa)
#===============
#PSW 3
#==============

model <- glm(group2 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach, data=aa, family = "binomial")

aa$psvalue <- predict(model, type="response")

aa$w2 <- ifelse(aa$group2 == 1, 1/aa$psvalue, 1/(1-aa$psvalue))

fit <- glm(severe ~ group2+ Age+ Sex+ stage4+ Comorbidity2+ ASAscore+ Approach,
           data = aa, family=binomial() , weights = (w2))

#=============
#multivariate
table(aa$severe)
fit <- glm(severe ~ group2+ Age+ Sex+ stage4+ Comorbidity2+ ASAscore+ Approach,
           data = aa, family="quasibinomial" , weights = weights)



fit <- glm(CD31 ~ group4 + Age + Sex + stage4 +  Comorbidity2 + ASAscore + Approach + Resectionextent,
           data = aa_ps, family=quasibinomial , weights = weights)



summary(fit)
exp(fit$coefficients)
exp(confint(fit))

odds.ratio(fit, level=0.95)

fit <- glm(total ~ group2, data=aa, family = "binomial")




fit <- glm(severe ~ group2+ Age+ Sex+ stage4+ Comorbidity2+ ASAscore+ Approach,
           data = aa, family="quasibinomial" , weights = weights)


fit <- glm(severe ~ group4+Age+ Sex+ stage4+ Comorbidity2+ ASAscore+ Approach, data=bb, family = "quasibinomial", weights = overlap)
summary(fit)
exp(fit$coefficients)
exp(confint(fit))


cc$No.ofharvestedLN <- as.numeric(cc$No.ofharvestedLN)

# Weighted t test

weighted_ttest(Age ~ group2 + weights, cc)

cc1 <- cc %>% filter(group2 == "Distal")
cc2 <- cc %>% filter(group2 == "PPG")
weighted_sd(cc1$Age, cc1$overlap)
weighted_sd(cc2$No.ofharvestedLN, cc2$overlap)

# weighted chi-square test

weighted_chisqtest(Sex ~ group2 + weights, aa)
wtd.table(aa$Sex, aa$group2, weights=aa$w2)

cc_severe <- cc %>% filter (Approach != "open gastrectomy")
cc_open <- cc %>% filter (Approach == "open gastrectomy")



wtd.chi.sq(aa$Sex, aa$group2,weight= aa$overlap)

xtabs(overlap ~  Sex+ group2, data=cc) 

cc$o2 <- cc$overlap*100

table(aa$total)













# Number of harvested LN nodes
names(bb)
str(bb$No.ofharvestedLN)
bb$No.ofharvestedLN <- as.numeric(bb$No.ofharvestedLN)
summ(bb$No.ofharvestedLN, by = bb$group2)


# Clavien-Dindo classification
table(bb$ClavienDioclassification)


# Severe morbidity
table(bb$CD31, bb$group2)


# Severe morbidity with minimal invasion
table(bb$Approach)


cc <- bb %>% filter(Approach != "open gastrectomy")

table(cc$CD31, cc$group2)


# Severe morbidity with open

dd <- bb %>% filter(Approach == "open gastrectomy")

table(dd$CD31, dd$group2)


# mortality
table(bb$Postoperativemortality)

bb$death <- if_else(bb$Postoperativemortality == 0, 0, 1)



#=============================
table(bb$Complication1, bb$group2)

table(bb$Approach)

summary(bb)

sum(is.na(bb$Age))

cc <- bb %>% filter (!is.na(Age))
cc <- aa %>%  filter(!is.na(Age) & !is.na(Resectionextent) & !is.na(Approach))
sum(is.na(aa$Resectionextent))



fit <- multinom(group4 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent, data=aa_ps)
summary(fit)
# calculate weight

sum(is.na(cc$Approach))

cc <- cc %>% filter (!is.na(Age))

W.out <- weightit(group2 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach  ,
                  data = cc, method = "glm")
W.out #print the output

summary(W.out)
# Check the balance after weighting

bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))

# Love plot
tiff("D:/csv/love_plot4.tiff", width = 1800, height = 1800, res = 300, compression = "lzw")

love.plot(W.out,stars = "raw", threshold = .1, abs=T,sample.names = c("Original", "Weighted"), position = "top")

dev.off()


postscript("D:/csv/love_plot4.eps", horizontal = FALSE, onefile = FALSE, paper = "special",  height = 8, width = 8)
love.plot(W.out,stars = "raw", threshold = .1, abs=T,sample.names = c("Original", "Weighted"), position = "top")
dev.off()



# save the weight value to dataset
aa_ps$weights <- W.out$weights

summ(aa_ps$weights)


fit <- glm(CD31 ~ group4 + Age + Sex + stage4 +  Comorbidity2 + ASAscore + Approach + Resectionextent,
           data = aa_ps, family=quasibinomial , weights = weights)



summary(fit)
exp(fit$coefficients)
exp(confint(fit))

odds.ratio(fit, level=0.95)



# Weighted t test

weighted_ttest(Age ~ group4 + weights, aa_ps)

cc1 <- cc %>% filter(group2 == "Distal")
cc2 <- cc %>% filter(group2 == "PPG")
weighted_sd(cc1$No.ofharvestedLN, cc1$weights)
weighted_sd(cc2$No.ofharvestedLN, cc2$weights)

# weighted chi-square test

weighted_chisqtest(Sex ~ group4 + weights, aa_ps)
wtd.table(aa_ps$Sex, aa_ps$group4, weights=aa_ps$w2)

cc_severe <- cc %>% filter (Approach != "open gastrectomy")
cc_open <- cc %>% filter (Approach == "open gastrectomy")

library(weights)


model.2way <- lm( No.ofharvestedLN~ group4, aa_ps, weights = w2)

summary(model.2way)


chisq.test(aa_ps$Sex, aa_ps$group4)

wtable(aa_min$Sex, aa_min$group4, w = aa_min$w2)


#group4 chi square test

wtd.chi.sq(aa_ps$c14, aa_ps$group4,weight= aa_ps$w2)



xtabs(w2 ~  c14+ group4, data=aa_ps) 


#minimal invasive group
table(aa$Approach, aa$group4)

aa_min <- aa_ps %>% filter(Approach != "open gastrectomy")
aa_open <- aa_ps %>% filter(Approach == "open gastrectomy")







# Extent of gastric resection
# 1: distal gastrectomy, 2: total, 3: proximal, 4: PPG, 5:wedge resection
str(aa$Resectionextent)
tab1(aa$Resectionextent)



aa$Resectionextent <- factor(aa$Resectionextent, levels = c(1:5), labels = c("Distal", "Total", "Proximal", "PPG", "Wedge"))
# aa$Resectionextent <- if_else(aa$Resectionextent == "4" | aa$Resectionextent =="5"|aa$Resectionextent =="3", 3, aa$Resectionextent )


# aa$Resectionextent <- factor(aa$Resectionextent, levels = c(1:3), labels = c("distal gastrectomy", "total gastrectomy", "others"))


# complication: 
table(aa$CD31)
#table(aa$Complicationswithinpostoperative30days)
#tab1(aa$Complicationswithinpostoperative30days)

#aa$Complicationswithinpostoperative30days <- factor(aa$Complicationswithinpostoperative30days, levels = c(0,1), labels = c("No", "Yes"))
#sum(is.na(aa$Complicationswithinpostoperative30days))

#aa$cx <- if_else(aa$Complicationswithinpostoperative30days == "No", 0 ,1)

aa$CD31 <- if_else(aa$CD31 == 31, 1, 0)


# Clavien-Dindo classification
# tab1(aa$ClavienDioclassification)


# Group of hospital


tab1(aa$group4)
aa$group4 <- if_else(aa$group4 == 4, 0, 
                     if_else(aa$group4 ==3, 1, 
                             if_else(aa$group4 ==2, 2, 3)))
aa$group4 <- factor(aa$group4, levels = c(0:3), labels = c("group D","group C","group B","group A"))

names(aa)
str(aa$group4)
bb <- aa

Age + Sex + stage4 +  Comorbidity2 + ASAscore + Approach + Resectionextent


cc <- aa %>%  filter(!is.na(Age) & !is.na(Resectionextent) & !is.na(Approach))
sum(is.na(aa$Resectionextent))

#=======================
setwd("D:/Rwork/")

render("jsh_MD.Rmd")
#=======================
sum(is.na(aa$CD31))

#=======================
# PS Weight
#=======================
summary(aa)
aa_ps <- aa %>% filter(!is.na(Age) & !is.na(Resectionextent)  )
#  519 patients excluded)

table(aa_ps$Resectionextent)

sum(is.na(aa_ps$group4))
# Make model

ps_model1 <- group4 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent

# PS_weight
ps_model1_w <- SumStat(ps.formula = ps_model1,  data = aa_ps, weight = c("overlap"))

ps_model1_w$ess

summary(ps_model1_w, weighted.var = TRUE, metric = "ASD")


# Plotting PS_weight
plot(ps_model1_w, type = "balance")



# Saving weights
ps_weights <- ps_model1_w$ps.weights

ps_weights <- ps_weights %>%  select(-group4)

aa_ps <- aa_ps %>%  cbind(ps_weights)
names(aa_ps)




fit0 <- glm(CD31 ~ group4, data=aa, family = "binomial")
summary(fit0)
exp(fit0$coefficients)
logistic.display(fit0)

fit1 <- glm(CD31 ~ group4 + Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent, data=aa_ps, family = "binomial", weights = overlap)
summary(fit1)
exp(fit1$coefficients)
logistic.display(fit1)




fit2 <- glm(CD31 ~ group4 , data=aa_ps, weights = overlap, family = "binomial")
summary(fit2)
exp(fit1$coefficients)
logistic.display(fit2)







#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Distal vs PPG
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

bb <- aa %>% filter(Resectionextent %in% c("Distal", "PPG") )

bb$group2 <- if_else(bb$Resectionextent == "Distal", 1, 
                   if_else (bb$Resectionextent == "PPG", 4, NA))

tab1(bb$group2)


tabpct (bb$group2, bb$Sex)

tab1(bb$ClavienDioclassification)

cc <- bb %>% filter (!is.na(Age))


# Number of harvested LN nodes
names(bb)
str(bb$No.ofharvestedLN)
bb$No.ofharvestedLN <- as.numeric(bb$No.ofharvestedLN)
summ(bb$No.ofharvestedLN, by = bb$group2)


# Clavien-Dindo classification
table(bb$ClavienDioclassification)


# Severe morbidity
table(bb$CD31, bb$group2)


# Severe morbidity with minimal invasion
table(bb$Approach)


cc <- bb %>% filter(Approach != "open gastrectomy")

table(cc$CD31, cc$group2)


# Severe morbidity with open

dd <- bb %>% filter(Approach == "open gastrectomy")

table(dd$CD31, dd$group2)


# mortality
table(bb$Postoperativemortality)

bb$death <- if_else(bb$Postoperativemortality == 0, 0, 1)



#=============================
table(bb$Complication1, bb$group2)

table(bb$Approach)

summary(bb)

sum(is.na(bb$Age))

cc <- bb %>% filter (!is.na(Age))
cc <- aa %>%  filter(!is.na(Age) & !is.na(Resectionextent) & !is.na(Approach))
sum(is.na(aa$Resectionextent))



fit <- multinom(group4 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach + Resectionextent, data=aa_ps)
summary(fit)
# calculate weight

sum(is.na(cc$Approach))

cc <- cc %>% filter (!is.na(Age))

W.out <- weightit(group2 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach  ,
                  data = cc, method = "glm")
W.out #print the output

summary(W.out)
# Check the balance after weighting

bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))

# Love plot
tiff("D:/csv/love_plot4.tiff", width = 1800, height = 1800, res = 300, compression = "lzw")

love.plot(W.out,stars = "raw", threshold = .1, abs=T,sample.names = c("Original", "Weighted"), position = "top")

dev.off()


postscript("D:/csv/love_plot4.eps", horizontal = FALSE, onefile = FALSE, paper = "special",  height = 8, width = 8)
love.plot(W.out,stars = "raw", threshold = .1, abs=T,sample.names = c("Original", "Weighted"), position = "top")
dev.off()



# save the weight value to dataset
aa_ps$weights <- W.out$weights

summ(aa_ps$weights)


fit <- glm(CD31 ~ group4 + Age + Sex + stage4 +  Comorbidity2 + ASAscore + Approach + Resectionextent,
           data = aa_ps, family=quasibinomial , weights = weights)



summary(fit)
exp(fit$coefficients)
exp(confint(fit))

odds.ratio(fit, level=0.95)



# Weighted t test

weighted_ttest(Age ~ group4 + weights, aa_ps)

cc1 <- cc %>% filter(group2 == "Distal")
cc2 <- cc %>% filter(group2 == "PPG")
weighted_sd(cc1$No.ofharvestedLN, cc1$weights)
weighted_sd(cc2$No.ofharvestedLN, cc2$weights)

# weighted chi-square test

weighted_chisqtest(Sex ~ group4 + weights, aa_ps)
wtd.table(aa_ps$Sex, aa_ps$group4, weights=aa_ps$w2)

cc_severe <- cc %>% filter (Approach != "open gastrectomy")
cc_open <- cc %>% filter (Approach == "open gastrectomy")

library(weights)


model.2way <- lm( No.ofharvestedLN~ group4, aa_ps, weights = w2)

summary(model.2way)


chisq.test(aa_ps$Sex, aa_ps$group4)

wtable(aa_min$Sex, aa_min$group4, w = aa_min$w2)


#group4 chi square test

wtd.chi.sq(aa_ps$c14, aa_ps$group4,weight= aa_ps$w2)



xtabs(w2 ~  c14+ group4, data=aa_ps) 


#minimal invasive group
table(aa$Approach, aa$group4)

aa_min <- aa_ps %>% filter(Approach != "open gastrectomy")
aa_open <- aa_ps %>% filter(Approach == "open gastrectomy")


# Group 4 after weight (w2), mean sd
library(Hmisc)

am1 <- Hmisc::wtd.mean(aaa$No.ofharvestedLN, weights=aaa$w2)
as1 <- sqrt(wtd.var(aaa$No.ofharvestedLN, weights=aaa$w2))

am2 <- Hmisc::wtd.mean(aab$No.ofharvestedLN, weights=aab$w2)
as2 <- sqrt(wtd.var(aab$No.ofharvestedLN, weights=aab$w2))

am3 <- Hmisc::wtd.mean(aac$No.ofharvestedLN, weights=aac$w2)
as3 <- sqrt(wtd.var(aac$No.ofharvestedLN, weights=aac$w2))

am4 <- Hmisc::wtd.mean(aad$No.ofharvestedLN, weights=aad$w2)
as4 <- sqrt(wtd.var(aad$No.ofharvestedLN, weights=aad$w2))

print(c(am1, as1, am2, as2, am3, as3, am4, as4))

chisq.test(aa$Cx1_1, aa$group4)


names(aa_ps)
table(aa_ps$Cx1_1, aa_ps$group4)
aa_ps$death <- if_else(aa_ps$ClavienDioclassification == "50", 1, 0)


table(aa_ps$CD31)

aa_ps <- aa_ps %>% mutate(comp = if_else(CD31 ==1, Complication1,0))

tabpct(aa_ps$comp, aa_ps$group4)




aaa <- aa_ps %>% filter(group4 == "group A")
aab <- aa_ps %>% filter(group4 == "group B")
aac <- aa_ps %>% filter(group4 == "group C")
aad <- aa_ps %>% filter(group4 == "group D")











names(aa_ps)

table(aa_ps$group4, weight= aa_ps$weights)


count(x = aa_ps, group4, wt=w2)

names(aa_ps)

names(bb)

aa_ps$w2 <- aa_ps$weights/10





aa_ps$c1 <- if_else(aa_ps$comp == 1, 1, 0)
aa_ps$c2 <- if_else(aa_ps$comp == 2, 1, 0)
aa_ps$c3 <- if_else(aa_ps$comp == 3, 1, 0)
aa_ps$c4 <- if_else(aa_ps$comp == 4, 1, 0)
aa_ps$c5 <- if_else(aa_ps$comp == 5, 1, 0)
aa_ps$c6 <- if_else(aa_ps$comp == 6, 1, 0)
aa_ps$c7 <- if_else(aa_ps$comp == 7, 1, 0)
aa_ps$c8 <- if_else(aa_ps$comp == 8, 1, 0)
aa_ps$c9 <- if_else(aa_ps$comp == 9, 1, 0)
aa_ps$c10 <- if_else(aa_ps$comp == 10, 1, 0)
aa_ps$c11 <- if_else(aa_ps$comp == 11, 1, 0)
aa_ps$c12 <- if_else(aa_ps$comp == 12, 1, 0)
aa_ps$c13 <- if_else(aa_ps$comp == 13, 1, 0)
aa_ps$c14 <- if_else(aa_ps$comp == 14, 1, 0)



# weighted logistic
fit1 <- glm(CD31 ~ group2+ Age+ Sex+ BMI+ Tumorsizecm+ stage4+ Comorbidity2+ ASAscore+ ECOGscore+ hostpialstayday+ Approach+ No.ofharvestedLN,
            data = cc, family="quasibinomial" , weights = weights)



summary(fit1)
exp(fit1$coefficients)





fit2 <- glm(CD31 ~ group2 + Age + Sex + stage4 +  Comorbidity2 + ASAscore + Approach,
            data = cc, family="quasibinomial" , weights = weights)



summary(fit2)
exp(fit2$coefficients)



confint(fit1)














cc$psvalue <- predict(fit, type = "response")
cc$weight.ATE <- ifelse(cc$group2 == 1, 1/cc$psvalue, 1/(1-cc$psvalue))

fit1 <- glm (CD31 ~ group2 + Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach, data=cc, family = "binomial", weights = overlap2)
exp(fit1$coefficients)
summary(fit1)


fit2 <- glm (CD31 ~ group2 + Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach, data=cc, family = "binomial", weights = weight.ATE)
exp(fit2$coefficients)
summary(fit2)
library(cobalt)
library(WeightIt)
w.out1 <- WeightIt::weightit(
  group2 ~ Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach,
  data = cc, estimand = "ATE", method = "ps")

set.cobalt.options(binary = "std")

love.plot(w.out1)

















# Make model
ps_model1 <- group2 ~ Age + Sex + stage4 +  Comorbidity2 + ASAscore + Approach

# PS_weight
ps_model1_w <- SumStat(ps.formula = ps_model1, data = cc, weight = c("overlap"))

ps_model1_w$ess

summary(ps_model1_w, weighted.var = TRUE, metric = "ASD")


# Plotting PS_weight
plot(ps_model1_w, type = "balance")


# Saving weights
ps_weights <- ps_model1_w$ps.weights

ps_weights <- ps_weights %>%  select(-group2)

cc <- cc %>%  cbind(ps_weights)
cc$overlap2 <- cc$overlap * 100

tabpct(bb$ASAscore, bb$group2)


# Weighted t test
library(sjstats)
weighted_ttest(No.ofharvestedLN ~ group2 + overlap2, cc)
weighted_sd(cc1$No.ofharvestedLN, cc1$overlap2)

# weighted chi-square test
library(questionr)
weighted_chisqtest(Comorbidity2 ~ group2 + weight.ATE, cc)
wtd.table(cc$Comorbidity2, cc$group2, weights=cc$weight.ATE)


cc_c2 <- cc %>% filter (Complication1 == 2)
table(bb$Approach)


table(cc$Complication1)

tabpct(bb$Complication1, bb$group2)


cc1 <- cc %>% filter(group2 == 1)
cc2 <- cc %>% filter(group2 == 2)

summ(bb$BMI, by=bb$group2)

plot(ps_model1_w, type = "balance")


# Logistic
fit <- glm (CD31 ~ group2 + Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach, data=cc, family = "binomial")
logistic.display(fit)

fit1 <- glm (CD31 ~ group2 + Age + Sex + stage4 + Comorbidity2 + ASAscore + Approach, data=cc, family = "binomial", weights = overlap2)
exp(fit1$coefficients)
summary(fit1)


bb<- read.csv("D:/csv/jsh3.csv")
names(bb)
tabpct(aa$Complication1, aa$group4)
