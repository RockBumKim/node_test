data <- cbind(c(100, 225), c(58, 45));
summary(epi2x2(data));
install.packages("fmsb")
library(fmsb)


rd <- read.csv("D:/csv/risk_diff.csv")

ich1 <- riskdifference(rd$ich_a[1], rd$ich_p[1], rd$t1_a[1], rd$t1_p[1], CRC=TRUE)
ich2 <- riskdifference(rd$mi_a[2], rd$mi_p[2], rd$t1_a[2], rd$t1_p[2], CRC=TRUE)
print(ich1)
print(ich2)


# Sensertivity analysis : leave-one-out
sens <- read.csv("D:/csv/asa_sens.csv")

names(sens)

library(meta)

cv <-  metabin(CV_A1, Total_A, CV_P1, Total_P,
               data = sens, studlab = paste(Study, Year), sm = "RR")

mace <-  metabin(MACE_A1, Total_A, MACE_P1, Total_P,
                 data = sens, studlab = paste(Study, Year), sm = "RR")
mi <- metabin(MI_A1, Total_A, MI_P1, Total_P,
        data = sens, studlab = paste(Study, Year), sm = "RR")

se <- metabin(SE_A1, Total_A, SE_P1, Total_P,
        data = sens, studlab = paste(Study, Year), sm = "RR")

is <- metabin(IS_A1, Total_A, IS_P1, Total_P,
        data = sens, studlab = paste(Study, Year), sm = "RR")

hs <- metabin(HS_A1, Total_A, HS_P1, Total_P,
        data = sens, studlab = paste(Study, Year), sm = "RR")

ich <- metabin(ICH_A1, Total_A, ICH_P1, Total_P,
        data = sens, studlab = paste(Study, Year), sm = "RR")

gi <- metabin(GI_A1, Total_A, GI_P1, Total_P,
        data = sens, studlab = paste(Study, Year), sm = "RR")

mb <- metabin(MB_A1, Total_A, MB_P1, Total_P,
              data = sens, studlab = paste(Study, Year), sm = "RR")


metainf(mb)
metainf(mb, pooled = "random")



# Trim and fill
m2 <- metabin(MB_A1, Total_A, MB_P1, Total_P,
              data = sens, studlab = paste(Study, Year), sm = "RR")

cvt <- trimfill(cv)
macet <- trimfill(mace)
mit <- trimfill(mi)
set <- trimfill(se)
ist <- trimfill(is)
hst <- trimfill(hs)
icht <- trimfill(ich)
git <- trimfill(gi)
mbt <- trimfill(mb)

tiff("D:/csv/cv_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(cvt, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/mace_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(macet, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/mi_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(mit, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/se_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(set, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/is_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(ist, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/hs_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(hst, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/ich_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(icht, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/gi_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(git, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

tiff("D:/csv/mb_funnel.tiff", width = 6, height = 6, units = 'in', res = 600, compression = 'lzw')
funnel(mbt, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()



# emf files
library(devEMF)
emf("D:/csv/cv_funnel.emf")
funnel(cvt, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/mace_funnel.emf")
funnel(macet, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/mi_funnel.emf")
funnel(mit, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/se_funnel.emf")
funnel(set, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/is_funnel.emf")
funnel(ist, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/hs_funnel.emf")
funnel(hst, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/ich_funnel.emf",)
funnel(icht, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/gi_funnel.emf")
funnel(git, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()

emf("D:/csv/mb_funnel.emf")
funnel(mbt, cex = 1.5, studlab = T, cex.studlab = 0.8, pos.studlab = 4)
dev.off()




# Egger's test

metabias(cv, method.bias = "linreg")

metabias(mace, method.bias = "linreg")

metabias(mi, method.bias = "linreg")

metabias(se, method.bias = "linreg")

metabias(is, method.bias = "linreg", k.min=9)

metabias(hs, method.bias = "linreg", k.min = 9)

metabias(ich, method.bias = "linreg")

metabias(gi, method.bias = "linreg")

metabias(mb, method.bias = "linreg")


sens2 <- sens[complete.cases(sens), ]







group <- c(0, 0,0,1,1,1)
outcome <- c("Benefit (MACE)", "Harm (Major bleeding)", "net RD", "Benefit (MACE)", "Harm (Major bleeding)", "net RD")
value <- c(-3.97, 4.66, 0.69, -3.31, 11.35, 8.04)

data <- data.frame(group, outcome, value)
str(data)
data$group <- factor(data$group, levels = c(0,1), labels = c("Western population", "East Asian population"))


library(ggplot2)

figure_net <- ggplot(data, aes(x=group, y=value, fill = outcome)) + geom_bar(stat="identity", position = position_dodge(0.8), width = 0.7) + theme_classic() +
  scale_y_continuous(name="Risk Difference (per 1000 persons)", breaks = c(-5:12,2))+
  theme(legend.position = "top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=12, face="bold"))+
  scale_fill_manual(values = c("#56B4E9", "#D55E00", "#009E73"))+
  geom_hline(yintercept=0) +
  geom_text(aes(label=value), color="black",position = position_dodge(0.8),size=4.5)+
  geom_vline(xintercept=1.5, colour="grey", lty="dashed", size=1)+
  annotate("segment", x=1.5, xend=1.5, y=1, yend=3.5, size=1.0, colour="black", arrow=arrow())+
  annotate("segment", x=1.5, xend=1.5, y=-1, yend=-3.5, size=1.0, colour="black", arrow=arrow())+
  annotate("text", x=1.5, y=-4.0, label="Benefit", size=4)+
  annotate("text", x=1.5, y=4.0, label="Harm", size=4)+
  annotate("text", x=1.2, y=13.0, label="net NNT = 1449", size=4)+
  annotate("text", x=2.2, y=13.0, label="net NNT = 124", size=4)
  
  

ggsave("D:/csv/figure_net.emf", plot=figure_net)







funnel.meta(m2, xlim = c(0.5, 4), studlab = TRUE)


metabias(m1, method.bias = "linreg")

tf <- trimfill(m1)
summary(tf)



m1 <- metabin(MB_A1, Total_A, MB_P1, Total_P,
              data = sens, sm = "OR")



tf1
funnel(tf1)
funnel(tf1, pch = ifelse(tf1$trimfill, 1, 16), level = 0.9, random = T)



data(Olkin1995)
m1 <- metabin(ev.exp, n.exp, ev.cont, n.cont,
              data = Olkin1995, subset = c(41, 47, 51, 59),
              studlab = paste(author, year),
              sm = "RR", method = "I")

oldpar <- par(mfrow = c(2, 2))

# Standard funnel plot
#
funnel(m1)

# Funnel plot with confidence intervals, common effect estimate and
# contours
#
cc <- funnel(m1, common = TRUE,
             level = 0.95, contour = c(0.9, 0.95, 0.99))$col.contour
legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"), fill = cc)

# Contour-enhanced funnel plot with user-chosen colours
#
funnel(m1, common = TRUE,
       level = 0.95, contour = c(0.9, 0.95, 0.99),
       col.contour = c("darkgreen", "green", "lightgreen"),
       lwd = 2, cex = 2, pch = 16, studlab = TRUE, cex.studlab = 1.25)
legend(0.05, 0.05,
       c("0.1 > p > 0.05", "0.05 > p > 0.01", "< 0.01"),
       fill = c("darkgreen", "green", "lightgreen"))

par(oldpar)

#
# Use log odds ratios on x-axis
#
funnel(tf1, backtransf = FALSE)
funnel(tf1, pch = ifelse(tf1$trimfill, 1, 16), level = 0.9, random = FALSE,
       backtransf = FALSE)

trimfill(m1$TE, m1$seTE, sm = m1$sm)


library(metafor)

res <- rma(yi, vi, data=sens, measure="OR")

### carry out trim-and-fill analysis
taf <- trimfill(res)

### draw funnel plot with missing studies filled in
funnel(taf, legend=TRUE)




forest(metainf(m1))
forest(metainf(m1), layout = "revman5")
forest(metainf(m1, pooled = "random"))

metainf(m1, sortvar = study)
metainf(m1, sortvar = 7:1)

m2 <- update(m1, title = "Fleiss1993bin meta-analysis", backtransf = FALSE)
metainf(m2)

data(Fleiss1993cont)
m3 <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
               data = Fleiss1993cont, sm = "SMD")
metainf(m3)