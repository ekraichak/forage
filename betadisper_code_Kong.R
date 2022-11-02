us <- read.csv("C:\\Users\\Daonurai\\Desktop\\Andaman\\R1_UnderGrowth.csv", header = T)
before.fire <- us %>% filter(Plot != "106", Plot != "110", Fire.appear == "Fire")
names(before.fire)
mtrx.TR.IVI <- picante::sample2matrix(before.fire[, c(3, 13, 8)])
names(mtrx.TR.IVI)
nrow(mtrx.TR.IVI)

envi <- read.csv("C:\\Users\\Daonurai\\Desktop\\Andaman\\Environmental_all.csv", header = TRUE)
t.env <- envi %>% filter(Fire.appear == "Fire", ccpoint == "nocc", Plot != "106", Plot != "110")
names(t.env)
nrow(t.env)

set.seed(7)
v.dist.c <- vegdist(mtrx.TR.IVI, method = "bray", na.rm = T) 
dispersion <- betadisper(v.dist.c, group = t.env$Fire)
permutest(dispersion)

## S3 method for class 'betadisper'
df <- eigenvals(dispersion)
summary(dispersion)
sum(dispersion$eig)
## Percentage of axes
ss <- df/sum(df) * 100

pch = rep(16)
pch[t.env[, 7] == "before"] = 3
pch[t.env[, 7] == "after"] = 4

windowsFonts(font = windowsFont("Times New Roman"))
font <- par(family ="font")

plot(dispersion, axes = c(1, 2), cex = 0.5, col = c("#2E5440", "#D2A010", "#762a83"),
     lty = "solid", lwd = 3, 
     hull = FALSE, ellipse = TRUE, pch = c(3, 8, 4, 10), #pch = text(dispersion$vectors)
     segments = TRUE, seg.col = brewer.pal(4, "Pastel2"), seg.lty = 1, seg.lwd = 2,
     label = F, label.cex = 5, main = "Dispersion",
     xlab = paste("PCoA1 - ", round(ss[1], 2), "%", sep = ""),
     ylab = paste("PCoA2 - ", round(ss[2], 2), "%", sep = ""),
     xlim = c(-0.5, 0.5), ylim = c(-0.7, 0.4)) ##sd ellipse

leg.text <- c("before", "after")
legend("bottomright", leg.text , pch = c(3, 8, 4, 10), bty = "n", cex = 1.1,
       col = c("#2E5440", "#D2A010"), 
       text.col = c("#2E5440", "#D2A010"))