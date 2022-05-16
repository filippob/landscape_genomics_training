# 0. Libraries and ad hoc functions required ------------------------------
# library("devtools")
# devtools::install_github("https://github.com/royfrancis/pophelper")
# BiocManager::install("qvalue")
# install.packages("rgdal")

require("here")
require("LEA")
require("maps")
require("lfmm")
require("qqman")
require("vegan")
require("rgdal")
require("qvalue")
require("StAMPP")
require("scales")
require("raster")
require("robust")
require("ggplot2")
require("bigsnpr")
require("parallel")
require("fuzzySim")
require("corrplot")
require("adegenet")
require("pophelper")
require("gridExtra")
require("data.table")
require("ggVennDiagram")

workdir <- here("scripts")
message("Moving to the working directory: ", workdir)
setwd(workdir)

## support functions
source("support_functions.r")

# 1. Quality control of the molecular dataset -----------------------------

fam_raw <- read.table("goat.fam")
table(fam_raw$V1)

####
## Quality control with Plink 2.0
####
plink2 <- download_plink2(dir = "./", overwrite = FALSE, verbose = TRUE)

# You after the first call, you have to clean-up results in order to call
# snp_plinkQC again
tryCatch(
  {snp_plinkQC(

    # path to plink
    plink.path = plink2,

    # prefix for the input file
    prefix.in = "./goat",

    # type of the input file (binary)
    file.type = "--bfile",

    # prefix for the output file
    prefix.out = "./goat_qced",

    # threshold for the minor allele frequency
    maf = 0.1,

    # missingness allowed per SNP
    geno = 0.01,

    # missingness allowed per individual
    mind = 0.01,

    # Hardy-Weinberg equilibrium exact test
    hwe = 1e-50,

    # should the analysis be restricted to autosomes only?
    # since the default is human, we need to set this option to FALSE for not
    # to remove chromosomes >22
    autosome.only = FALSE,

    # options for plink (KING-robust kinship estimator to screen for related samples
    # 0.0884 is to remove up to second-degree relations)
    extra.options = "--chr-set 29 --allow-extra-chr --king-cutoff 0.0884",
    verbose = TRUE
  )}, error = function(e) {
    warning(
      paste("You have already called snp_plinkQC once, Please remove ",
            "goat_qced output files before calling this function again: ", e)
      )
  }
)

fam_qced <- read.table("goat_qced.fam")
table(fam_qced$V1)

####
## Associating each individual to its geographical coordinates
####

geo_coord <- read.table("geo_coordinates.txt",header = T)
View(geo_coord)
geo_coord <- geo_coord[match(fam_qced$V2, geo_coord$smarter_id), ]

# 2. Population structure  ------------------------------------------------

####
## Preparing the right format
####

system(paste0(plink2, " --bfile goat_qced --export A --chr-set 29 --allow-extra-chr --out goat_qced"))
goat_qced <- fread("./goat_qced.raw", header = T)
write.table(goat_qced, "goat_qced.raw", col.names=T, row.names=F, quote=FALSE, sep=" ")
popstruct_input <- read.PLINK("goat_qced.raw", parallel=F, sep="\t")

goat_qced <- goat_qced[, -c(1:6)]
goat_qced <- as.matrix(goat_qced)
write.lfmm(goat_qced, "goat_qced.lfmm")
goat_qced <- lfmm2geno("goat_qced.lfmm")

####
## Principal component analysis (PCA)
####

pca <- glPca(popstruct_input, loadings = F)
print(pca)
eig.perc <- 100*pca$eig/sum(pca$eig)
head(cumsum(eig.perc))

# PCA scatterplot
scatter(pca, posi="topright", xax = 1, yax = 2, label = popstruct_input@pop, )
title("PCA axes 1-2")

jpeg("01-PCA-1vs2.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(pca$scores[, 1:2], t="n", xlab=paste0("PC1 (", round(eig.perc[1], 1), "%)"), ylab=paste0("PC2 (", round(eig.perc[2], 1), "%)"))
abline(h=0, v=0, col="gray")
myCol <- colorplot(pca$scores,pca$scores, add.plot=T, transp=T, alpha=0.5, xaxt="n", yaxt="n", cex=0)
text(x=pca$scores[, 1], y=pca$scores[, 2], labels=as.character(popstruct_input@pop), col=myCol)
add.scatter.eig(pca$eig[1:10],2,1,2, posi="topright", inset=.05, ratio=.3)
dev.off()

jpeg("02-PCA-1vs3.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(pca$scores[, c(1,3)], t="n", xlab=paste0("PC1 (", round(eig.perc[1], 1), "%)"),
     ylab=paste0("PC3 (", round(eig.perc[3], 1), "%)"))
abline(h=0, v=0, col="gray")
text(x=pca$scores[, 1], y=pca$scores[, 3], labels=as.character(popstruct_input@pop), col=myCol)
add.scatter.eig(w = pca$eig[1:10], xax = 1, yax = 3, posi="topright", inset=.05, ratio=.3)
dev.off()

jpeg("03-PCA-map.jpeg", width = 7, height = 7, units = 'in', res = 800)
map("italy", fill = T, col = "lightgray", border="gray")
points(geo_coord$longitude, geo_coord$latitude, pch=geo_coord$symbol,
       col=myCol, cex=1.5)
legend("bottomleft", pch=c(0,1,2,5), legend=c("ORO","BIO","ARG", "ASP"))
title(main = paste("PCA (", round(sum(eig.perc[1:3]), 1), "% var. explained)", sep=""))
dev.off()

####
## Discriminant analysis of principal components (DAPC)
####

# K-means analysis on the principal components
# Max K to test in K-means analysis
maxk <- 6
grp <- find.clusters(popstruct_input, pca.select = "percVar",
                     perc.pca = 99, max.n.clust = maxk, choose.n.clust = TRUE)

# 1st DAPC run to individuate the optimal number of principal components to retain
# for not to incur into overfitting issues
dapc <- dapc(popstruct_input, grp$grp, pca.select = "percVar", perc.pca = 99, n.da = 1)
print(dapc)

jpeg("04-DAPC-ascore.jpeg", width = 7, height = 7, units = 'in', res = 800)
ascore <- optim.a.score(dapc)
dev.off()

# 2nd DAPC run with optimal number of PCs retained
dapc <- dapc(popstruct_input, grp$grp, n.pca=ascore$best, n.da=length(grp$size) - 1)
print(dapc)

jpeg("05-DAPC-discriminant function.jpeg", width = 7, height = 7, units = 'in', res = 800)
scatter(dapc, bg = "white", legend = F, cleg = 0.6, solid = 0.4, col="lightgray")
dev.off()

# Plotting DAPC results on the map
geo_coord$grp <- dapc$grp
geo_coord$grp_symbol <- NA
geo_coord$grp_symbol[which(dapc$grp == 1)] <- 15
geo_coord$grp_symbol[which(dapc$grp == 2)] <- 17
dapc_col <- colorRampPalette(c("red", "gold", "lightblue", "blue"))
dapc_col <- dapc_col(10)[as.numeric(cut(dapc$ind.coord, breaks = 10))]

jpeg("06-DAPC-map.jpeg", width = 7, height = 7, units = 'in', res = 800)
map("italy", fill = T, col = "lightgray", border="gray")
points(geo_coord$longitude, geo_coord$latitude, pch=geo_coord$grp_symbol,
       col=dapc_col, cex=1.5)
legend("bottomleft", legend = sort(as.character(unique(geo_coord$grp))),
       pch=sort(unique(geo_coord$grp_symbol)))
title(main = "DAPC")
dev.off()

####
##  Fst estimation
####
popstruct_input@pop <- dapc$assign
fst <- stamppFst(popstruct_input, nboots=100, percent=95, nclusters=detectCores()-1)
print(fst)

####
##  Gene flow estimation
####

# The calculation is based on the Wright's formula (Wright 1949),
# where N is the effective population size of the concerned population, m the
# migration rate, and Fst the fixation index between the pair of populations considered.
# Nm is an indirect estimate of gene flow. It can be interpreted as the number of
# individuals migrating into a population each generation that actually contribute
# to its gene pool. It reflects the relative balance between gene flow and genetic drift and,
# more specifically, the current/recent importance of gene flow relative to drift
# (Slatkin & Barton, 1989). Generally, Nm>1 is considered as an evidence of
# sufficient gene flow between the populations to prevent substantial differentiation
# due to genetic drift. Uncritical applications of the Wright's formula discouraged
# in contexts where departures from th"e implicit assumptions of the infinite island model
# are expected to occur (Whitlock & Mccauley, 1999).
Nm(Fst = fst$Fsts[2,1])

####
## Sparse non-negative matrix factorization
####

# Input file for sNMF (.geno)
obj.snmf <- snmf(
  goat_qced,
  K=1:maxk,
  repetitions=5,
  entropy=T,
  ploidy=2,
  project="new",
  CPU = detectCores()-1
)

snmf_entropy <- calculate_snmf_entropy(snmf_res = obj.snmf, maxk = maxk)
head(snmf_entropy)

# Cross-entropy
jpeg("07-sNMF-cross-entrpy.jpeg", width = 7, height = 7, units = 'in', res = 800)
boxplot(CE~K, data=snmf_entropy,
        xlab="N.er of ancestral populations", ylab="Cross-entropy")
dev.off()

jpeg("08-sNMF-barplot.jpeg", width = 7, height = 7, units = 'in', res = 800)
snmf_barplot(snmf_res = obj.snmf, K = 2, dapc_assign = dapc$assign, pop = fam_qced$V1)
dev.off()

# 3. Environmental dataset ------------------------------------------------

env <- list.files(path = "./", pattern = ".tif", all.files=TRUE, full.names=FALSE)
env <- stack(paste0("./", env))
print(env[[1]])
plot(env[[1]])

jpeg("09-Envar-map.jpeg", width = 10, height = 10, units = 'in', res = 800)
plot(env)
dev.off()

# extract environmental information in correspondance of the individual coordinates
env <- extract(env, geo_coord[, c("longitude", "latitude")])
View(env)

# investigating collinearity
env_cor <- cor(env)

jpeg("10-Envar-cor.jpeg", width = 7, height = 7, units = 'in', res = 800)
corrplot(
  env_cor, order = "original", type = "upper", diag=T, tl.cex = 1,
  tl.col = c(rep("red", 11), rep("blue", 8), "forestgreen"), addCoef.col = "white",
  number.cex = 0.6, number.font = 1
)
dev.off()

# Variance inflation factor (VIF) associated to each variable in the dataset
env_vif <- multicol(vars = as.data.frame(env))
env_vif

# Threshold for VIF
vif.thr <- 10

# Obtaining a dataset with max(VIF) <vif.thr
while (length(which(env_vif$VIF > vif.thr)) >= 1) {
  t <- rownames(env_vif)[1]
  t <- which(colnames(env) == t)
  env <- env[, -t]
  env_vif <- multicol(vars = as.data.frame(env))
  print(env_vif); cat("\n")
}

env_cor <- cor(env)
jpeg("11-Envar-selected.jpeg", width = 7, height = 7, units = 'in', res = 800)
corrplot(
  env_cor, order = "original", type = "upper", diag=T, tl.cex = 1.5,
  tl.col = c(rep("red", 2), rep("blue", 2), "forestgreen"), addCoef.col = "white",
  number.cex = 1.5, number.font = 1
)
dev.off()

jpeg("12-Envar-selected-hist.jpeg", width = 7, height = 7, units = 'in', res = 800)
par(mfrow=c(2, 3))
for (i in 1:ncol(env)) {
  hist(env[,i], main="", xlab=colnames(env)[i])
  title(colnames(env)[i])
}
dev.off()

# 4. Gene-environment association analysis --------------------------------

####
## LFMM analysis
####

# imputation of missing genotypes (required prior to running LFMM)
# Here, not necessary because we are working with 0% missing data
# best <- which.min(cross.entropy(obj.snmf, K=2))
# impute(obj.snmf, "goat_qced.lfmm", method="mode", K=2, run=best)

# Input files and parameters
Y <- read.lfmm("goat_qced.lfmm")
K <- 2

# Latent factors mixed models
lfmm <- lfmm_ridge(Y=Y, X=env, K=K)
str(lfmm)

# p-values for all SNP-environment associations
lfmm_pvalues <- lfmm_test(Y=Y, X=env, lfmm=lfmm, calibrate="gif")
lfmm_pvalues <- lfmm_pvalues$calibrated.pvalue
rownames(lfmm_pvalues) <- read.table("./goat_qced.bim",h=F,stringsAsFactors=F)$V2
dim(lfmm_pvalues);head(lfmm_pvalues)

# Histograms of p-values
jpeg("13-Lfmm-pvlaues.jpeg", width = 7, height = 5, units = 'in', res = 800)
par(mfrow = c(2,3), mar = c(3.5,3.5,3,0.5), mgp = c(2.5,0.8,0))
for (i in 1:ncol(lfmm_pvalues)) {
  hist(lfmm_pvalues[, i], breaks=20, xlab="p-values", xlim=c(0, 1),
       main=colnames(lfmm_pvalues)[i], col="darkgray", border="darkgray")
}
dev.off()

# Quantile-Quantile plots
jpeg("14-Lfmm-qqplot.jpeg", width = 7, height = 5, units = 'in', res = 800)
par(mfrow = c(2,3), mar = c(3.5,4,3,0.5), mgp = c(2.5,0.8,0))
for (i in 1:ncol(lfmm_pvalues)) {
  qqplot(rexp(length(lfmm_pvalues[,i]), rate=log(10)), -log10(lfmm_pvalues[,i]),
         xlab="Expected quantile", ylab=expression(-log[10] ~ p-values),
         pch=19, cex=0.4, main=colnames(lfmm_pvalues)[i])
  abline(0, 1)
}
dev.off()

####
## false discovery rate (FDR) control
####

# q-values calculation
qobj <- qvalue(lfmm_pvalues[, "ELEV"])
hist(qobj) +
  ggtitle(paste0("p-value density histogram - ", colnames(lfmm_pvalues)[5]))
plot(qobj)

# q-values need to be computed by environmental variable
lfmm_qvalues <- lfmm_qvalue(pv = lfmm_pvalues, bim = "./goat_qced.bim")
head(lfmm_qvalues)

# Histograms of q-values
jpeg("15-Lfmm-qvlaues.jpeg", width = 7, height = 5, units = 'in', res = 800)
par(mfrow = c(2,3), mar = c(3.5,3.5,3,0.5), mgp = c(2.5,0.8,0))
for (i in 1:ncol(lfmm_pvalues)) {
  hist(lfmm_qvalues[, i], breaks=20, xlab="q-values", xlim=c(0, 1),
       main=colnames(lfmm_qvalues)[i], col="darkgray", border="darkgray")
}
dev.off()

# LFMM results
lfmm_res <- lfmm_qvalue_cut(qvalues=lfmm_qvalues, nenv=5, cutoff=0.2, bim.info=T)
lfmm_res

# Manhattan plots
jpeg("16-Lfmm-manhattan.jpeg", width = 7, height = 5, units = 'in', res = 800)
par(mfrow = c(length(unique(lfmm_res$ENV)), 1),
    oma = c(1, 1, 1, 1), mar = c(4, 5, 2.7, 1),
    mgp = c(2.5, 0.8, 0))

for (i in 1:length(unique(lfmm_res$ENV))) {
  P <- lfmm_pvalues[, which(unique(lfmm_res$ENV)[i] == colnames(lfmm_pvalues))]
  man.df <- cbind.data.frame(P, lfmm_qvalues[, c("SNP", "CHR", "BP")])
  manhattan(man.df, genomewideline=F, suggestiveline=F,
            main=paste0("Lfmm (K=", K,") - ", colnames(lfmm_pvalues)[which(unique(lfmm_res$ENV)[i] == colnames(lfmm_pvalues))]),
            highlight=lfmm_res$SNP[which(lfmm_res$ENV==unique(lfmm_res$ENV)[i])],
            col=c(alpha("gray70", 0.7), alpha("gray30", 0.7)), xlab = "Chromosome",
            cex.lab=1.5, cex.axis=1, cex.main=1.5, font.main=3)
  box()
}
dev.off()

####
## Redundancy analysis
####

# Defining the population structure variables from LFMM
PS <- scale(lfmm$U)
PS <- as.data.frame(PS)
colnames(PS) <- c("PS1", "PS2")
rownames(PS) <- read.table("goat_qced.fam")$V2
head(PS)

# response matrix
colnames(Y) <- read.table("goat_qced.bim")$V2
rownames(Y) <- read.table("goat_qced.fam")$V2
Y[1:6, 1:3]

# Standardization of the environmental variables
env <- scale(env, center=TRUE, scale=TRUE)

# Recovering scaling coefficients (for later use)
scale_env <- attr(env, 'scaled:scale')
center_env <- attr(env, 'scaled:center')

# Environmental matrix
env <- as.data.frame(env)
row.names(env) <- read.table("goat_qced.fam")$V2
head(env)

# Table gathering all variables
Variables <- data.frame(PS, env)
head(Variables)
Variables_cor <- cor(Variables)

jpeg("17-pRDA-predictor-cor.jpeg", width = 7, height = 5, units = 'in', res = 800)
corrplot(
  Variables_cor, order = "original", type = "upper",
  diag=T, tl.cex = 0.7, tl.col = c(
    rep("gray10", 2), rep("red", 2), rep("blue", 2),
    "forestgreen"), addCoef.col = "white", number.cex = 0.7, number.font = 1
)
dev.off()

# Genome scan using pRDA
RDA_env <- rda(Y ~ BIO03 + BIO08 + BIO13 + BIO15 + ELEV + Condition(PS1 + PS2),  Variables)
print(RDA_env)

# adjusting the proportion of variance explained by the environmental predictors
RsquareAdj(RDA_env)

# eigenvalues for the constrained axes
summary(eigenvals(RDA_env, model = "constrained"))
jpeg("18-pRDA-screeplot.jpeg", width = 7, height = 5, units = 'in', res = 800)
screeplot(RDA_env, main="Eigenvalues of constrained axes")
dev.off()

# Statistical significance of the RDA model using F-statistics
# Null hypothesis: no linear relationship exists between the SNP data and the environmental predictors
signif.full <- anova.cca(RDA_env, parallel=detectCores()-1) # default is permutation=999
signif.full
# Full model is significant

# Now, we can check each constrained axis for significance to determine which
# constrained axes we should investigate for candidate loci (for this test,
# each constrained axis is tested using all previous constrained axes as conditions)
# please do not run this line it will take a very long time!
# signif.axis <- anova.cca(RDA_env, by="axis", parallel=getOption("mc.cores"))
# signif.axis
# We should keep the first three axes

# Plot the RDA: SNPs are in red, individuals in black, the blue vectors are the environmental predictors.
# The relative arrangement of these items in the ordination space reflects their
# relationship with the ordination axes, which are linear combinations of the predictor variables.
jpeg("19-pRDA-1vs2.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(RDA_env, type="n", scaling=3)
# the SNPs
points(RDA_env, display="species", pch=20, cex=0.7, col="gray32", scaling=3)
# the individuals
points(RDA_env, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=geo_coord$color)
# the predictors
text(RDA_env, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topleft", legend=c("Orobica", "Bionda dell'Adamello", "Argentata", "Aspromontana"), bty="n", col="gray32", pch=21, cex=1, pt.bg=c("royalblue", "forestgreen", "orange", "gold"))
dev.off()

# axes 1 & 3
jpeg("20-pRDA-1vs3.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(RDA_env, type="n", scaling=3, choices=c(1,3))
points(RDA_env, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices=c(1,3))
points(RDA_env, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=geo_coord$color, choices=c(1,3))
text(RDA_env, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("Orobica", "Bionda dell'Adamello", "Argentata", "Aspromontana"), bty="n", col="gray32", pch=21, cex=1, pt.bg=c("royalblue", "forestgreen", "orange", "gold"))
dev.off()

# Candidate SNPs involved in local adaptation can be identified by their loading
# in the ordination space. Here, we will derive the SNP loadings from the three significant
# constrained axes only:
load.rda <- as.data.frame(scores(RDA_env, choices=c(1:3), display="species"))
load.rda$SNP <- read.table("goat_qced.bim")$V2
load.rda$CHR <- read.table("goat_qced.bim")$V1
load.rda$BP <- read.table("goat_qced.bim")$V4
head(load.rda)

jpeg("21-pRDA-loadings-hist.jpeg", width = 5, height = 7, units = 'in', res = 800)
par(mfrow=c(3,1))
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
dev.off()

# Identification based on +/- 3.5 SD from the mean
cand1 <- rda_outliers(load_rda = load.rda, n_rda = 1, n_sign = 3, z = 3.5);nrow(cand1)
cand2 <- rda_outliers(load_rda = load.rda, n_rda = 2, n_sign = 3, z = 3.5);nrow(cand2)
cand3 <- rda_outliers(load_rda = load.rda, n_rda = 3, n_sign = 3, z = 3.5);nrow(cand3)
ncand <- nrow(cand1) + nrow(cand2) + nrow(cand3)
ncand

# Correlation between each candidate SNP and the environmental variables
# and return the name of the environmental variable most correlated
rda_res <- rda_outliers_env(cand1 = cand1, cand2 = cand2, cand3 = cand3, env = env, Y = Y)

# Plot the SNPs
# Assigning a color code to the candidate SNPs based on the association with the
# env. variables
cand_col <- rep(NA, nrow(rda_res))
cand_col[which(rda_res$ENV == "BIO03")] <- "orange"
cand_col[which(rda_res$ENV == "BIO08")] <- "red"
cand_col[which(rda_res$ENV == "BIO13")] <- "royalblue"
cand_col[which(rda_res$ENV == "BIO15")] <- "lightblue"
cand_col[which(rda_res$ENV == "ELEV")] <- "forestgreen"

# assign a color code to the neutral SNPs
snp_col <- rownames(RDA_env$CCA$v)
names(snp_col) <- snp_col
snp_col[1:length(snp_col)] <- "ghostwhite"

# axes 1 & 2
jpeg("22-pRDA-1vs2-SNPs.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(RDA_env, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
snp_col[1:length(snp_col)] <- "ghostwhite"
points(RDA_env, display="species", pch=21, cex=1, col="gray32", bg=snp_col, scaling=3)
snp_col[1:length(snp_col)] <- rgb(0,1,0, alpha=0)
for (i in 1:length(rda_res$SNP)) {
  snp_col[match(rda_res$SNP[i],names(snp_col))] <- cand_col[i]
}
points(RDA_env, display="species", pch=21, cex=1, col=rgb(0,1,0, alpha=0), bg=snp_col, scaling=3)
text(RDA_env, scaling=3, display="bp", col="#0868ac", cex=1)
legend("topleft", legend=c("BIO03","BIO08","BIO13","BIO15","ELEV"), bty="n", col="gray32",
       pch=21, cex=1, pt.bg=c("orange", "red", "royalblue", "lightblue", "forestgreen"))
dev.off()

# axes 1 & 3
jpeg("23-pRDA-1vs3-SNPs.jpeg", width = 7, height = 7, units = 'in', res = 800)
plot(RDA_env, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
snp_col[1:length(snp_col)] <- "ghostwhite"
points(RDA_env, display="species", pch=21, cex=1, col="gray32", bg=snp_col, scaling=3, choices=c(1,3))
snp_col[1:length(snp_col)] <- rgb(0,1,0, alpha=0)
for (i in 1:length(rda_res$SNP)) {
  snp_col[match(rda_res$SNP[i],names(snp_col))] <- cand_col[i]
}
points(RDA_env, display="species", pch=21, cex=1, col=rgb(0,1,0, alpha=0), bg=snp_col, scaling=3, choices=c(1,3))
text(RDA_env, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("topleft", legend=c("BIO03","BIO08","BIO13","BIO15","ELEV"), bty="n", col="gray32",
       pch=21, cex=1, pt.bg=c("orange", "red", "royalblue", "lightblue", "forestgreen"))
dev.off()

# 5. Adaptive landscape and genomic offset --------------------------------

# How many - if any - oultiers in common between pRDA and LFMM?
list_outliers <- list(pRDA=unique(rda_res$SNP), LFMM = unique(lfmm_res$SNP))
jpeg("24-Venn.jpeg", width = 7, height = 7, units = 'in', res = 800)
ggVennDiagram(list_outliers, category.names = c("partial RDA", "LFMM"), lty="solid", size=0.2) +
  scale_fill_gradient2(low = "white", high = 'gray40') +
  scale_color_manual(values = c("grey", "grey", "grey", "grey")) +
  guides(fill = "none")
dev.off()

# Shared outliers between RDA and LFMM
shared_outliers <- table(c(unique(rda_res$SNP), LFMM = unique(lfmm_res$SNP)))
shared_outliers <- names(shared_outliers[which(shared_outliers > 1)])
res <- data.frame(shared_outliers=shared_outliers)
for (i in 1:nrow(res)) {
  tmp1 <- which(res$shared_outliers[i] == rda_res$SNP)
  tmp1 <- rda_res$ENV[tmp1]
  if(length(tmp1>1)) tmp1 <- paste(tmp1, collapse = ",")
  tmp2 <- which(res$shared_outliers[i] == lfmm_res$SNP)
  tmp2 <- lfmm_res$ENV[tmp2]
  if(length(tmp2>1)) tmp2 <- paste(tmp2, collapse = ",")
  res$pRDA[i] <- tmp1
  res$LFMM[i] <- tmp2
}
print(res)

# Adaptively enriched RDA
RDA_outliers <- rda(Y[,res$shared_outliers] ~ BIO08 + BIO15 + ELEV,  Variables)
anova.cca(RDA_outliers, by="axis", parallel=detectCores()-1)

# Adaptive index
adaptive_landscape <- adaptive_index(
  RDA = RDA_outliers, K = 2,
  env_pres = stack(raster("BIO08.tif"), raster("BIO15.tif"), raster("ELEV.tif")),
  method = "loadings", scale_env = scale_env[c("BIO08","BIO15","ELEV")],
  center_env = center_env[c("BIO08","BIO15","ELEV")]
  )

jpeg("25-Adaptive landscape.jpeg", width = 8, height = 3, units = 'in', res = 800)
par(mfrow=c(1, 3))
plot(RDA_outliers, type="n", scaling=3)
points(RDA_outliers, display="species", pch=20, cex=1.7, col="gray32", scaling=3)
text(RDA_outliers, scaling=3, display="bp", col="#0868ac", cex=1)
plot(adaptive_landscape$RDA1, legend=F, axes=FALSE, bty="n", col=topo.colors(100))
legend("topright", legend="RDA1", bty="n", cex = 2)
plot(adaptive_landscape$RDA2, legend=FALSE, axes=FALSE, bty="n", col=topo.colors(100))
legend("topright", legend="RDA2", bty="n", cex = 2)
mtext(text = "Adaptive index - present conditions", side = 3, line = -2, outer = T, at = 0.5, cex=1.5, font=2)
dev.off()

# Genomic offset
ras_2080 <- stack(raster("BIO08-2080.grd"), raster("BIO15-2080.grd"), raster("ELEV.tif"))
names(ras_2080) <- c("BIO08", "BIO15", "ELEV")
adaptive_landscape_2080 <- adaptive_index(
  RDA = RDA_outliers, K = 2, env_pres = ras_2080, method = "loadings",
  scale_env = scale_env[c("BIO08","BIO15","ELEV")], center_env = center_env[c("BIO08","BIO15","ELEV")]
)

genomic_offset_2080_RDA1 <- (adaptive_landscape_2080$RDA1 - adaptive_landscape$RDA1)^2
genomic_offset_2080_RDA2 <- (adaptive_landscape_2080$RDA2 - adaptive_landscape$RDA2)^2

jpeg("26-Genomic offset.jpeg", width = 12, height = 7, units = 'in', res = 800)
par(mfrow=c(1, 2), mar=c(3, 0, 5, 5))
plot(genomic_offset_2080_RDA1, legend=T, axes=FALSE, col=gray.colors(10^4))
contour(genomic_offset_2080_RDA1, levels=c(10^-3.5, 2.5, 7),
        labels=NULL, add=T, labcex=0, col=c("green", "gold", "red"))
legend("topright", title = "RDA1", legend=c("7","2.5","~0"), lty=1,
       col=c("red", "gold", "green"),
       bty="n", cex = 1.5)
plot(genomic_offset_2080_RDA2, legend=T, axes=FALSE, col=gray.colors(10^4))
contour(genomic_offset_2080_RDA2, levels=c(10^-3.5, 0.5, 2), labels=NULL, add=T, labcex=0,
        col=c("green", "gold", "red"))
legend("topright", title = "RDA2", legend=c("2","0.5","~0"), lty=1, col=c("red", "gold", "green"),
       bty="n", cex = 1.5)
mtext(text = "Genomic offset - 2080", side = 3, line = -3, outer = T, at = 0.5, cex=2.5, font=1)
dev.off()
