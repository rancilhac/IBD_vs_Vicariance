##### 1 - function adapted to perform IBD tests with single localities

# a function to test whether two linear regressions between genetic and geographic distances are statistically similar
# Loïs Rancilhac, 01/22
# Adapted from prabclus::regeqdist (Hausdorf & Hennig 2020) 


diffreg <- function(dmx1 , dmx2 , dmy1, dmy2){
  # dmx1, dmx2 = matrices of geographic distances 1 & 2
  # dmy1, dmy2 = matrices of genetic distances 1 & 2
  
  dmxc <- dmyc <- jr <- lmfit <- xvi <- yvi <- list()
  nc <- sediff <- coefdiff <- pval <- condition <- numeric(0)
  dmxc[[1]] <- dmx1
  dmxc[[2]] <- dmx2
  dmyc[[1]] <- dmy1
  dmyc[[2]] <- dmy2
  
  groups <- c(1,2)
  
  for (i in 1:2) {
    jr[[i]] <- list()
    # N ind in matrix i
    nc[i] <- sum(dim(dmxc[[i]])[1])
    # sous matrice geo en vecteur
    xvi[[i]] <- as.vector(as.dist(dmxc[[i]]))
    # sous matrice genet en vecteur
    yvi[[i]] <- as.vector(as.dist(dmyc[[i]]))
  }
  
  xall <- c(xvi[[1]], xvi[[2]])
  xcenter <- mean(xall)
  
  clm <- jackpseudo <- jackestcl <- jackvarcl <- list()
  jackse <- jackest <- tstat <- tdf <- numeric(0)
  for (i in 1:2) {
    jackpseudo[[i]] <- list()
    jackestcl[[i]] <- jackvarcl[[i]] <- numeric(0)
    for (j in 1:2) jackpseudo[[i]][[j]] <- numeric(0)
  }
  for (i in 1:2) {
    #sous matrice geo pondérée
    xvi[[i]] <- xvi[[i]] - xcenter
    lmfit[[i]] <- lm(yvi[[i]] ~ xvi[[i]])
    mm <- model.matrix(~xvi[[i]])
    condition[i] <- kappa(mm)
    clm[[i]] <- coef(lmfit[[i]])
  }
  for (i in 1:2) for (j in 1:2) jr[[i]][[j]] <- bootstrap::jackknife(1:sum(dim(dmxc[[i]])[1]), regdist, dmx = dmxc[[i]], dmy = dmyc[[i]], xcenter = xcenter, param = j)
  
  for (j in 1:2) {
    for (i in 1:2) if (is.na(jr[[i]][[j]]$jack.se)) 
      computable <- FALSE
    coefdiff[j] <- clm[[1]][j] - clm[[2]][j]
    for (i in 1:2) {
      for (k in 1:nc[i]) jackpseudo[[i]][[j]][k] <- nc[i] * 
          clm[[i]][j] - (nc[i] - 1) * jr[[i]][[j]]$jack.values[k]
      jackestcl[[i]][j] <- mean(jackpseudo[[i]][[j]])
      jackvarcl[[i]][j] <- var(jackpseudo[[i]][[j]])
    }
    jackest[j] <- jackestcl[[1]][j] - jackestcl[[2]][j]
    jackse[j] <- sqrt(jackvarcl[[1]][j]/nc[1] + jackvarcl[[2]][j]/nc[2])
    tstat[j] <- jackest[j]/jackse[j]
    tdf[j] <- jackse[j]^4/((jackvarcl[[1]][j]/nc[1])^2/(nc[1] - 
                                                          1) + (jackvarcl[[2]][j]/nc[2])^2/(nc[2] - 1))
    if (tstat[j] > 0) 
      pval[j] <- 2 * (pt(tstat[j], tdf[j], lower.tail = FALSE))
    else pval[j] <- 2 * (pt(tstat[j], tdf[j]))
  }
  
  out <- list(pval = pval, coefdiff = coefdiff, condition = condition, 
              lmfit = lmfit, jr = jr, xcenter = xcenter, tstat = tstat, 
              tdf = tdf, jackest = jackest, jackse = jackse, jackpseudo = jackpseudo, 
              groups = groups)
  class(out) <- "regeqdist"
  out
}

##### 2- Tests to detect barriers to gene flow for inland lineages (for more details about the groups, cf. Table 3 & 4 of the manuscript)
##### Matrices of geographic and genetic distances are available in the same github repository

library(prabclus)

####### 2.1 - Grp1 = N-Atl Grp2 = Rif

geographic_coordinates <- read.csv("1-Coordinates_AtlN_Rif.csv", sep="\t", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

ncol <- max(count.fields("./1-Genet_dist_AtlN_Rif.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./1-Genet_dist_a_AtlN_Rif.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="1-9_loci_Alleles_AtlN_Rif.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

# Create plots and do reg tests using a distances

genet_dist_within_AtlN_a <- genetic_distances_a[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Rif_a <- genetic_distances_a[which(popmap == 2), which(popmap == 2)]
geo_dist_within_AtlN <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Rif <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between <- genetic_distances_a[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_a, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_a, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

plot(genet_dist_within_AtlN_a ~ geo_dist_within_AtlN, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Rif_a ~ geo_dist_within_Rif, pch=19, col="red")
abline(lm(c(genet_dist_within_AtlN_a) ~ c(geo_dist_within_AtlN)))
abline(lm(c(genet_dist_within_Rif_a) ~ c(geo_dist_within_Rif)), col="red")
points(genet_dist_between ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between) ~ c(geo_dist_between)), col="blue")

# Create plots and do reg tests using shared alleles distances

genet_dist_within_AtlN_s <- genetic_distances_s[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Rif_s <- genetic_distances_s[which(popmap == 2), which(popmap == 2)]
geo_dist_within_AtlN <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Rif <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_s <- genetic_distances_s[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_s, grouping=popmap, groups=c(1,2))
pv <- H01_reg$pval[2]*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_s, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H031_reg$pval, "\n", "p-value H032 ", H032_reg$pval, sep="")}

plot(genet_dist_within_AtlN_s ~ geo_dist_within_AtlN, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Rif_s ~ geo_dist_within_Rif, pch=19, col="red")
abline(lm(c(genet_dist_within_AtlN_s) ~ c(geo_dist_within_AtlN)), lty=3, lwd=3)
abline(lm(c(genet_dist_within_Rif_s) ~ c(geo_dist_within_Rif)), col="red", lty=3, lwd=3)
points(genet_dist_between_s ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue", lwd=2)

######## 2.2 - Grp1 = Rif Grp2 = MidAtl

# Import geographic coordinates and calculate distances

geographic_coordinates <- read.csv("2-Coordinates_Rif_MidAtl.csv", sep=",", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

ncol <- max(count.fields("./2-Genet_dist_Rif_MidAtl.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./2-Genet_dist_a_Rif_MidAtl.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="2-9_loci_Alleles_Rif_MidAtl.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

# Create plots and do reg tests using a distances

genet_dist_within_MidAtl_a <- genetic_distances_a[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Rif_a <- genetic_distances_a[which(popmap == 2), which(popmap == 2)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Rif <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_a <- genetic_distances_a[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_a, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_a, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

plot(genet_dist_within_MidAtl_a ~ geo_dist_within_MidAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Rif_a ~ geo_dist_within_Rif, pch=19, col="red")
abline(lm(c(genet_dist_within_MidAtl_a) ~ c(geo_dist_within_MidAtl)))
abline(lm(c(genet_dist_within_Rif_a) ~ c(geo_dist_within_Rif)), col="red")
points(genet_dist_between_a ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between_a) ~ c(geo_dist_between)), col="blue")

# Create plots and do reg tests using s distances

genet_dist_within_MidAtl_s <- genetic_distances_s[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Rif_s <- genetic_distances_s[which(popmap == 2), which(popmap == 2)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Rif <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_s <- genetic_distances_s[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_s, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_s, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

plot(genet_dist_within_MidAtl_s ~ geo_dist_within_MidAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Rif_s ~ geo_dist_within_Rif, pch=19, col="red")
abline(lm(c(genet_dist_within_MidAtl_s) ~ c(geo_dist_within_MidAtl)), lty=3, lwd=3)
abline(lm(c(genet_dist_within_Rif_s) ~ c(geo_dist_within_Rif)), col="red", lty=3, lwd=3)
points(genet_dist_between_s ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue", lwd=2)

######## 2.3 - Grp1 = MidAtl Grp2 = Mara

# Import geographic coordinates and calculate distances

geographic_coordinates <- read.csv("3-Coordinates_MidAtl_Mara.csv", sep=",", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

ncol <- max(count.fields("./genet_dist_Mid-Atl_Mara.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./3-Genet_dist_MidAtl_Mara.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples
write.table(genetic_distances_a, "./Genet_dist_a_MidAtl_Mara.csv", sep=";")

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="3-9_loci_Alleles_MidAtl_Mara.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

# Create plots and do reg tests using a distances

genet_dist_within_MidAtl_a <- genetic_distances_a[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Mara_a <- genetic_distances_a[which(popmap == 2), which(popmap == 2)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Mara <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_a <- genetic_distances_a[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_a, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_a, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

plot(genet_dist_within_MidAtl_a ~ geo_dist_within_MidAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Mara_a ~ geo_dist_within_Mara, pch=19, col="red")
abline(lm(c(genet_dist_within_MidAtl_a) ~ c(geo_dist_within_MidAtl)))
abline(lm(c(genet_dist_within_Mara_a) ~ c(geo_dist_within_Mara)), col="red")
points(genet_dist_between_a ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between_a) ~ c(geo_dist_between)), col="blue")

# Create plots and do reg tests using s distances

genet_dist_within_MidAtl_s <- genetic_distances_s[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Mara_s <- genetic_distances_s[which(popmap == 2), which(popmap == 2)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Mara <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_s <- genetic_distances_s[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_s, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_s, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

plot(genet_dist_within_MidAtl_s ~ geo_dist_within_MidAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Mara_s ~ geo_dist_within_Mara, pch=19, col="red")
abline(lm(c(genet_dist_within_MidAtl_s) ~ c(geo_dist_within_MidAtl)), lty=3, lwd=3)
abline(lm(c(genet_dist_within_Mara_s) ~ c(geo_dist_within_Mara)), col="red", lty=3, lwd=3)
points(genet_dist_between_s ~ geo_dist_between, pch="+", col="blue", cex=1.5)
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue", lwd=2)

####### 2.4 Grp1 = Mara Grp2 = SAtl

geographic_coordinates <- read.csv("4-Coordinates_Mara_SAtl.csv", sep=",", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

ncol <- max(count.fields("./4-Genet_dist_a_SAtl_Mara.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./4-Genet_dist_a_SAtl_Mara.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="4-9_loci_Alleles_S_Atl_Mara.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

# Create plots and do reg tests using a distances

genet_dist_within_SAtl_a <- genetic_distances_a[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Mara_a <- genetic_distances_a[which(popmap == 2), which(popmap == 2)]
geo_dist_within_SAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Mara <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_a <- genetic_distances_a[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_a, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_a, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

plot(genet_dist_within_SAtl_a ~ geo_dist_within_SAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Mara_a ~ geo_dist_within_Mara, pch=19, col="red")
abline(lm(c(genet_dist_within_SAtl_a) ~ c(geo_dist_within_SAtl)))
abline(lm(c(genet_dist_within_Mara_a) ~ c(geo_dist_within_Mara)), col="red")
points(genet_dist_between_a ~ geo_dist_between, col="blue", pch="+", cex=1.5)
abline(lm(c(genet_dist_between_a) ~ c(geo_dist_between)), col="blue")

# Create plots and do reg tests using s distances

genet_dist_within_SAtl_s <- genetic_distances_s[which(popmap == 1), which(popmap == 1)]
genet_dist_within_Mara_s <- genetic_distances_s[which(popmap == 2), which(popmap == 2)]
geo_dist_within_SAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_Mara <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_s <- genetic_distances_s[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

H01_reg <- regeqdist(geographic_distances, genetic_distances_s, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_s, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}


plot(genet_dist_within_SAtl_s ~ geo_dist_within_SAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_Mara_s ~ geo_dist_within_Mara, pch=19, col="red")
abline(lm(c(genet_dist_within_SAtl_s) ~ c(geo_dist_within_SAtl)), lty=3, lwd=3)
abline(lm(c(genet_dist_within_Mara_s) ~ c(geo_dist_within_Mara)), col="red", lty=3, lwd=3)
points(genet_dist_between_s ~ geo_dist_between, col="blue", pch="+", cex=1.5)
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue", lwd=2)

####### 2.5 Grp1 = MidAtl Grp2 = SAtl

# Import geographic coordinates and calculate distances

geographic_coordinates <- read.csv("5-Coordinates_S-Atl_Mid-Atl.csv", sep=",", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

ncol <- max(count.fields("./5-Genet_dist_a_Mid-Atl_S-Atl.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./5-Genet_dist_a_Mid-Atl_S-Atl.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="5-9_loci_Alleles_SAtl_MidAtl.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

# Create plots and do reg tests using a distances

genet_dist_within_SAtl_a <- genetic_distances_a[which(popmap == 1), which(popmap == 1)]
genet_dist_within_MidAtl_a <- genetic_distances_a[which(popmap == 2), which(popmap == 2)]
geo_dist_within_SAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_a <- genetic_distances_a[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

plot(genet_dist_within_SAtl_a ~ geo_dist_within_SAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_MidAtl_a ~ geo_dist_within_MidAtl, pch=19, col="red")
abline(lm(c(genet_dist_within_SAtl_a) ~ c(geo_dist_within_SAtl)))
abline(lm(c(genet_dist_within_MidAtl_a) ~ c(geo_dist_within_MidAtl)), col="red")
points(genet_dist_between_a ~ geo_dist_between, pch="+", cex=1.5, col="blue")
abline(lm(c(genet_dist_between_a) ~ c(geo_dist_between)), col="blue")

H01_reg <- regeqdist(geographic_distances, genetic_distances_a, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_a, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

# Create plots and do reg tests using shared allele distances

genet_dist_within_SAtl_s <- genetic_distances_s[which(popmap == 1), which(popmap == 1)]
genet_dist_within_MidAtl_s <- genetic_distances_s[which(popmap == 2), which(popmap == 2)]
geo_dist_within_SAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 2), which(popmap == 2)]

genet_dist_between_s <- genetic_distances_s[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

plot(genet_dist_within_SAtl_s ~ geo_dist_within_SAtl, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
points(genet_dist_within_MidAtl_s ~ geo_dist_within_MidAtl, pch=19, col="red")
abline(lm(c(genet_dist_within_SAtl_s) ~ c(geo_dist_within_SAtl)), lty=3, lwd=3)
abline(lm(c(genet_dist_within_MidAtl_s) ~ c(geo_dist_within_MidAtl)), col="red", lty=3, lwd=3)
points(genet_dist_between_s ~ geo_dist_between, pch="+", cex=1.5, col="blue")
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue", lwd=2)

H01_reg <- regeqdist(geographic_distances, genetic_distances_s, grouping=popmap, groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_s, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap, rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

####### 2.6 Grp1 = MidAtl Grp2 = SYah

# Import geographic coordinates and calculate distances

geographic_coordinates <- read.csv("9_loci_Alleles_Mid_Atl_v_S-Yah_with_coord.csv", sep=",", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

ncol <- max(count.fields("./genet_dist_MidAtl_SYah.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./genet_dist_MidAtl_SYah.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="6-9_loci_Alleles_MidAtl_SYah.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

# Create plots and do reg tests using a distances

genet_dist_within_MidAtl_a <- genetic_distances_a[which(popmap == 1), which(popmap == 1)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]

genet_dist_between_a <- genetic_distances_a[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

plot(genetic_distances_a ~ geographic_distances, col="red", pch=19, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), ylab="Genetic distance", xlab="log(Geographical distance)")
abline(lm(c(genetic_distances_a) ~ c(geographic_distances)), col="red")
points(genet_dist_within_MidAtl_a ~ geo_dist_within_MidAtl, pch=19)
abline(lm(c(genet_dist_within_MidAtl_a) ~ c(geo_dist_within_MidAtl)))

points(genet_dist_between_a ~ geo_dist_between, pch="+", cex=1.5, col="blue")
abline(lm(c(genet_dist_between_a) ~ c(geo_dist_between)), col="blue")

H_reg <- diffreg(geo_dist_within_MidAtl, geographic_distances, genet_dist_within_MidAtl_a, genetic_distances_a)

# Create plots and do reg tests using shared allele distances

genet_dist_within_MidAtl_s <- genetic_distances_s[which(popmap == 1), which(popmap == 1)]
geo_dist_within_MidAtl <- geographic_distances[which(popmap == 1), which(popmap == 1)]

genet_dist_between_s <- genetic_distances_s[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]

plot(genetic_distances_s ~ geographic_distances, col="red", pch=19, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), ylab="Genetic distance", xlab="log(Geographical distance)")
abline(lm(c(genetic_distances_s) ~ c(geographic_distances)), col="red", lty=3, lwd=3)
points(genet_dist_within_MidAtl_s ~ geo_dist_within_MidAtl, pch=19)
abline(lm(c(genet_dist_within_MidAtl_s) ~ c(geo_dist_within_MidAtl)), lty=3, lwd=3)

points(genet_dist_between_s ~ geo_dist_between, pch="+", cex=1.5, col="blue")
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue")

H_reg <- diffreg(geo_dist_within_MidAtl, geographic_distances, genet_dist_within_MidAtl_s, genetic_distances_s)


##### 3- Tests to detect barriers to gene flow for coastal lineages (for more details about the groups, cf. Table 2 & 3 of the manuscript)
##### Matrices of geographic and genetic distances are available in the same github repository

geographic_coordinates <- read.csv("7-Coordinates_All_Atlantic_groups.csv", sep=",", header=F)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4:6]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Calculate genetic distances 2: shared allele distances calculated in prabclus

genepop_file <- alleleconvert(file="7-9_loci_Alleles_All_Atlantic_groups.gen", format.in = "genepop", firstcolname = T)
shared_allele_dist <- alleleinit(allelematrix=genepop_file, rows.are.individuals = T)
genetic_distances_s <- shared_allele_dist$distmat

ncol <- max(count.fields("./7-Genet_dist_a_All_Atlantic_groups.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./7-Genet_dist_a_All_Atlantic_groups.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples

genet_dist_within_AtlS <- genetic_distances[which(popmap[,2] == 1), which(popmap[,2] == 1)]
geo_dist_within_AtlS <- geographic_distances[which(popmap[,2] == 1), which(popmap[,2] == 1)]

genet_dist_between <- genetic_distances[which(popmap == 1), which(popmap == 2)]
geo_dist_between <- geographic_distances[which(popmap == 1), which(popmap == 2)]


plot(genet_dist_within_AtlS ~ geo_dist_within_AtlS, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances), max(genetic_distances)), pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
abline(lm(c(genet_dist_within_AtlS) ~ c(geo_dist_within_AtlS)))
points(genet_dist_between ~ geo_dist_between, pch=3, col="green")
abline(lm(c(genet_dist_between) ~ c(geo_dist_between)), col="green")

## Testing N-Atl versus all other Atlantic groups (cf. Table 2)

# Create plots and do reg tests using shared allele distances

genet_dist_within_NAtl_s <- genetic_distances_s[which(popmap[,3] == 2), which(popmap[,3] == 2)]
geo_dist_within_NAtl <- geographic_distances[which(popmap[,3] == 2), which(popmap[,3] == 2)]
genet_dist_within_SAtl_s <- genetic_distances_s[which(popmap[,3] == 1), which(popmap[,3] == 1)]
geo_dist_within_SAtl <- geographic_distances[which(popmap[,3] == 1), which(popmap[,3] == 1)]

genet_dist_between_s <- genetic_distances_s[which(popmap[,3] == 1), which(popmap[,3] == 2)]
geo_dist_between <- geographic_distances[which(popmap[,3] == 1), which(popmap[,3] == 2)]

plot(genet_dist_within_NAtl_s ~ geo_dist_within_NAtl, col="red", pch=19, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_s), max(genetic_distances_s)), ylab="Genetic distance", xlab="log(Geographical distance)")
abline(lm(c(genet_dist_within_NAtl_s) ~ c(geo_dist_within_NAtl)), col="red", lty=3, lwd=3)
points(genet_dist_within_SAtl_s ~ geo_dist_within_SAtl, pch=19)
abline(lm(c(genet_dist_within_SAtl_s) ~ c(geo_dist_within_SAtl)), lty=3, lwd=3)

points(genet_dist_between_s ~ geo_dist_between, pch="+", cex=1.5, col="blue")
abline(lm(c(genet_dist_between_s) ~ c(geo_dist_between)), col="blue", lwd=2)

H01_reg <- regeqdist(geographic_distances, genetic_distances_s, grouping=popmap[,2], groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_s, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap[,2], rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_s, grouping=popmap[,2], rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

# Create plots and do reg tests using a distances

genet_dist_within_NAtl_a <- genetic_distances_a[which(popmap[,3] == 2), which(popmap[,3] == 2)]
geo_dist_within_NAtl <- geographic_distances[which(popmap[,3] == 2), which(popmap[,3] == 2)]
genet_dist_within_SAtl_a <- genetic_distances_a[which(popmap[,3] == 1), which(popmap[,3] == 1)]
geo_dist_within_SAtl <- geographic_distances[which(popmap[,3] == 1), which(popmap[,3] == 1)]

genet_dist_between_a <- genetic_distances_a[which(popmap[,3] == 1), which(popmap[,3] == 2)]
geo_dist_between <- geographic_distances[which(popmap[,3] == 1), which(popmap[,3] == 2)]

plot(genet_dist_within_NAtl_a ~ geo_dist_within_NAtl, col="red", pch=19, xlim=c(min(geographic_distances), max(geographic_distances)), ylim=c(min(genetic_distances_a), max(genetic_distances_a)), ylab="Genetic distance", xlab="log(Geographical distance)")
abline(lm(c(genet_dist_within_NAtl_a) ~ c(geo_dist_within_NAtl)), col="red")
points(genet_dist_within_SAtl_a ~ geo_dist_within_SAtl, pch=19)
abline(lm(c(genet_dist_within_SAtl_a) ~ c(geo_dist_within_SAtl)))

points(genet_dist_between_a ~ geo_dist_between, pch="+", cex=1.5, col="blue")
abline(lm(c(genet_dist_between_a) ~ c(geo_dist_between)), col="blue")

H01_reg <- regeqdist(geographic_distances, genetic_distances_a, grouping=popmap[,2], groups=c(1,2))
pv <- min(H01_reg$pval)*2
if(pv >= 0.05){ H02_reg <- regdistbetween(geographic_distances, genetic_distances_a, grouping=popmap)
cat("p-value H01 ", pv, "\n", "p-value H02 ", H02_reg$pval, sep="")} else { 
  H03_reg1 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap[,2], rgroup=1)
  H03_reg2 <- regdistbetweenone(geographic_distances, genetic_distances_a, grouping=popmap[,2], rgroup=2)
  cat("p-value H01 ", pv, "\n", "p-value H031 ", H03_reg1$pval, "\n", "p-value H032 ", H03_reg2$pval, sep="")}

## Testing the rest of Atlantic coast (from M-Bou to S-Atl) gradually

# using s genetic distance

genet_dist_no_Natl_s <- genetic_distances_s[which(popmap[,1] != 1), which(popmap[,1] != 1)]
genet_dist_no_Natl_a <- genetic_distances_a[which(popmap[,1] != 1), which(popmap[,1] != 1)]
geo_dist_no_Natl <- geographic_distances[which(popmap[,1] != 1), which(popmap[,1] != 1)]

genet_dist_no_Natl_no_MBOU_s <- genetic_distances_s[which(popmap[,1] %in% c(3,4,5,6)), which(popmap[,1] %in% c(3,4,5,6))]
genet_dist_no_Natl_no_MBOU_a <- genetic_distances_a[which(popmap[,1] %in% c(3,4,5,6)), which(popmap[,1] %in% c(3,4,5,6))]
geo_dist_no_Natl_no_MBOU <- geographic_distances[which(popmap[,1] %in% c(3,4,5,6)), which(popmap[,1] %in% c(3,4,5,6))]

genet_dist_no_Natl_no_MBOU_no_BLADMEHD_s <- genetic_distances_s[which(popmap[,1] %in% c(3,4,5)), which(popmap[,1] %in% c(3,4,5))]
genet_dist_no_Natl_no_MBOU_no_BLADMEHD_a <- genetic_distances_a[which(popmap[,1] %in% c(3,4,5)), which(popmap[,1] %in% c(3,4,5))]
geo_dist_no_Natl_no_MBOU_no_BLADMEHD <- geographic_distances[which(popmap[,1] %in% c(3,4,5)), which(popmap[,1] %in% c(3,4,5))]

genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_s <- genetic_distances_s[which(popmap[,1] %in% c(4,5)), which(popmap[,1] %in% c(4,5))]
genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_a <- genetic_distances_a[which(popmap[,1] %in% c(4,5)), which(popmap[,1] %in% c(4,5))]
geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU <- geographic_distances[which(popmap[,1] %in% c(4,5)), which(popmap[,1] %in% c(4,5))]

genet_dist_only_SAtl_s <- genetic_distances_s[which(popmap[,1] == 5), which(popmap[,1] == 5)]
genet_dist_only_SAtl_a <- genetic_distances_a[which(popmap[,1] == 5), which(popmap[,1] == 5)]
geo_dist_only_SAtl <- geographic_distances[which(popmap[,1] == 5), which(popmap[,1] == 5)]

plot(genet_dist_no_Natl_s ~ geo_dist_no_Natl, col="red", pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
abline(lm(c(genet_dist_no_Natl_s) ~ c(geo_dist_no_Natl)), col="red", lty=3, lwd=3)
points(genet_dist_no_Natl_no_MBOU_s ~ geo_dist_no_Natl_no_MBOU, col="blue", pch=19)
abline(lm(c(genet_dist_no_Natl_no_MBOU_s) ~ c(geo_dist_no_Natl_no_MBOU)), col="blue", lty=3, lwd=3)
points(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_s ~ geo_dist_no_Natl_no_MBOU_no_BLADMEHD, col="orange", pch=19)
abline(lm(c(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_s) ~ c(geo_dist_no_Natl_no_MBOU_no_BLADMEHD)), col="orange", lty=3, lwd=3)
points(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_s ~ geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, col="darkgreen", pch=19)
abline(lm(c(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_s) ~ c(geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU)), col="darkgreen", lty=3, lwd=3)
points(genet_dist_only_SAtl_s ~ geo_dist_only_SAtl, col="black", pch=19)
abline(lm(c(genet_dist_only_SAtl_s) ~ c(geo_dist_only_SAtl)), col="black", lty=3, lwd=3)

# Regression lines comparisons. The order of the tests refers to Table 3, coastal lineages section
H_reg1 <- diffreg(geo_dist_no_Natl, geo_dist_no_Natl_no_MBOU, genet_dist_no_Natl_s, genet_dist_no_Natl_no_MBOU_s)
H_reg2 <- diffreg(geo_dist_no_Natl, geo_dist_no_Natl_no_MBOU_no_BLADMEHD, genet_dist_no_Natl_s, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_s)
H_reg3 <- diffreg(geo_dist_no_Natl, geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, genet_dist_no_Natl_s, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_s)
H_reg4 <- diffreg(geo_dist_no_Natl, geo_dist_only_SAtl, genet_dist_no_Natl_s, genet_dist_only_SAtl_s)

H_reg5 <- diffreg(geo_dist_no_Natl_no_MBOU, geo_dist_no_Natl_no_MBOU_no_BLADMEHD, genet_dist_no_Natl_no_MBOU_s, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_s)
H_reg6 <- diffreg(geo_dist_no_Natl_no_MBOU, geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, genet_dist_no_Natl_no_MBOU_s, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_s)
H_reg7 <- diffreg(geo_dist_no_Natl_no_MBOU, geo_dist_only_SAtl, genet_dist_no_Natl_no_MBOU_s, genet_dist_only_SAtl_s)

H_reg8 <- diffreg(geo_dist_no_Natl_no_MBOU_no_BLADMEHD, geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_s, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_s)
H_reg9 <- diffreg(geo_dist_no_Natl_no_MBOU_no_BLADMEHD, geo_dist_only_SAtl, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_s, genet_dist_only_SAtl_s)

H_reg10 <- diffreg(geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, geo_dist_only_SAtl, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_s, genet_dist_only_SAtl_s)

# using a genetic distance

plot(genet_dist_no_Natl_a ~ geo_dist_no_Natl, col="red", pch=19, ylab="Genetic distance", xlab="log(Geographical distance)")
abline(lm(c(genet_dist_no_Natl_a) ~ c(geo_dist_no_Natl)), col="red")
points(genet_dist_no_Natl_no_MBOU_a ~ geo_dist_no_Natl_no_MBOU, col="blue", pch=19)
abline(lm(c(genet_dist_no_Natl_no_MBOU_a) ~ c(geo_dist_no_Natl_no_MBOU)), col="blue")
points(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_a ~ geo_dist_no_Natl_no_MBOU_no_BLADMEHD, col="orange", pch=19)
abline(lm(c(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_a) ~ c(geo_dist_no_Natl_no_MBOU_no_BLADMEHD)), col="orange")
points(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_a ~ geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, col="darkgreen", pch=19)
abline(lm(c(genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_a) ~ c(geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU)), col="darkgreen")
points(genet_dist_only_SAtl_a ~ geo_dist_only_SAtl, col="black", pch=19)
abline(lm(c(genet_dist_only_SAtl_a) ~ c(geo_dist_only_SAtl)), col="black")

# Regression lines comparisons. The order of the tests refers to Table 3, coastal lineages section
H_reg1 <- diffreg(geo_dist_no_Natl, geo_dist_no_Natl_no_MBOU, genet_dist_no_Natl_a, genet_dist_no_Natl_no_MBOU_a)
H_reg2 <- diffreg(geo_dist_no_Natl, geo_dist_no_Natl_no_MBOU_no_BLADMEHD, genet_dist_no_Natl_a, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_a)
H_reg3 <- diffreg(geo_dist_no_Natl, geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, genet_dist_no_Natl_a, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_a)
H_reg4 <- diffreg(geo_dist_no_Natl, geo_dist_only_SAtl, genet_dist_no_Natl_a, genet_dist_only_SAtl_a)

H_reg5 <- diffreg(geo_dist_no_Natl_no_MBOU, geo_dist_no_Natl_no_MBOU_no_BLADMEHD, genet_dist_no_Natl_no_MBOU_a, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_a)
H_reg6 <- diffreg(geo_dist_no_Natl_no_MBOU, geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, genet_dist_no_Natl_no_MBOU_a, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_a)
H_reg7 <- diffreg(geo_dist_no_Natl_no_MBOU, geo_dist_only_SAtl, genet_dist_no_Natl_no_MBOU_a, genet_dist_only_SAtl_a)

H_reg8 <- diffreg(geo_dist_no_Natl_no_MBOU_no_BLADMEHD, geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_a, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_a)
H_reg9 <- diffreg(geo_dist_no_Natl_no_MBOU_no_BLADMEHD, geo_dist_only_SAtl, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_a, genet_dist_only_SAtl_a)

H_reg10 <- diffreg(geo_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU, geo_dist_only_SAtl, genet_dist_no_Natl_no_MBOU_no_BLADMEHD_no_SBOU_a, genet_dist_only_SAtl_a)
