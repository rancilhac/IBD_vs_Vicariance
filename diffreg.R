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
