
# helper functions

Nm <- function(Fst) {
  Nm <- (1/4)*((1/Fst) - 1)
  return(Nm)
}

calculate_snmf_entropy <- function(snmf_res=NULL, maxk=NULL) {
  for (i in 1:maxk) assign(paste0("ce_k", i), cross.entropy(snmf_res, K=i))
  tmp <- ls()[grep(ls(), pattern="ce_k")]
  for (i in 1:maxk) {
    ce_ki <- as.data.frame(get(tmp[i]))
    ki <- colnames(ce_ki)
    ki <- substr(ki,5,nchar(ki))
    ce_ki$V2 <- ki
    colnames(ce_ki) <- c("CE", "K")
    assign(paste0("ce_k", ki), ce_ki)
    rm(ce_ki)
  }
  ce <- data.frame()
  for (i in 1:maxk) ce <- rbind.data.frame(ce, get(tmp[i]))
  ce$K <- as.numeric(ce$K)
  return(ce)
}

snmf_barplot <- function(snmf_res=NULL, K=NULL, dapc_assign=NULL, pop=NULL) {
  best.run <- as.data.frame(cross.entropy(object = snmf_res, K = K))
  best.run <- which(best.run[, 1] == min(best.run[, 1]))
  print(paste0("best run: ", best.run))
  dapc_assign <- as.character(dapc_assign)
  K_snmf <- as.data.frame(Q(snmf_res, K = K, run = best.run))
  K_snmf$Grp <- dapc_assign
  K_snmf$Pop <- pop
  K_snmf <- K_snmf[order(K_snmf$Grp), ]

  dapc_assign <- as.data.frame(K_snmf$Grp)
  colnames(dapc_assign) <- "Cl."
  pop <- as.data.frame(K_snmf$Pop)
  colnames(pop) <- "Pop."
  tab_label <- cbind.data.frame(pop, dapc_assign)

  K_snmf <- K_snmf[, -c(ncol(K_snmf)-1,ncol(K_snmf))]
  qlist <- list(K_snmf)
  names(qlist) <- paste0("K=", K)
  qlist <- as.qlist(qlist)
  qlist <- alignK(qlist)
  snmf_plot <- plotQ(qlist, basesize = 11,
                     # clustercol = c("blue", "red", "lig"),
                     sortind = "Cluster2", sharedindlab = FALSE, showsubtitle = T,
                     subtitlelab = "Global ancestry coefficients", showlegend = F,
                     showsp = T, splab = names(qlist), grplab = tab_label,
                     grplabsize = 4, linesize = 0.5, pointsize = 4,
                     returnplot = T, exportplot = F,
                     showyaxis = T, panelspacer = 0.25, panelratio = c(3, 1))

  grid.arrange(snmf_plot$plot[[1]])
}

lfmm_qvalue <- function(pv = NULL, bim = NULL) {

  res <- as.data.frame(matrix(rep(NA, nrow(pv)*ncol(pv)), nrow(pv)))
  colnames(res) <- colnames(pv)
  rownames(res) <- rownames(pv)

  for (i in 1:ncol(pv)) {
    q <- qvalue::qvalue(pv[,i])
    q <- q$qvalues
    res[,i] <- q
    rm(q)
  };rm(i)

  if (!is.null(bim)) {
    bim <- data.table::fread(bim)
    if (length(which(rownames(pv) == bim$V2)) == nrow(pv)) {
      res$SNP <- bim$V2
      res$CHR <- bim$V1
      res$BP <- bim$V4
    } else {
      cat("Error: markers are not in the same order in 'pv' and 'bim'.")
      cat("Information cannot be matched properly: please check order of markers prior to run 'lfmm_qvalue'.")
    }
  }

  return(res)

}

lfmm_qvalue_cut <- function (qvalues = NULL, nenv = NULL, cutoff = NULL, bim.info = TRUE, use.bim = FALSE, bim = NULL) {

  if (bim.info == TRUE) {
    res <- list()
    e <- qvalues[, 1:nenv]
    s <- qvalues[, (nenv+1):(nenv+3)]
    for (i in 1:ncol(e)) {
      j <- which(e[, i] < cutoff)
      if (length(j) == 0) {
        t <- data.frame("ENV" = colnames(e)[i], "SNP" = NA, "CHR" = NA, "BP" = NA, "qval" = NA)
        res[[i]] <- t
      } else {
        t <- data.frame("ENV" = colnames(e)[i], "SNP" = s[j, 1], "CHR" = s[j, 2], "BP" = s[j, 3], "qval" = e[j, i])
        res[[i]] <- t
      }
    }
    res <- do.call(rbind.data.frame, res)
    res <- na.omit(res)
    return(res)
  }

  if (bim.info == FALSE & use.bim == FALSE) {
    res <- list()
    for (i in 1:ncol(qvalues)) {
      j <- which(qvalues[, i] < cutoff)
      if (length(j) == 0) {
        t <- data.frame("ENV" = colnames(qvalues)[i], "SNP" = NA, "qval" = NA)
        res[[i]] <- t
      } else {
        t <- data.frame("ENV" = colnames(qvalues)[i], "SNP" = rownames(qvalues)[j], "qval" = qvalues[j, i])
        res[[i]] <- t
      }
    }
    res <- do.call(rbind.data.frame, res)
    res <- na.omit(res)
    return(res)
  }

  if (bim.info == FALSE & use.bim == TRUE) {

    cat("Attention: markers in 'qvalues' and 'bim' are assumed to be in the same order.")
    cat("Please check this condition before using 'lfmm_qvalue_cut' with 'bim.info = FALSE' and 'use.bim = TRUE'.")

    bim <- data.table::fread(bim)
    qvalues$SNP <- bim$V2
    qvalues$CHR <- bim$V1
    qvalues$BP <- bim$V4

    res <- list()
    e <- qvalues[, 1:nenv]
    s <- qvalues[, (nenv+1):(nenv+3)]
    for (i in 1:ncol(e)) {
      j <- which(e[, i] < cutoff)
      if (length(j) == 0) {
        t <- data.frame("ENV" = colnames(e)[i], "SNP" = NA, "CHR" = NA, "BP" = NA, "qval" = NA)
        res[[i]] <- t
      } else {
        t <- data.frame("ENV" = colnames(e)[i], "SNP" = s[j, 1], "CHR" = s[j, 2], "BP" = s[j, 3], "qval" = e[j, i])
        res[[i]] <- t
      }
    }
    res <- do.call(rbind.data.frame, res)
    res <- na.omit(res)
    return(res)
  }
}

rda_outliers <- function(load_rda=NULL, n_rda=NULL, n_sign=NULL, z=NULL){
  lims <- mean(load_rda[, n_rda]) + c(-1, 1) * z * sd(load_rda[, n_rda])
  tmp <- which(load_rda[, n_rda] < lims[1] | load_rda[, n_rda] > lims[2])
  tmp <- load_rda[tmp, c(n_rda, (n_sign+1):ncol(load_rda))]
}

rda_outliers_env <- function(cand1=NULL, cand2=NULL, cand3=NULL, env=NULL, Y=NULL) {
  cand1 <- cbind.data.frame(rep(1,times=nrow(cand1)), cand1)
  cand2 <- cbind.data.frame(rep(2,times=nrow(cand2)), cand2)
  cand3 <- cbind.data.frame(rep(3,times=nrow(cand3)), cand3)
  colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","loading", "SNP", "CHR", "BP")
  cand <- rbind.data.frame(cand1, cand2, cand3)
  cand$SNP <- as.character(cand$SNP)
  foo <- matrix(nrow=nrow(cand), ncol=ncol(env))
  colnames(foo) <- colnames(env)
  for (i in 1:length(cand$SNP)) {
    nam <- cand$SNP[i]
    snp.gen <- Y[,nam]
    foo[i,] <- apply(env,2,function(x) cor(x,snp.gen))
  }
  foo <- abs(foo)
  foo <- as.data.frame(foo)
  foo$ENV <- NA
  for(i in 1:length(cand$SNP)) {
    t <- which(foo[i, 1:(ncol(foo)-1)] == max(foo[i, 1:(ncol(foo)-1)]))
    foo$ENV[i] <- colnames(foo)[t]
  }
  cand <- cbind.data.frame(cand,foo)
  rownames(cand)<-NULL
  cand <- cand[, c("ENV", "SNP", "CHR", "BP", "axis", "loading")]
  return(cand)
}

# Credits to Capblancq, T., & Forester, B. R. (2021). Redundancy analysis: A Swiss Army Knife for landscape genomics. Methods in Ecology and Evolution, 12(12), 2298-2309.
adaptive_index <- function(RDA, K, env_pres, range = NULL, method = "loadings", scale_env, center_env){

  # Formatting environmental rasters for projection
  var_env_proj_pres <- as.data.frame(rasterToPoints(env_pres[[row.names(RDA$CCA$biplot)]]))

  # Standardization of the environmental variables
  var_env_proj_RDA <- as.data.frame(scale(var_env_proj_pres[,-c(1,2)], center_env[row.names(RDA$CCA$biplot)], scale_env[row.names(RDA$CCA$biplot)]))

  # Predicting pixels genetic component based on RDA axes
  Proj_pres <- list()
  if(method == "loadings"){
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(apply(var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], 1, function(x) sum( x * RDA$CCA$biplot[,i])))), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }

  # Prediction with RDA model and linear combinations
  if(method == "predict"){
    pred <- predict(RDA, var_env_proj_RDA[,names(RDA$CCA$biplot[,i])], type = "lc")
    for(i in 1:K){
      ras_pres <- rasterFromXYZ(data.frame(var_env_proj_pres[,c(1,2)], Z = as.vector(pred[,i])), crs = crs(env_pres))
      names(ras_pres) <- paste0("RDA_pres_", as.character(i))
      Proj_pres[[i]] <- ras_pres
      names(Proj_pres)[i] <- paste0("RDA", as.character(i))
    }
  }

  # Mask with the range if supplied
  if(!is.null(range)){
    Proj_pres <- lapply(Proj_pres, function(x) mask(x, range))
  }

  # Returning projections for current climates for each RDA axis
  return(Proj_pres = Proj_pres)
}
