
if (!all(c("sommer", "lme4")%in% installed.packages())) {
  if(!"sommer" %in% installed.packages()) {install.packages("sommer")
    library(sommer)}
  if(!"lme4" %in% installed.packages()) {install.packages("lme4")
    library(lme4)}
} else {library(sommer); library(lme4)}

#### Analise de modelos mistos para ensaio com delineamento de blocos aumentados
#### e completos - lme4 package

analyzeTrial.lme4FD <- function(x){
    modfit <- lmer(Y ~ control + block_number + (1 | accession_name:new), data=x, REML = T)
    return(modfit)
}

analyzeTrial.lme4 <- function(x){
  if(any(x$studyDesign == "DBA")){
    modfit <- lmer(y ~ check + rep + (1 | clone:new), data=x, REML = T)
    return(modfit)
  } else {
    modfit <- lmer(y ~ rep + (1 | clone), data=x, REML = T)
    return(modfit)
  }
}

analyzeTrialrdMod.lme4 <- function(x, trait){
  modfit <- lm(Y ~ control + block_number, data=x)
  return(modfit)
}

Deviance.MM <- function(x, y) {
  anova(x, y)
}

getVarComp.lme4 <- function(model){
  as.data.frame(VarCorr(model))[,c("grp","vcov")] %>%
    rename(c("grp" = "effect", "vcov" = "VarEstimate")) %>%
    spread(key = effect, value = VarEstimate, fill=NA)
}


getVarOutput.lme4 <- function(modfit){
  var_label <- row.names(modfit)
  if (any(is.na(modfit))){
    varcomp <- data.frame(NA, NA)
  } else{
    varcomp <- data.frame(matrix(unlist(modfit), nrow=nrow(modfit), byrow=T))
  }
  colnames(varcomp) <- c("gen_var", "res_var")
  varcomp <- round(varcomp[,c("gen_var","res_var")], 3)
  var.out <- cbind(var_label, varcomp)
  return(var.out)
}



getBLUPs.lme4 <- function(model){
  as.data.frame(ranef(model)$"accession_name:new" +
                  fixef(model)["(Intercept)"]) %>%
    rename(c("(Intercept)" = names(model)))
     # spread(key = Clone, value = BLUP, fill=NA)
}


##### #### Analise de modelos mistos para ensaio com delineamento de blocos aumentados
#### e completos - sommer package
analyzeTrial.sommer <- function(x){
  if(any(x$studyDesign == "CET")){
    #### Blocos Aumentados
    modfit <- mmer(y ~ check + rep,
                   random = ~ clone:new,
                   data = x,
                   tolparinv = 1e-6,
                   verbose = F,
                   method = "EM",
                   getPEV = T)
  } else {
    #### Blocos Completos
    modfit <- mmer(y ~ rep,
                   random = ~ clone,
                   data=x,
                   tolparinv = 1e-6,
                   verbose = F,
                   method = "EM",
                   getPEV = T)
  }
  return(modfit)
}


getVarComp.sommer <- function(model){
  unlist(model$sigma) %>%
    data.frame(effect = c("Clone", "Residual"), VarEstimate = .) %>%
    spread(key = effect, value = VarEstimate, fill=NA)
}


analyzeTrial.sommerConj <- function(x){
  #### Blocos Completos
  if(length(unique(x$trial)) > 1){
    modfit <- mmer(y ~ 1,
                   random = ~ repTrial + clone + trial:clone,
                   data = x,
                   tolparinv = 1e-5,
                   verbose = F,
                   method = "EM",
                   getPEV = T)
  } else {
    modfit <- mmer(y ~ 1,
                   random = ~ repTrial + clone,
                   data = x,
                   tolparinv = 1e-5,
                   verbose = F,
                   method = "EM",
                   getPEV = T)
  }
  return(modfit)
}

analyzeTrial.lme4Conj <- function(x){
  #### Blocos Completos
  if(length(unique(x$LocYear)) > 1){
    modfit <- lmer(y ~ (1|repTrial) + (1|clone) + (1|LocYear:clone),
                   data = x,
                   REML = T)
  } else {
    modfit <- lmer(y ~ (1|repTrial) + (1|clone),
                   data = x,
                   REML = T)
  }
  return(modfit)
}
