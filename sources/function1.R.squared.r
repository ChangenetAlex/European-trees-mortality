rsquared.ACInfo = "rsquared.AC \n A R function to calculate various components  \n Require spaMM package \n \n Input data: Your model (HLfit) \n
Return various components:
The family of your model
The link function
The method used to calculate r2 (cf Nakagawa)
The margina and conditionnal R2
The loglikelihood and the AIC
The IntraClassCorrelation coef
The Proportion change in variance "
bannerBreak = "\n*********************************************************************\n"
cat(paste0(bannerBreak, rsquared.ACInfo, bannerBreak, "\n"))

library(spaMM)
rsquared.AC <- function(model, method = "trigamma") {
  if (is.null(method))
    method <- "trigamma"
  
  link <- model$family$link
  
  family. <- model$family$family[1]
  
  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)
  
  X <- model[["X.pv"]] # For spaMM #xtraire la cor matrix
  
  sigmaF <-
    var(as.vector(spaMM::fixef(model) %*% t(X))) #lambda Full
  
  #sigma <- unclass(model[["lambda"]]) # Pourquoi pas le random model ? take variance random effect lambd random = alpha = sigmaL
  sigma <-
    sum(as.numeric(unclass(model[["lambda"]]))) #ptet en faire la somme sigma random effect
  
  data <- model$data
  
  
  if (family. == "poisson") {
    if (!method %in% c("delta", "lognormal", "trigamma"))
      stop("Unsupported method!")
    
    #sigmaL <- sum(GetOLRE(sigma, modele, X, data, RE = "RE")) ## ? residual error... Juste variance sigma
    #sigmaE <- sum(GetOLRE(sigma, modele, X, data, RE = "OLRE")) ## ? OLRE additive dispersion compoenent (part of the RE) #pbm with poisson
    # Serta rien on a déjà lambda. C'est lambda e
    
    rand <-
      onlyBars(model[["predictor"]]) # extraire just les random effects
    f <-
      paste(model[["predictor"]][[2]], " ~ 1 + ", onlyBars(model[["predictor"]], slopes = FALSE)) # Model null
    if (length(grep("Matern", model[["predictor"]], fixed = T)) >= 1) {
      f <- sub("(1 | lat", "Matern(1 | lat", f, fixed = T)
    } else {
      f <- f
    }
    nullmodel <-
      suppressWarnings(spaMM::fitme(
        formula(f),
        family = poisson(link = link),
        data = data
      ))
    as.character(rand)
    #lambda <- attr(VarCorr(nullmodele), "sc")^2
    lambda <-
      as.numeric(exp(fixef(nullmodel) + 0.5 * sum(as.numeric(
        unclass(nullmodel[["lambda"]])
      ))))
    #lambda <- exp(fixef(nullmodele)[1] + (sigmaL + sigmaE)/2)
    
    omega <- 1 #disp full model
    omegaN <- 1 #disp null model
    
    #Value of omega to extract when overdisp > 1
    
    if (link == "mu^0.5") {
      sigmaD <- 0.25 * omega
      sigmaDN <- 0.25 * omegaN
    }
    else {
      #square root link
      
      if (link == "log") {
        nu <- omega / lambda
        nuN <- omegaN / lambda
        
        if (method == "delta") {
          sigmaD <- nu
          sigmaDN <- nuN
        }
        
        if (method == "lognormal") {
          sigmaD <- log(1 + nu)
          sigmaDN <- log(1 + nuN)
        }
        
        if (method == "trigamma") {
          sigmaD <- trigamma(1 / nu)
          sigmaDN <- trigamma(1 / nuN)
        }
        
      } else
        stop("Unsupported link function!")
      
    }
    
  } else
    
    if (family. == "binomial") {
      if (method == "trigamma")
        method <- "delta"
      
      if (!method %in% c("theoretical", "delta"))
        stop("Unsupported method!")
      
      #sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "all")) #Variance of random effect (different de nakagawa car based on the real model
      
      # Rajout d'un calcul de sigma L ou garder sigma tuot court c'est bien egale a sigma tout court
      
      sigmaE <- 0 #Additive disp
      
      if (link == "logit") {
        if (method == "theoretical")
          sigmaD <- pi^2 / 3 #okay
        
        if (method == "delta") {
          sigmaDN <- pi^2 / 3
          
          rand <-
            onlyBars(model[["predictor"]]) # extraire just les random effects
          f <-
            paste(model[["predictor"]][[2]], " ~ 1 + ", onlyBars(model[["predictor"]], slopes = FALSE)) # Model null
          if (length(grep("Matern", model[["predictor"]], fixed = T)) >=
              1) {
            f <- sub("(1 | lat", "Matern(1 | lat", f, fixed = T)
          } else {
            f <- f
          }
          nullmodel <-
            suppressWarnings(spaMM::fitme(
              formula(f),
              family = binomial(link = link),
              data = data
            ))
          
          vt <-
            sum(as.numeric(nullmodel[["lambda"]])) # somme des effets aléatoire du full
          
          # check a la fin le sigmaL et si o peut l'obtenir
          pmean <- plogis(as.numeric(fixef(nullmodel)) - 0.5 * vt *
                            tanh(as.numeric(fixef(nullmodel)) * (1 + 2 * exp(-0.5 * vt)) /
                                   6))
          
          sigmaD <- 1 / (pmean * (1 - pmean)) #result presque idem
          
        }
        
      } else if (link == "cloglog") {
        if (method == "theoretical")
          sigmaD <- pi^2 / 6 #okay
        
        if (method == "delta") {
          sigmaDN <- pi^2 / 6
          
          rand <-
            onlyBars(model[["predictor"]]) # extraire just les random effects
          f <-
            paste(model[["predictor"]][[2]], " ~ 1 + ", onlyBars(model[["predictor"]], slopes = FALSE)) # Model null
          if (length(grep("Matern", model[["predictor"]], fixed = T)) >=
              1) {
            f <- sub("(1 | lat", "Matern(1 | lat", f, fixed = T)
          } else {
            f <- f
          }
          nullmodel <-
            suppressWarnings(spaMM::fitme(
              formula(f),
              family = binomial(link = link),
              data = data
            ))
          
          vt <-
            sum(as.numeric(nullmodel[["lambda"]])) # somme des effets aléatoire du full
          
          # check a la fin le sigmaL et si o peut l'obtenir
          pmean <- plogis(as.numeric(fixef(nullmodel)) - 0.5 * vt *
                            tanh(as.numeric(fixef(nullmodel)) * (1 + 2 * exp(-0.5 * vt)) /
                                   6))
          
          sigmaD <-
            pmean / ((log(1 - pmean)) ^ 2 * (1 - pmean)) #result presque idem
          
        }
        
        
      }
      
      
    } else
      
      if (family. == "Gamma") {
        if (!method %in% c("delta", "lognormal", "trigamma"))
          stop("Unsupported method!")
        
        rand <-
          onlyBars(model[["predictor"]]) # extraire just les random effects
        f <-
          paste(model[["predictor"]][[2]], " ~ 1 + ", onlyBars(model[["predictor"]], slopes = FALSE)) # Model null
        if (length(grep("Matern", model[["predictor"]], fixed = T)) >= 1) {
          f <- sub("(1 | lat", "Matern(1 | lat", f, fixed = T)
        } else {
          f <- f
        }
        nullmodel <-
          suppressWarnings(spaMM::fitme(
            formula(f),
            family = binomial(link = link),
            data = data
          ))
        
        #sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "all")) #effet random = sigma normal
        
        sigmaE <- 0
        
        #lambda <- attr(VarCorr(model), "sc")^2 Method to get variance residuals avec lme4.
        lambda <- model$phi #check that is the real value
        lambdaN <- nullmodel$phi #check that is the real value
        
        omega <- 1
        
        if (link == "log") {
          nu <- omega / lambda
          nuN <- omega / lambdaN
          
          if (method == "delta") {
            sigmaD <- 1 / nu
            sigmaDN <- 1 / nuN
          }
          if (method == "lognormal") {
            sigmaD <- log(1 + 1 / nu)
            sigmaDN <- log(1 + 1 / nuN)
          }
          if (method == "trigamma") {
            sigmaD <- trigamma(nu)
            sigmaDN <- trigamma(nuN)
          }
          
        } else
          
          if (link == "inverse") {
            
          } else
            stop("Unsupported link function!")
        
      } else
        
        if (family. == "negbin") {
          if (!method %in% c("delta", "lognormal", "trigamma"))
            stop("Unsupported method!")
          
          
          rand <-
            onlyBars(model[["predictor"]]) # extraire just les random effects
          f <-
            paste(model[["predictor"]][[2]], " ~ 1 + ", onlyBars(model[["predictor"]], slopes = FALSE)) # Model null
          if (length(grep("Matern", model[["predictor"]], fixed = T)) >=
              1) {
            f <- sub("(1 | lat", "Matern(1 | lat", f, fixed = T)
          } else {
            f <- f
          }
          nullmodel <-
            suppressWarnings(spaMM::fitme(
              formula(f),
              family = negbin(link = link),
              data = data
            ))
          
          
          #lambda <- attr(VarCorr(nullmodele), "sc")^2
          lambda <-
            as.numeric(exp(fixef(nullmodel) + 0.5 * sum(as.numeric(
              unclass(nullmodel[["lambda"]])
            ))))
          theta <-
            environment(model[["family"]][["linkfun"]])[["shape"]]
          thetaN <-
            environment(nullmodel[["family"]][["linkfun"]])[["shape"]]
          
          if (link == "log") {
            nu <- (1 / lambda) + (1 / theta)
            nuN <- (1 / lambda) + (1 / thetaN)
            
            if (method == "delta") {
              sigmaD <- nu
              sigmaDN <- nuN
            }
            
            if (method == "lognormal") {
              sigmaD <- log(1 + nu)
              sigmaDN <- log(1 + nuN)
            }
            
            if (method == "trigamma") {
              sigmaD <- trigamma(nu ^ (-1))
              sigmaDN <- trigamma(nuN ^ (-1))
            }
            
          } else
            stop("Unsupported link function!")
          
        } else
          stop("Unsupported family!")
  
  if (length(grep("Matern", model[["predictor"]], fixed = T)) >= 1) {
    mar <- (sigmaF) / (sigmaF + sigmaD + sigma) #modified
    con <- (sigmaF + sigma) / (sigmaF + sigma + sigmaD) #modified
    marAIC <- AIC(model)[[1]]
    hLikelihood <- model$APHLs$p_v
    ICCraw1 <-
      as.numeric(unclass(nullmodel[["lambda"]][[1]])) / sum(as.numeric(unclass(nullmodel[["lambda"]]))) +
      sigmaDN # pour binomial = ICC theorique
    ICCadj1 <-
      as.numeric(unclass(model[["lambda"]][[1]])) / sum(as.numeric(unclass(model[["lambda"]]))) +
      sigmaD #pourbinomial = ICC obs
    ICCraw2 <-
      as.numeric(unclass(nullmodel[["lambda"]][[2]])) / sum(as.numeric(unclass(nullmodel[["lambda"]]))) +
      sigmaDN #idem
    ICCadj2 <-
      as.numeric(unclass(model[["lambda"]][[2]])) / sum(as.numeric(unclass(model[["lambda"]]))) +
      sigmaD #idem
    PCV1 <-
      1 - as.numeric(unclass(model[["lambda"]][[1]])) / as.numeric(unclass(nullmodel[["lambda"]][[1]]))
    PCV2 <-
      1 - as.numeric(unclass(model[["lambda"]][[2]])) / as.numeric(unclass(nullmodel[["lambda"]][[2]]))
    PCVran <-
      1 - sum(as.numeric(unclass(model[["lambda"]]))) / sum(as.numeric(unclass(nullmodel[["lambda"]])))
    PCVObs <- 1 - sigmaD / sigmaDN # Si = 0 => poisson ou gamma
    
    list(
      family = family.,
      link = link,
      method = method,
      Marginal = mar,
      Conditional = con,
      Lik = hLikelihood,
      AICm = marAIC,
      ICCraw1 = ICCraw1,
      ICCadj1 = ICCadj1,
      PCVran = PCVran,
      PCVObs = PCVObs,
      ICCraw2 = ICCraw2,
      ICCadj2 = ICCadj2,
      PCV1 = PCV1,
      PCV2 = PCV2
    )
  } else {
    mar <- (sigmaF) / (sigmaF + sigmaD + sigma) #modified
    con <- (sigmaF + sigma) / (sigmaF + sigma + sigmaD) #modified
    marAIC <- AIC(model)[[1]]
    hLikelihood <- model$APHLs$p_v
    ICCraw1 <-
      as.numeric(unclass(nullmodel[["lambda"]][[1]])) / sum(as.numeric(unclass(nullmodel[["lambda"]]))) +
      sigmaDN # pour binomial = ICC theorique
    ICCadj1 <-
      as.numeric(unclass(model[["lambda"]][[1]])) / sum(as.numeric(unclass(model[["lambda"]]))) +
      sigmaD #pourbinomial = ICC obs
    PCVran <-
      1 - sum(as.numeric(unclass(model[["lambda"]]))) / sum(as.numeric(unclass(nullmodel[["lambda"]])))
    PCVObs <- 1 - sigmaD / sigmaDN # Si = 0 => poisson ou gamma
    
    list(
      family = family.,
      link = link,
      method = method,
      Marginal = mar,
      Conditional = con,
      Lik = hLikelihood,
      AICm = marAIC,
      ICCraw1 = ICCraw1,
      ICCadj1 = ICCadj1,
      PCVran = PCVran,
      PCVObs = PCVObs,
      ICCraw2 = NA,
      ICCadj2 = NA,
      PCV1 = NA,
      PCV2 = NA
    )
  }
}
