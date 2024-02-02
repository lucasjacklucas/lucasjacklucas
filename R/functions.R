latentvar <- function(data.name){
  # load packages
  library(R2jags)
  library(tidyverse)
  data <- data.name
  # correlation test function
  ctest.function <- function(x){
    test <- length(table(sign(cor(x, use="pair"))>0))>1
    if (test==TRUE) {
      print("Some negative correlations among items. Using Model 2.")
    } else {
      print("All correlations positive. Using Model 1")
    }}

  # remove if all missing
  cols <- ncol(data)
  data$test <- apply(data[,2:cols],MARGIN=1,function(x) sum(is.na(x)))
  data <- data %>% filter(test!=cols-1) %>% dplyr::select(-test)
  items <- data %>% dplyr::select(-1)

  # set up data for latent model
  matrix <- items %>%
    mutate(across(everything(),function(x)c(scale(x)))) %>%
    as.matrix()
  list.data <- list(
    y = matrix,
    nrow = nrow(matrix),
    K = ncol(matrix)
  )
  K <- ncol(matrix)
  ctest.function(items)

  #model1
  mod1 <- function(){
    for (i in 1:nrow){
      latent[i] ~ dnorm(0,1)
      for(k in 1:K){
        y[i,k] ~ dnorm(mu[i,k], tau[k])
        mu[i,k] <- b[k]*latent[i]
      }
    }
    for(i in 1:K){
      b[i] ~ dnorm(0,.1);T(0, )
      tau[i]~dgamma(1, .1)
    }
  }

  #model2
  mod2 <- function(){
    for (i in 1:nrow){
      latent[i] ~ dnorm(0,1)
      for(k in 1:K){
        y[i,k] ~ dnorm(mu[i,k], tau[k])
        mu[i,k] <- b[k]*latent[i]
      }
    }
    b[1] <- 1
    tau[1] ~ dgamma(1,.1)
    for(i in 2:K){
      b[i] ~ dnorm(0,.1)
      tau[i]~dgamma(1, .1)
    }
  }
  model.choice <- ifelse(length(table(sign(cor(items, use="pair"))>0))>1,mod2,mod1)

  #run jags
  latent.params <- c("latent", "b")
  fit.latent <- jags(data=list.data, inits=NULL, latent.params, model.file=model.choice, n.chains=2, n.iter=4000, n.burnin=2000)
  x <- print(fit.latent)
  x.out<-data.frame(x$summary)

  # betas
  x.out$row <- rownames(x.out)
  varnames <- colnames(items)
  betas <- x.out %>% filter(grepl('b', row)) %>%
    rename(lwr = 3, med = 5, upr = 7) %>%
    mutate(var = varnames) %>%
    dplyr::select(var, lwr, med, upr)
  beta.plot <- ggplot(betas, aes(x=reorder(var, med), y=med, ymin=lwr, ymax=upr)) +
    geom_pointrange() + coord_flip() + ylab("Discrimination Parameter") + xlab("")

  # convergence
  convergence <- x.out %>% dplyr::select(row, Rhat, n.eff)

  # extract latent values
  beginning <- K+2
  end <- nrow(x.out)
  latent <- x.out[beginning:end,]
  data$latent <- latent$X50.
  data$latent_se <- latent$sd
  data$latent_lwr <- latent$X2.5.
  data$latent_upr <- latent$X97.5.
  data <- data %>% dplyr::select(1, contains("latent"))

  # extract all draws
  draws <- data.frame(t(fit.latent[["BUGSoutput"]][["sims.list"]][["latent"]]))
  item_id <- data[,1]
  alldraws <- cbind(item_id,draws)

  # create results object
  list(results = data,
       betas = betas,
       beta.plot = beta.plot,
       convergence = convergence,
       alldraws = alldraws,
       model.selection = model.choice)
}

latent.mod1 <- function(data.name){
  # load packages
  library(R2jags)
  library(tidyverse)
  data <- data.name
  # remove if all missing
  cols <- ncol(data)
  data$test <- apply(data[,2:cols],MARGIN=1,function(x) sum(is.na(x)))
  data <- data %>% filter(test!=cols-1) %>% dplyr::select(-test)
  items <- data %>% dplyr::select(-1)

  # set up data for latent model
  matrix <- items %>%
    mutate(across(everything(),function(x)c(scale(x)))) %>%
    as.matrix()
  list.data <- list(
    y = matrix,
    nrow = nrow(matrix),
    K = ncol(matrix)
  )
  K <- ncol(matrix)

  #model1
  mod1 <- function(){
    for (i in 1:nrow){
      latent[i] ~ dnorm(0,1)
      for(k in 1:K){
        y[i,k] ~ dnorm(mu[i,k], tau[k])
        mu[i,k] <- b[k]*latent[i]
      }
    }
    for(i in 1:K){
      b[i] ~ dnorm(0,.1);T(0, )
      tau[i]~dgamma(1, .1)
    }
  }

  #run jags
  latent.params <- c("latent", "b")
  fit.latent <- jags(data=list.data, inits=NULL, latent.params, model.file=mod1, n.chains=2, n.iter=4000, n.burnin=2000)
  x <- print(fit.latent)
  x.out<-data.frame(x$summary)

  # betas
  x.out$row <- rownames(x.out)
  varnames <- colnames(items)
  betas <- x.out %>% filter(grepl('b', row)) %>%
    rename(lwr = 3, med = 5, upr = 7) %>%
    mutate(var = varnames) %>%
    dplyr::select(var, lwr, med, upr)
  beta.plot <- ggplot(betas, aes(x=reorder(var, med), y=med, ymin=lwr, ymax=upr)) +
    geom_pointrange() + coord_flip() + ylab("Discrimination Parameter") + xlab("")

  # convergence
  convergence <- x.out %>% dplyr::select(row, Rhat, n.eff)

  # extract latent values
  beginning <- K+2
  end <- nrow(x.out)
  latent <- x.out[beginning:end,]
  data$latent <- latent$X50.
  data$latent_se <- latent$sd
  data$latent_lwr <- latent$X2.5.
  data$latent_upr <- latent$X97.5.
  data <- data %>% dplyr::select(1, contains("latent"))

  # extract all draws
  draws <- data.frame(t(fit.latent[["BUGSoutput"]][["sims.list"]][["latent"]]))
  item_id <- data[,1]
  alldraws <- cbind(item_id,draws)

  # create results object
  list(results = data,
       betas = betas,
       beta.plot = beta.plot,
       convergence = convergence,
       alldraws = alldraws)
}

latent.mod2 <- function(data.name){
  # load packages
  library(R2jags)
  library(tidyverse)
  data <- data.name

  # remove if all missing
  cols <- ncol(data)
  data$test <- apply(data[,2:cols],MARGIN=1,function(x) sum(is.na(x)))
  data <- data %>% filter(test!=cols-1) %>% dplyr::select(-test)
  items <- data %>% dplyr::select(-1)

  # set up data for latent model
  matrix <- items %>%
    mutate(across(everything(),function(x)c(scale(x)))) %>%
    as.matrix()
  list.data <- list(
    y = matrix,
    nrow = nrow(matrix),
    K = ncol(matrix)
  )
  K <- ncol(matrix)

  #model2
  mod2 <- function(){
    for (i in 1:nrow){
      latent[i] ~ dnorm(0,1)
      for(k in 1:K){
        y[i,k] ~ dnorm(mu[i,k], tau[k])
        mu[i,k] <- b[k]*latent[i]
      }
    }
    b[1] <- 1
    tau[1] ~ dgamma(1,.1)
    for(i in 2:K){
      b[i] ~ dnorm(0,.1)
      tau[i]~dgamma(1, .1)
    }
  }
  model.choice <- ifelse(length(table(sign(cor(items, use="pair"))>0))>1,mod2,mod1)

  #run jags
  latent.params <- c("latent", "b")
  fit.latent <- jags(data=list.data, inits=NULL, latent.params, model.file=mod2, n.chains=2, n.iter=4000, n.burnin=2000)
  x <- print(fit.latent)
  x.out<-data.frame(x$summary)

  # betas
  x.out$row <- rownames(x.out)
  varnames <- colnames(items)
  betas <- x.out %>% filter(grepl('b', row)) %>%
    rename(lwr = 3, med = 5, upr = 7) %>%
    mutate(var = varnames) %>%
    dplyr::select(var, lwr, med, upr)
  beta.plot <- ggplot(betas, aes(x=reorder(var, med), y=med, ymin=lwr, ymax=upr)) +
    geom_pointrange() + coord_flip() + ylab("Discrimination Parameter") + xlab("")

  # convergence
  convergence <- x.out %>% dplyr::select(row, Rhat, n.eff)

  # extract latent values
  beginning <- K+2
  end <- nrow(x.out)
  latent <- x.out[beginning:end,]
  data$latent <- latent$X50.
  data$latent_se <- latent$sd
  data$latent_lwr <- latent$X2.5.
  data$latent_upr <- latent$X97.5.
  data <- data %>% dplyr::select(1, contains("latent"))

  # extract all draws
  draws <- data.frame(t(fit.latent[["BUGSoutput"]][["sims.list"]][["latent"]]))
  item_id <- data[,1]
  alldraws <- cbind(item_id,draws)

  # create results object
  list(results = data,
       betas = betas,
       beta.plot = beta.plot,
       convergence = convergence,
       alldraws = alldraws)
}

linreg <- function(regdat){
  library(broom)
  library(tidyverse)
  library(stargazer)
  regdat <- regdat
  my_formula <- function(colPosition, trainSet){
    dep_part<- paste(colnames(trainSet)[colPosition],"~",sep=" ")
    ind_part<- paste(colnames(trainSet)[-colPosition],collapse=" + ")
    dt_formula<- as.formula(paste(dep_part,ind_part,sep=" "))
    return(dt_formula)
  }
  model.info <- my_formula( 1, regdat)
  mod <- lm(model.info, data=regdat)
  tidymod <- tidy(mod, conf.int=T) %>% mutate(sig = ifelse(p.value<0.05,1,0)) %>%
    filter(term!="(Intercept)")
  dvname <- colnames(regdat[1])
  modplot <- ggplot(tidymod, aes(x=term,y=estimate,ymin=conf.low, ymax=conf.high)) +
    geom_hline(yintercept=0, linetype="dotted") +
    xlab("") + ylab("") +
    geom_pointrange() + coord_flip() +
    theme_minimal() + theme(panel.border = element_rect(fill=NA, linewidth=0.5),
                            legend.position = "none") +
    ggtitle(paste("DV =", dvname, sep=""))
  print(summary(regdat))
  print(model.info)
  stargazer(mod, type="text")
  list(model = mod, tidymodel = tidymod, modelplot = modplot)
  plot(modelplot)
}
