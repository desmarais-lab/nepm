

# function to add the network effect variables
make_network_effects <- function(pdat,y,id,time,max_time_out){
  # pdat is a panel dataset
  # y is the character name of the numeric dependent variable
  # id is the cahracter name of the character unit/node id
  # time is the character name of the time period, as sequential integers

  # list unique nodes and possible ties among them
  node_id <- sort(unique(pdat[,id]))
  utimes <- sort(unique(pdat[,time]))

  possible_ties <- rbind(t(combn(node_id,2)), t(combn(node_id,2))[,c(2,1)])

  max_time <- max(table(pdat[,id]))
  possible_senders <- names(table(pdat[,id]))[which(table(pdat[,id]) >=
                                                      (max_time-max_time_out))]
  possible_ties <- possible_ties[is.element(possible_ties[,1],possible_senders),]

  # make a matrix to store lagged tie sender y values
  tie_x <- matrix(0,nrow(pdat),nrow(possible_ties))
  utimes <- sort(unique(pdat[,time]))

  # replace time with time sequence
  tseq <- (1:length(utimes))[match(pdat[,time],utimes)]

  # loop over times and sender nodes to set tie effect variables
  for(t in 1:length(utimes)){
    for(v in node_id){
      # find the value to be stored
      value <- pdat[which(pdat[,id]==v & tseq == (t-1)),y]
      # if v was not in data at t-1, set as NA
      if(length(value) < 1) value <- NA
      tie_x[which(tseq==t),which(possible_ties[,1] == v)] <- value
    }
  }

  colnames(tie_x) <- paste(possible_ties[,1],"_",possible_ties[,2],sep="")

  for(i in 1:ncol(tie_x)){
    tie_x[which(pdat[,id] != possible_ties[i,2]), i] <- 0
  }

  data.frame(pdat,tie_x)

}


#' A function to simulate a single time period from nepr
#'
#' @param ytm1 is a nodes x 1 vector of values of y in the previous time period
#' @param x is a nodes x k matrix where k is the number of covariates
#' @param beta is a (k + 1) x 1 vector of regression coefficients
#' @param gamma is a nties x 1 vector of tie effects
#' @param ties is a nties x 2 matrix where the columns are in order of sender/receiver
#' @param nodes, character vector of all of the vertex names. ids in ties should match
#' @param sig is the standard deviation of the error term
#' @return numeric vector the same length as ytm1 corresponding to the simulated outcomes.
#' @export
#' @examples
#' library(nepm)
#'
#' ## generative process for NEPM
#' # set some parameters
#' times <- 50 # time periods in data
#' nodes <- 50 # number of nodes
#' nties <- 20 # number of ties
#'
#' # simulate some data
#' set.seed(1001416)
#' node_id <- paste("n",1:nodes,sep="")
#' possible_ties <- rbind(t(combn(node_id,2)), t(combn(node_id,2))[,c(2,1)])
#' ties <- possible_ties[sample(1:nrow(possible_ties),nties),] # select the ties
#'
#' # make model parameter values
#' gamma <- .45 # set parameter to force stationarity
#' beta <- c(.5,-.5)
#' sig <- 1
#'
#' # generate a covariate that is not time-varying (easy to make time varying)
#' x <- rnorm(nodes)
#'
#'
#'
#' # matrix of simulated time periods
#' ytm1 <- rnorm(nodes,sd=sig)
#' panel_data <- NULL
#' for(t in 1:(1000+times)){
#'   yt <- simulate_time_period(ytm1,x,beta,gamma,ties,1)
#'   panel_data <- rbind(panel_data,t(yt))
#'   ytm1 <- yt
#' }
#'
#' # take last time steps as final data
#' panel_data <- panel_data[1001:(1000+times),]
#'
#' # long data
#' long_data <- NULL
#' for(i in 1:ncol(panel_data)){
#'   dati <- data.frame(id = rep(paste("n",i,sep=""),nrow(panel_data)),y=panel_data[,i],
#'                      time=1:nrow(panel_data),x=x[i])
#'   long_data <- rbind(long_data,dati)
#'
#' }
#'
#' system.time(nepm_test <- nepm(long_data,x_names="x",
#'                               y= "y",
#'                               id ="id",
#'                               time = "time"))
#'
#' nepm_estimate <- lm(nepm_test$nepm_formula,data=nepm_test$new_pdat)
simulate_time_period <- function(ytm1,x,beta,gamma,nodes,ties,sig){

  # make adjacency matrix of network effects
  amat <- matrix(0,length(nodes),length(nodes))
  rownames(amat) <- nodes
  colnames(amat) <- nodes
  amat[ties] <- gamma

  amat <- t(amat)

  # network effect
  net_effs <- amat%*%cbind(ytm1)

  # new y
  cbind(1,x)%*%cbind(beta) + net_effs + rnorm(length(nodes),sd=sig)

}





#' A function to run NEPM
#' @import abess foreach doParallel doRNG DoubleML glmnet
#' @param pdat The panel dataset as a dataframe.
#' @param x_names Character vector giving the names of the covariates. Should be column names in pdat.
#' @param y Character name of the dependent variable in pdat.
#' @param id Character name of the unit/node id in pdat. Should be a character variable.
#' @param time Character name of the numeric time variable in pdat.
#' @param max_time_out Integer, number of time periods a node can be out of the data in the beginning.
#' @param test_edges Logical, indicator of whether to test the significance of the edges (as a block)
#' @param nperm integer, number of permutations to use in the edge permutation test
#' @param ncore integer, number of cores to use in bootstrap
#' @return list with a character vector of edges inferred, a dataframe that can be
#' used to run nepm, a formula that combines the edges and the covariates in the model,
#' and other arguments used in the original call to nepm.
#' @export
#' @examples
#' library(nepm)
#'
#' ## generative process for NEPM
#' # set some parameters
#' times <- 50 # time periods in data
#' nodes <- 50 # number of nodes
#' nties <- 20 # number of ties
#'
#' # simulate some data
#' set.seed(1001416)
#' node_id <- paste("n",1:nodes,sep="")
#' possible_ties <- rbind(t(combn(node_id,2)), t(combn(node_id,2))[,c(2,1)])
#' ties <- possible_ties[sample(1:nrow(possible_ties),nties),] # select the ties
#'
#' # make model parameter values
#' gamma <- .45 # set parameter to force stationarity
#' beta <- c(.5,-.5)
#' sig <- 1
#'
#' # generate a covariate that is not time-varying (easy to make time varying)
#' x <- rnorm(nodes)
#'
#' # generate new time period given previous time point
#' sim_yt <- function(ytm1,x,beta,gamma,ties,sig){
#'   # ytm1 is a nodes x 1 vector
#'   # x is a nodes x k matrix where k is the number of covariates
#'   # beta is a (k + 1) x 1 vector of regression coefficients
#'   # gamma is a nties x 1 vector of tie effects
#'   # ties is a nties x 2 matrix where the columns are in order of sender/receiver
#'   # sigsq is the standard deviation of the error term
#'
#'   nodes <- length(ytm1)
#'
#'   # make adjacency matrix of network effects
#'   amat <- matrix(0,nodes,nodes)
#'   rownames(amat) <- paste("n",1:nodes,sep="")
#'   colnames(amat) <- paste("n",1:nodes,sep="")
#'   amat[ties] <- gamma
#'
#'   amat <- t(amat)
#'
#'   # network effect
#'   net_effs <- amat%*%cbind(ytm1)
#'
#'   # new y
#'   cbind(1,x)%*%cbind(beta) + net_effs + rnorm(nodes,sd=sig)
#'
#' }
#'
#'
#' # matrix of simulated time periods
#' ytm1 <- rnorm(nodes,sd=sig)
#' panel_data <- NULL
#' for(t in 1:(1000+times)){
#'   yt <- sim_yt(ytm1,x,beta,gamma,ties,1)
#'   panel_data <- rbind(panel_data,t(yt))
#'   ytm1 <- yt
#' }
#'
#' # take last time steps as final data
#' panel_data <- panel_data[1001:(1000+times),]
#'
#' # long data
#' long_data <- NULL
#' for(i in 1:ncol(panel_data)){
#'   dati <- data.frame(id = rep(paste("n",i,sep=""),nrow(panel_data)),y=panel_data[,i],
#'                      time=1:nrow(panel_data),x=x[i])
#'   long_data <- rbind(long_data,dati)
#'
#' }
#'
#' system.time(nepm_sim <- nepm(long_data,x_names="x",
#'                               y= "y",
#'                               id ="id",
#'                               time = "time"))
#'
#' nepm_estimate <- lm(nepm_sim$nepm_formula,data=nepm_sim$new_pdat)
nepm <- function(pdat,x_names,main_covs=x_names,y,id,time,max_time_out = 0,ncore=2){

  y_name <- y

  net_eff_data <- make_network_effects(pdat,y,id,time,max_time_out)

  x <- cbind(net_eff_data[,x_names],as.matrix(net_eff_data[,-(1:ncol(pdat))]))

  colnames(x)[1:length(x_names)] <- x_names

  y <- net_eff_data[,y_name]
  yx_full <- data.frame(y,x)

  x_sd <- apply(x,2,sd,na.rm=T)
  x <- x[,which(x_sd > 0)]
  yx <- data.frame(y,x)

  yx <- na.omit(yx)

  yx_sd <- apply(yx,2,sd)

  yx <- yx[,yx_sd > 0]

  #ridge_cv <- glmnet::cv.glmnet(as.matrix(yx[,-1]),yx[,1], alpha = 0)

  # Best lambda value
  #best_lambda <- ridge_cv$lambda.min

  abess_res <- abess::abess(yx[,-1],yx[,1])
                            #support.size = 0:(length(unique(pdat[,id])) + length(x_names)))
                            #lambda)

  var_names <- abess::extract(abess_res)$support.vars

  edges <- var_names[!is.element(var_names,names(pdat))]

  yx_id <- na.omit(net_eff_data[,c(id,y_name,x_names,edges)])

  obj_dml_data = DoubleML::DoubleMLData$new(na.omit(yx_full[,c("y",x_names,edges)]),
                                            y_col = "y", d_cols = main_covs)

  learner = mlr3::lrn("regr.cv_glmnet")
  ml_l_sim = learner$clone()
  ml_m_sim = learner$clone()


  obj_dml_plr_sim = DoubleML::DoubleMLPLR$new(obj_dml_data, ml_l = ml_l_sim, ml_m = ml_m_sim)
  obj_dml_plr_sim$fit(store_models=T)

  if(length(edges) > 0){

    n_est <- length(obj_dml_plr_sim$models$ml_l[[1]][[1]])
    n_vars <- length(obj_dml_plr_sim$models$ml_l)
    all_nz <- NULL
    for(i in 1:n_est){
      for(j in 1:n_vars){
        beta_param <- coef(obj_dml_plr_sim$models$ml_l[[j]][[1]][[i]]$model)
        beta_names <- row.names(beta_param)
        beta_param <- as.numeric(beta_param)
        beta_nm <- (beta_names[beta_param != 0])[-1]
        beta_nz <- (beta_param[beta_param != 0])[-1]
        names(beta_nz) <- beta_nm
        all_nz <- c(all_nz,beta_nz)
      }
    }

    #nz_ests <- table(all_nz)
    edges <-all_nz[is.element(names(all_nz),edges)]

  }

  return(list(edges=edges,new_pdat = net_eff_data[,union(names(pdat),names(edges))],
              nepm_formula = as.formula(paste(y_name,"~",
                                              paste(c(x_names,names(edges)),collapse="+")
                                              ,sep="")),
              x_names=x_names,
              y=y_name,
              id=id,
              time=time,
              max_time_out=max_time_out,doubleML_est = obj_dml_plr_sim))

}




#' A function to run a permutation test for the null of no edges in NEPM.
#' @import abess foreach doParallel doRNG
#' @param nepm_result an object resulting from nepm()
#' @param nperm integer, number of permutations to use in the edge permutation test
#' @param ncore integer, number of cores to use in bootstrap
#' @return a p-value for the test of the null hypothesis of no edges
#' as well as the null distribution of the number of edges inferred.
#' @export
#' @examples
#' library(nepm)
#'
#' ## generative process for NEPM
#' # set some parameters
#' times <- 50 # time periods in data
#' nodes <- 50 # number of nodes
#' nties <- 20 # number of ties
#'
#' # simulate some data
#' set.seed(1001416)
#' node_id <- paste("n",1:nodes,sep="")
#' possible_ties <- rbind(t(combn(node_id,2)), t(combn(node_id,2))[,c(2,1)])
#' ties <- possible_ties[sample(1:nrow(possible_ties),nties),] # select the ties
#'
#' # make model parameter values
#' gamma <- .45 # set parameter to force stationarity
#' beta <- c(.5,-.5)
#' sig <- 1
#'
#' # generate a covariate that is not time-varying (easy to make time varying)
#' x <- rnorm(nodes)
#'
#' # generate new time period given previous time point
#' sim_yt <- function(ytm1,x,beta,gamma,ties,sig){
#'   # ytm1 is a nodes x 1 vector
#'   # x is a nodes x k matrix where k is the number of covariates
#'   # beta is a (k + 1) x 1 vector of regression coefficients
#'   # gamma is a nties x 1 vector of tie effects
#'   # ties is a nties x 2 matrix where the columns are in order of sender/receiver
#'   # sigsq is the standard deviation of the error term
#'
#'   nodes <- length(ytm1)
#'
#'   # make adjacency matrix of network effects
#'   amat <- matrix(0,nodes,nodes)
#'   rownames(amat) <- paste("n",1:nodes,sep="")
#'   colnames(amat) <- paste("n",1:nodes,sep="")
#'   amat[ties] <- gamma
#'
#'   amat <- t(amat)
#'
#'   # network effect
#'   net_effs <- amat%*%cbind(ytm1)
#'
#'   # new y
#'   cbind(1,x)%*%cbind(beta) + net_effs + rnorm(nodes,sd=sig)
#'
#' }
#'
#'
#' # matrix of simulated time periods
#' ytm1 <- rnorm(nodes,sd=sig)
#' panel_data <- NULL
#' for(t in 1:(1000+times)){
#'   yt <- sim_yt(ytm1,x,beta,gamma,ties,1)
#'   panel_data <- rbind(panel_data,t(yt))
#'   ytm1 <- yt
#' }
#'
#' # take last time steps as final data
#' panel_data <- panel_data[1001:(1000+times),]
#'
#' # long data
#' long_data <- NULL
#' for(i in 1:ncol(panel_data)){
#'   dati <- data.frame(id = rep(paste("n",i,sep=""),nrow(panel_data)),y=panel_data[,i],
#'                      time=1:nrow(panel_data),x=x[i])
#'   long_data <- rbind(long_data,dati)
#'
#' }
#'
#' system.time(nepm_sim <- nepm(long_data,x_names="x",
#'                               y= "y",
#'                               id ="id",
#'                               time = "time"))
#'
#' nepm_test <- test_nepm(nepm_sim)
test_nepm <- function(nepm_result,nperm=50,ncore=2){

  pdat <- nepm_result$new_pdat[,c(nepm_result$y,
                                  nepm_result$id,
                                  nepm_result$time,
                                  nepm_result$x_names)]

  y_name <- nepm_result$y

  net_eff_data <- make_network_effects(pdat,
                                       nepm_result$y,
                                       nepm_result$id,
                                       nepm_result$time,
                                       nepm_result$max_time_out)

  null_edges <- NULL

  y_neteff <- net_eff_data[,y_name]

  for(p in 1:nperm){

    net_eff_data[,y_name] <- sample(y_neteff,length(y_neteff))

    x <- cbind(net_eff_data[,nepm_result$x_names],as.matrix(net_eff_data[,-(1:ncol(pdat))]))

    colnames(x)[1:length(nepm_result$x_names)] <- nepm_result$x_names

    x_sd <- apply(x,2,sd,na.rm=T)
    x <- x[,which(x_sd > 0)]

    y <- net_eff_data[,y_name]

    yx <- data.frame(y,x)

    yx <- na.omit(yx)

    yx_sd <- apply(yx,2,sd)

    yx <- yx[,yx_sd > 0]

    abess_res <- abess::abess(yx[,-1],yx[,1],
                              support.size = 0:(length(unique(pdat[,nepm_result$id]))
                                                + length(nepm_result$x_names)))

    var_names <- abess::extract(abess_res)$support.vars

    edges <- var_names[!is.element(var_names,names(pdat))]

    yx_id <- na.omit(net_eff_data[,c(nepm_result$id,
                                     nepm_result$y,
                                     nepm_result$x_names,edges)])

    # screen_est <- SIS::SIS(as.matrix(yx[,-1]),yx[,1])

    # edges <- intersect(c(x_names,edges)[screen_est$ix],edges)

    if(length(edges) > 0){

      boot_iter <- 80
      boot_res <- NULL
      uid <- unique(yx_id[,nepm_result$id])

      cl <- parallel::makeCluster(ncore)
      doParallel::registerDoParallel(cl)

      `%dorng%` <- doRNG::`%dorng%`

      boot_res <- foreach::foreach(i = 1:boot_iter,.options.RNG=9202011) %dorng% {
        boot_id <- NULL
        for(u in uid){
          row_u <- which(yx_id[,nepm_result$id]==u)
          boot_id <- c(boot_id,sample(row_u,length(row_u),rep=T))
        }
        yx_id_b <- yx_id[boot_id,]
        x_b <- as.matrix(yx_id_b[,c(nepm_result$x_names,edges)])
        y_b <- yx_id_b[,nepm_result$y]
        cv.fit <- glmnet::cv.glmnet(x_b,y_b)
        eff_b <- rownames(coef(cv.fit, s = "lambda.min"))[as.numeric(coef(cv.fit, s = "lambda.min")) !=0]
        eff_b
      }

      parallel::stopCluster(cl)

      boot_res <- c(unlist(boot_res))

      eff_freq <- table(boot_res)

      est_edges <- names(eff_freq[which(eff_freq==boot_iter)])

      edges <- intersect(est_edges,edges)

    }

    null_edges <- c(null_edges,length(edges))

  }

  return(list(p_value = mean(null_edges >= length(nepm_result$edges)),
              null_dist = null_edges))

}



#' A function to compare forecast performance of nepm to lm
#'
#' @param pdat The panel dataset as a dataframe.
#' @param x_names Character vector giving the names of the covariates. Should be column names in pdat.
#' @param y Character name of the dependent variable in pdat.
#' @param id Character name of the unit/node id in pdat. Should be a character variable.
#' @param time Character name of the numeric time variable in pdat.
#' @param boot Indicator of whether to use bootstrapping to calculate uncertainty measures.
#' @param nboot Integer, number of bootstrap iterations.
#' @param max_time_out Integer, number of time periods a node can be out of the data in the beginning.
#' @param start_tim Integer, the time period at which to start the forecasting experiment. Defaults to half the time points.
#' @return list with a character vector of edges inferred, a dataframe that can be used to run nepm, and a formula that combines the edges and the covariates in the model.
#' @export
#' @examples
#' library(nepm)
#'
#' ## generative process for NEPM
#' # set some parameters
#' times <- 50 # time periods in data
#' nodes <- 50 # number of nodes
#' nties <- 20 # number of ties
#'
#' # simulate some data
#' set.seed(1001416)
#' node_id <- paste("n",1:nodes,sep="")
#' possible_ties <- rbind(t(combn(node_id,2)), t(combn(node_id,2))[,c(2,1)])
#' ties <- possible_ties[sample(1:nrow(possible_ties),nties),] # select the ties
#'
#' # make model parameter values
#' gamma <- .45 # set parameter to force stationarity
#' beta <- c(.5,-.5)
#' sig <- 1
#'
#' # generate a covariate that is not time-varying (easy to make time varying)
#' x <- rnorm(nodes)
#'
#' # generate new time period given previous time point
#' sim_yt <- function(ytm1,x,beta,gamma,ties,sig){
#'   # ytm1 is a nodes x 1 vector
#'   # x is a nodes x k matrix where k is the number of covariates
#'   # beta is a (k + 1) x 1 vector of regression coefficients
#'   # gamma is a nties x 1 vector of tie effects
#'   # ties is a nties x 2 matrix where the columns are in order of sender/receiver
#'   # sigsq is the standard deviation of the error term
#'
#'   nodes <- length(ytm1)
#'
#'   # make adjacency matrix of network effects
#'   amat <- matrix(0,nodes,nodes)
#'   rownames(amat) <- paste("n",1:nodes,sep="")
#'   colnames(amat) <- paste("n",1:nodes,sep="")
#'   amat[ties] <- gamma
#'
#'   amat <- t(amat)
#'
#'   # network effect
#'   net_effs <- amat%*%cbind(ytm1)
#'
#'   # new y
#'   cbind(1,x)%*%cbind(beta) + net_effs + rnorm(nodes,sd=sig)
#'
#' }
#'
#'
#' # matrix of simulated time periods
#' ytm1 <- rnorm(nodes,sd=sig)
#' panel_data <- NULL
#' for(t in 1:(1000+times)){
#'   yt <- sim_yt(ytm1,x,beta,gamma,ties,1)
#'   panel_data <- rbind(panel_data,t(yt))
#'   ytm1 <- yt
#' }
#'
#' # take last time steps as final data
#' panel_data <- panel_data[1001:(1000+times),]
#'
#' # long data
#' long_data <- NULL
#' for(i in 1:ncol(panel_data)){
#'   dati <- data.frame(id = rep(paste("n",i,sep=""),nrow(panel_data)),y=panel_data[,i],
#'                      time=1:nrow(panel_data),x=x[i])
#'   long_data <- rbind(long_data,dati)
#'
#' }
#'
#' fit_result <- forecast_comparison(long_data,x_names="x",
#'                               y= "y",
#'                               id ="id",
#'                               time = "time")
forecast_comparison <- function(pdat,x_names,y,id,time,
                                max_time_out = 0,start_time=NULL,ncore=2){

  utimes <- sort(unique(pdat[,time]))

  if(is.null(start_time)) start_time <- utimes[ceiling(length(utimes)*0.7)]

  y_name <- y

  net_eff_data <- make_network_effects(pdat,y,id,time,max_time_out)

  x <- cbind(net_eff_data[,x_names],as.matrix(net_eff_data[,-(1:ncol(pdat))]))

  colnames(x)[1:length(x_names)] <- x_names

  x_sd <- apply(x,2,sd,na.rm=T)
  x <- x[,which(x_sd > 0)]

  y <- net_eff_data[,y_name]

  tyx <- data.frame(t=net_eff_data[,time],id = net_eff_data[,id],y,x)

  names(tyx)[2] <- id
  names(tyx)[3] <- y_name

  tyx <- na.omit(tyx)

  test_times <- utimes[which(utimes >= start_time)]

  mae_nepm <- NULL
  mae_lm <- NULL

  for(t in test_times){

    set.seed(9202011)

    tyx_train <- tyx[which(tyx$t < t),]
    tyx_test <- tyx[which(tyx$t==t),]
    zsd <- which(apply(tyx_train,2,sd)==0)

    abess_res <- abess::abess(tyx_train[,-c(1,2,3,zsd)],tyx_train[,3])

    var_names <- abess::extract(abess_res)$support.vars

    edges <- var_names[!is.element(var_names,names(pdat))]

    yx <- cbind(tyx_train[,3],tyx_train[,var_names])

    if(length(edges) > 0){
      boot_iter <- 80
      boot_res <- NULL
      uid <- unique(tyx_train[,id])

      cl <- parallel::makeCluster(ncore)
      doParallel::registerDoParallel(cl)

      `%dorng%` <- doRNG::`%dorng%`

      boot_res <- foreach::foreach(i = 1:boot_iter,.options.RNG=9202011) %dorng% {
        boot_id <- NULL
        for(u in uid){
          row_u <- which(tyx_train[,id]==u)
          boot_id <- c(boot_id,sample(row_u,length(row_u),rep=T))
        }
        yx_id_b <- tyx_train[boot_id,]
        x_b <- as.matrix(yx_id_b[,c(x_names,edges)])
        y_b <- yx_id_b[,y_name]
        cv.fit <- glmnet::cv.glmnet(x_b,y_b)
        eff_b <- rownames(coef(cv.fit, s = "lambda.min"))[as.numeric(coef(cv.fit, s = "lambda.min")) !=0]
        eff_b
      }

      parallel::stopCluster(cl)

      boot_res <- c(unlist(boot_res))

      eff_freq <- table(boot_res)

      est_edges <- names(eff_freq[which(eff_freq==boot_iter)])

      edges <- intersect(est_edges,edges)

    }

    nepm_formula <- as.formula(paste(y_name,"~",paste(c(x_names,edges),collapse="+"),sep=""))

    lm_formula <- as.formula(paste(y_name,"~",paste(c(x_names),collapse="+"),sep=""))

    fit_nepm <- lm(nepm_formula,data=tyx_train)
    fit_lm <- lm(lm_formula,data=tyx_train)

    pred_nepm <- predict(fit_nepm,newdata=tyx_test)
    pred_lm <- predict(fit_lm,newdata=tyx_test)

    y_test <- tyx_test[,3]

    mae_nepm <- c(mae_nepm,mean(abs(y_test-pred_nepm)))
    mae_lm <- c(mae_lm,mean(abs(y_test-pred_lm)))

  }

  list(fit_nepm = mae_nepm,fit_lm = mae_lm)

}




