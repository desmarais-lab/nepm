

# function to add the network effect variables
make_network_effects <- function(pdat,y,id,time,max_time_out){
  # pdat is a panel dataset
  # y is the character name of the numeric dependent variable
  # id is the cahracter name of the character unit/node id
  # time is the character name of the time period, as sequential integers

  # list unique nodes and possible ties among them
  node_id <- sort(unique(pdat[,id]))
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
  amat <- matrix(0,nodes,nodes)
  rownames(amat) <- nodes
  colnames(amat) <- nodes
  amat[ties] <- gamma

  amat <- t(amat)

  # network effect
  net_effs <- amat%*%cbind(ytm1)

  # new y
  cbind(1,x)%*%cbind(beta) + net_effs + rnorm(nodes,sd=sig)

}





#' A function to run NEPM
#'
#' @param pdat The panel dataset as a dataframe.
#' @param x_names Character vector giving the names of the covariates. Should be column names in pdat.
#' @param y Character name of the dependent variable in pdat.
#' @param id Character name of the unit/node id in pdat. Should be a character variable.
#' @param time Character name of the numeric time variable in pdat.
#' @param max_time_out Integer, number of time periods a node can be out of the data in the beginning.
#' @param test_edges Logical, indicator of whether to test the significance of the edges (as a block)
#' @param nperm integer, number of permutations to use in the edge permutation test
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
nepm <- function(pdat,x_names,y,id,time,max_time_out = 0,
                 test_edges=F,nperm=100){

  set.seed(9202011)

  y_name <- y

  net_eff_data <- make_network_effects(pdat,y,id,time,max_time_out)

  x <- cbind(net_eff_data[,x_names],as.matrix(net_eff_data[,-(1:ncol(pdat))]))

  colnames(x)[1:length(x_names)] <- x_names

  x_sd <- apply(x,2,sd,na.rm=T)
  x <- x[,which(x_sd > 0)]

  y <- net_eff_data[,y_name]

  yx <- data.frame(y,x)

  yx <- na.omit(yx)

  abess_res <- abess::abess(yx[,-1],yx[,1],support.size = 0:(length(unique(pdat[,id])) + length(x_names)))

  var_names <- abess::extract(abess_res)$support.vars

  edges <- var_names[!is.element(var_names,names(pdat))]

  if(test_edges){
    null_edges <- NULL
    for(i in 1:nperm){
      y_p <- yx[,1]
      all_x_p <- yx[sample(1:nrow(yx),nrow(yx)),-1]
      abess_res_p <- abess::abess(all_x_p,y_p,support.size = 0:(length(unique(pdat[,id])) + length(x_names)))
      var_names <- abess::extract(abess_res_p)$support.vars
      edges_p <- var_names[!is.element(var_names,names(pdat))]
      null_edges <- c(null_edges,length(edges_p))
    }

    return(list(edges=edges,new_pdat = net_eff_data[,union(names(pdat),edges)],
         nepm_formula = as.formula(paste(y_name,"~",paste(c(x_names,edges),collapse="+"),sep="")),test_p=mean(null_edges >= length(edges))))

  }

  if(!test_edges){
    return(list(edges=edges,new_pdat = net_eff_data[,union(names(pdat),edges)],
       nepm_formula = as.formula(paste(y_name,"~",paste(c(x_names,edges),collapse="+"),sep=""))))
  }

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
forecast_comparison <- function(pdat,x_names,y,id,time,
                                max_time_out = 0,start_time=NULL){

  utimes <- sort(unique(pdat[,time]))

  if(is.null(start_time)) start_time <- utimes[ceiling(length(utimes)*0.7)]

  y_name <- y

  net_eff_data <- make_network_effects(pdat,y,id,time,max_time_out)

  x <- cbind(net_eff_data[,x_names],as.matrix(net_eff_data[,-(1:ncol(pdat))]))

  colnames(x)[1:length(x_names)] <- x_names

  x_sd <- apply(x,2,sd,na.rm=T)
  x <- x[,which(x_sd > 0)]

  y <- net_eff_data[,y_name]

  tyx <- data.frame(t=net_eff_data[,time],y,x)

  tyx <- na.omit(tyx)

  test_times <- utimes[which(utimes >= start_time)]

  mae_nepm <- NULL
  mae_lm <- NULL

  for(t in test_times){

    set.seed(9202011)

    tyx_train <- tyx[which(tyx$t < t),]
    tyx_test <- tyx[which(tyx$t==t),]
    zsd <- which(apply(tyx_train,2,sd)==0)

    abess_res <- abess::abess(tyx_train[,-c(1,2,zsd)],tyx_train[,2],support.size = 0:(length(unique(pdat[,id])) + length(x_names)))

    var_names <- abess::extract(abess_res)$support.vars

    edges <- var_names[!is.element(var_names,names(pdat))]

    nepm_formula <- as.formula(paste("y","~",paste(c(x_names,edges),collapse="+"),sep=""))

    lm_formula <- as.formula(paste("y","~",paste(c(x_names),collapse="+"),sep=""))

    fit_nepm <- lm(nepm_formula,data=tyx_train)
    fit_lm <- lm(lm_formula,data=tyx_train)

    pred_nepm <- predict(fit_nepm,newdata=tyx_test)
    pred_lm <- predict(fit_lm,newdata=tyx_test)

    y_test <- tyx_test[,2]

    mae_nepm <- c(mae_nepm,mean(abs(y_test-pred_nepm)))
    mae_lm <- c(mae_lm,mean(abs(y_test-pred_lm)))

  }

  list(fit_nepm = mae_nepm,fit_lm = mae_lm)

}




