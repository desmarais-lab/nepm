

# function to add the network effect variables
make_network_effects <- function(pdat,y,id,time){
  # pdat is a panel dataset
  # y is the character name of the numeric dependent variable
  # id is the cahracter name of the character unit/node id
  # time is the character name of the time period, as sequential integers

  # list unique nodes and possible ties among them
  node_id <- sort(unique(pdat[,id]))
  possible_ties <- rbind(t(combn(node_id,2)), t(combn(node_id,2))[,c(2,1)])

  # make a matrix to store lagged tie sender y values
  tie_x <- matrix(0,nrow(pdat),nrow(possible_ties))
  utimes <- sort(unique(pdat[,time]))

  # loop over times and sender nodes to set tie effect variables
  for(t in utimes){
    for(v in node_id){
      # find the value to be stored
      value <- pdat[which(pdat[,id]==v & pdat[,time] == (t-1)),y]
      # if v was not in data at t-1, set as NA
      if(length(value) < 1) value <- NA
      tie_x[which(pdat[,time]==t),which(possible_ties[,1] == v)] <- value
    }
  }

  colnames(tie_x) <- paste(possible_ties[,1],"_",possible_ties[,2],sep="")

  for(i in 1:ncol(tie_x)){
    tie_x[which(pdat[,id] != possible_ties[i,2]), i] <- 0
  }

  data.frame(pdat,tie_x)

}



#' A function to run NEPM
#'
#' @param pdat The panel dataset as a dataframe.
#' @param x_names Character vector giving the names of the covariates. Should be column names in pdat.
#' @param y Character name of the dependent variable in pdat.
#' @param id Character name of the unit/node id in pdat. Should be a character variable.
#' @param time Character name of the numeric time variable in pdat.
#' @param boot Indicator of whether to use bootstrapping to calculate uncertainty measures.
#' @param nboot Integer, number of bootstrap iterations.
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
#' system.time(nepm_test <- nepm(long_data,x_names="x",
#'                               y= "y",
#'                               id ="id",
#'                               time = "time"))
#'
#' nepm_estimate <- lm(nepm_test$nepm_formula,data=nepm_test$new_pdat)
nepm <- function(pdat,x_names,y,id,time,boot=F,nboot=500){
  require(abess)
  require(BMisc)

  y_name <- y

  net_eff_data <- make_network_effects(pdat,y,id,time)

  x <- cbind(net_eff_data[,x_names],as.matrix(net_eff_data[,-(1:ncol(pdat))]))

  colnames(x)[1:length(x_names)] <- x_names

  y <- net_eff_data$y

  yx <- na.omit(data.frame(y,x))

  abess_res <- abess::abess(yx[,-1],yx[,1],always.include = 1:length(x_names))

  var_names <- extract(abess_res)$support.vars

  edges <- var_names[!is.element(var_names,names(pdat))]

  list(edges=edges,new_pdat = net_eff_data[,union(names(pdat),edges)],
       nepm_formula = as.formula(paste(y_name,"~",paste(c(x_names,edges),collapse="+"),sep="")))

}





