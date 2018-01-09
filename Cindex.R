Cindex <- function (sigma = 0.1) {
  approxGrad <- function(x) { ## sigmoid function for gradient
    exp(x/sigma) / (sigma * (1 + exp(x/sigma))^2)
  }
  approxLoss <- function(x) { ## sigmoid function for loss
    1 / (1 + exp(x / sigma))
  }
  compute_weights <- function(y, w = 1){ ## compute weights
    ipcw_wow <- IPCweights(y[w != 0,])
    ipcw <- numeric(nrow(y))
    ipcw[w!=0] <- ipcw_wow
    survtime <- y[,1]
    n <- nrow(y)
    wweights <- matrix( (ipcw)^2, nrow = n, ncol = n)
    weightsj <- matrix(survtime, nrow = n, ncol = n)
    weightsk <- matrix(survtime, nrow = n, ncol = n, byrow = TRUE)
    weightsI <- (weightsj < weightsk) + 0
    wweights <- wweights * weightsI
    Wmat <- w %o% w
    wweights <- wweights * Wmat
    wweights <- wweights / sum(wweights)
    rm(weightsI); rm(weightsk); rm(weightsj)
    return(wweights)
  }
  ngradient = function(y, f, w = 1) { ## negative gradient
    if (!all(w %in% c(0,1)))
      stop(sQuote("weights"), " must be either 0 or 1 for family ",
           sQuote("UnoC"))
    survtime <- y[,1]
    event <- y[,2]
    if (length(w) == 1) w <- rep(1, length(event))
    if (length(f) == 1) {
      f <- rep(f, length(survtime))
    }
    n <- length(survtime)
    etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
    etak <- matrix(f, nrow = n, ncol = n)
    etaMat <- etak - etaj
    rm(etaj); rm(etak);
    weights_out <- compute_weights(y, w)
    M1 <- approxGrad(etaMat) * weights_out
    ng <- colSums(M1) - rowSums(M1)
    return(ng)
  }
  risk = function(y, f, w = 1) { ## empirical risk
    survtime <- y[,1]
    event <- y[,2]
    if (length(f) == 1) {
      f <- rep(f, length(y))
    }
    n <- length(survtime)
    etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
    etak <- matrix(f, nrow = n, ncol = n)
    etaMat <- (etak - etaj)
    rm(etaj); rm(etak);
    weights_out <- compute_weights(y, w)
    M1 <- approxLoss(etaMat) * weights_out
    return(- sum(M1))
  }
  Family(
    ngradient = ngradient,
    risk = risk,
    weights = "zeroone",
    offset = function(y, w = 1) {0},
    check_y = function(y) {
      if (!inherits(y, "Surv"))
        stop("response is not an object of class ", sQuote("Surv"),
             " but ", sQuote("family = UnoC()"))
      y},
    rclass = function(f){},
    name = paste("Concordance Probability by Uno")
  )
}
