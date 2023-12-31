#Gaussian Kernel
KSE <- function( x1, x2, l ) {
  KSE <- exp( -0.5 * sum((x1 - x2)/l)^2 )
}

#Kernel lineare (no-kernel)
Klin <- function( x1, x2 ) {
  Klin <- sum(x1*x2)
}

# funzione obiettivo 
eval_f <- function(x, delta, delta_star, nu ,m, xx) {
  #definizione dei parametri da ottimizzare 
  alpha <- x[1:(length(x)/2)]
  alpha_star <- x[-c(1:(length(x)/2))]
  
  stopifnot(length(alpha)==length(alpha_star))
  
  K <- matrix(1,length(alpha),length(alpha_star))
  for( i in 1:(m-1) )
    for( j in (i+1):m ) {
      
      # uso il Gaussian Kernel
      K[i,j] <- KSE(xx[i,],xx[j,], l=0.01 ) 
      #K[i,j] <- KSE(xx[i],xx[j], l=0.01 ) in caso di dati unidimensionali
      
      # uso il Kernel lineare
      # K[i,j] = Klin(X[i,],X[j,])
      
      K[j,i] <- K[i,j]
    }
  
  tot <- 0
  for( i in 1:m )
    for( j in 1:m )
      tot <- tot + (alpha[i]-alpha_star[i]) *
    (alpha[j]-alpha_star[j]) * K[i,j]
  
  tot <- 0.5 * tot  - delta*sum(alpha) + delta_star*sum(alpha_star)
  
  obj <- round(tot,8)
}

# vincolo d'uguaglianza
eval_g_eq <- function(x, delta, delta_star, nu ,m, xx) {
  alpha <- x[1:(length(x)/2)]
  alpha_star <- x[-c(1:(length(x)/2))]
  eqCons <- round(sum( alpha - alpha_star)-1,8)
}

#parametri settati per un'esempio
m <- 10 
nu <- 0.1
alpha <- round(rep(1/(nu*m),m),8)
alpha_star <- round(rep((1-nu)/(nu*m),m),8)
x0 <- c(alpha,alpha_star)

# check per verificare la feasibility iniziale:
if( all(x0>=0 & x0<=round(1/(nu*m),8) )) {
  cat(">Vincoli di disuguaglianza: OK!\n")
} else {
  stop("ERROR: Violazione dei vincoli di disuguaglianza!!!!!")
}

if( eval_g_eq(x0)==0 ) {
  cat("> Vincolo di uguaglianza: OK!\n")
} else {
  stop("ERROR: Violazione del vincolo di uguaglianza!!!!!")
}



cat("> Inizio ottimizzazione (ISRES)...\n")
res1 <- isres( x0=x0, fn=eval_f, lower=rep(0,2*m), 
               upper=rep(1/(nu*m),2*m), heq=eval_g_eq )
print(res1)

cat("> Inizio ottimizzazione (SLSQP)...\n")
res2 <- slsqp( x0=x0, fn=eval_f, gr=NULL, lower=rep(0,2*m), 
               upper=rep(1/(nu*m),2*m), heq=eval_g_eq )
print(res2)

# vincolo d'uguaglianza come due vincoli di disuguaglianza
eqAsIneqCons <- function( x ) {
  alpha <- x[1:(length(x)/2)]
  alpha_ <- x[-c(1:(length(x)/2))]
  S <- round(sum(alpha-alpha_),8)
  eqAsIneqCons <- c( S-1, 1-S ) 
}

cat("> Inizio ottimizzazione (COBYLA)...\n")
res3 <- cobyla( x0=x0, fn=eval_f, lower=rep(0,2*m),
                upper=rep(1/(nu*m),2*m), hin=eqAsIneqCons )
print(res3)

cat("> Inizio ottimizzazione (AUGLAG)...\n")
res4 <- auglag(x0=x0, fn=eval_f, gr = NULL , lower=rep(0,2*m), 
               upper=rep(1/(nu*m),2*m), heq=eval_g_eq, localsolver = "SLSQP") 
print(res4)

#Ottimizzazione degli iperparametri

obj.fun <- makeSingleObjectiveFunction(
  name = "svm",
  fn = function(x) {
    dati <- dati 
    m <- length(dati)
    delta <- delta
    alpha <- alpha
    alpha_star <- alpha_star
    alphas <- c(alpha, alpha_star)
    
    res <- slsqp( x0=x0, fn=eval_f, gr=NULL, 
                  lower=rep(0,2*m), upper=round(rep(1/(nu*m),2*m),8),
                  heq=eval_g_eq, nl.info = T, 
                  control = list(maxeval = 10), 
                  delta = delta,
                  delta_star = x[2],
                  nu = x[1],
                  m = m,
                  xx = dati)
    
    obj.fun <- res$value
  },
  par.set = makeParamSet(
    makeNumericParam("nu", lower = 0, upper = 1),
    makeNumericParam("delta_star", lower = 0, upper = 10)
  ),
  minimize = TRUE
)

des = generateDesign(n = 5, par.set = getParamSet(obj.fun),
                     fun = lhs::randomLHS)

surr.km = makeLearner("regr.km", predict.type = "se",
                      covtype = "matern3_2",
                      control = list(trace = FALSE))

control = makeMBOControl()
control = setMBOControlTermination(control, iters = 10)
control = setMBOControlInfill(control, crit = makeMBOInfillCritEI())

run2 = mbo(obj.fun, design = des, learner = surr.km,
           control = control, show.info = TRUE)
print(run) 