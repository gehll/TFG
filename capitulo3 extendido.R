
run_us <- bvar(varmat_us, lags = 2, 
               n_draw = 50000, n_burn = 25000,
               priors = priors, mh = mh)


summary(run_us)



## Estimación de los modelos VECM y BVECM
"
En este apartado vamos a analizar las series simuladas que son cointegradas, y los datos que recogen el LIBOR a 12 meses y los swaps OVERNIGHT a 9 meses. Primero, vamos a estimar el modeo VECM frecuentista a los datos simulados y después el modelo bayesiano para comparar los parámetros estimados con los *reales*. Recordar que los parámetros reales con los que se han simulado los datos siguen la forma de \ref{eq:simvec}. Se va a usar la librería `tsDyn`, y, en específico la función `VECM` para estimar el modelo frecuentista. Como  método de estimación se usará la máxima verosimilitud de Johansen ya que produce mejores estimaciones que mínimos cuadrados ordinarios en dos fases. Como sabemos que los datos simulados se han creado con una estructura autorregresiva de 2 retardos, lo especificaremos en los parámetros del modelo.
"

library(tsDyn)
set.seed(100)
library(mnormt)

innov<-rmnorm(1000, varcov=diag(2))
Bvecm <- rbind(c(-0.3, 0.2,-0.3, 0.2, 0.1), c(0.2, 0.4, -0.2, 0.1, 0.25))
datos_VEC <- VECM.sim(B=Bvecm,  n=1000, beta=1.5, lag=2,include="none", innov=innov)

vecm_datossim = VECM(datos_VEC, lag = 2, include = "none", estim = "ML")
summary(vecm_datossim)

"
Las estimaciones son muy precisas, sobretodo para el parámetro $\beta$ de la relación de largo plazo. En general todos los coeficientes estimados son muy cercanos a los reales reflejados en \ref{eq:simvec} y son significativos al 5% como era de esperar. 

Ahora, vamos a estimar estos coeficientes mediante el método bayesiano. La librería `bvartools` permite estimar modelos BVECM. Con la función `gen_vec()` se crea el objeto de clase *BVEC* donde se introduce la estructura del modelo que se quiere estimar (número de retardos, parte determinista, iteraciones del muestreador y número de burnin). Después, ese objeto que se crea se para a la función `add_priors()` donde se añaden las a prioris al objeto de tipo *BVEC*. En este caso, se van a utilizar a prioris no informativas para ver qué tal resultados da la estimación cuando no hay un componente importante de información a priori a la hora de estimar los parámetros. Después de crear el objeto *BVEC* con las a prioris definidas, se puede utilizar la función `draw_posterior()` para obtener las estimaciones a posteriori. En este caso, se ha desarrollado un código específico para la estimación a partir un ejemplo que hay dentro de las viñetas que contiene la librería. Se ha hecho esto porque la estimación a posteriori con la función `draw_posterior()` tardaba demasiado, y, teniendo en cuenta que más adelante se va a realizar un walk forward validantion, esto hacía que fuese computacionalmente imposible de hacer con la función `draw-posterior()`.
"

library(bvartools)

datos_VEC <- ts(datos_VEC)
data <- gen_vec(datos_VEC, p = 3, r=1, iterations = 5000, burnin = 2500)
data <- add_priors(data,
                   coint = list(v_i = 0, p_tau_i = 1),
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(df = 0, scale = .0001))

# Reset random number generator for reproducibility
set.seed(100)

# Obtain data matrices
y <- t(data$data$Y)
w <- t(data$data$W)
x <- t(data$data$X)

r <- data$model$rank # Set rank

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
k_gamma <- k * k_x # Total number of non-cointegration coefficients

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta

# Priors
a_mu_prior <- data$priors$noncointegration$mu # Prior means
a_v_i_prior <- data$priors$noncointegration$v_i # Inverse of the prior covariance matrix

v_i <- data$priors$cointegration$v_i
p_tau_i <- data$priors$cointegration$p_tau_i

sigma_df_prior <- data$priors$sigma$df # Prior degrees of freedom
sigma_scale_prior <- data$priors$sigma$scale # Prior covariance matrix
sigma_df_post <- tt + sigma_df_prior # Posterior degrees of freedom

# Initial values
beta <- matrix(0, k_w, r)
beta[1:r, 1:r] <- diag(1, r)

sigma_i <- diag(1 / .0001, k)

g_i <- sigma_i

iterations <- data$model$iterations # Number of iterations of the Gibbs sampler
burnin <- data$model$burnin # Number of burn-in draws
draws <- iterations + burnin # Total number of draws

# Data containers
draws_alpha <- matrix(NA, k_alpha, iterations)
draws_beta <- matrix(NA, k_beta, iterations)
draws_pi <- matrix(NA, k * k_w, iterations)
draws_gamma <- matrix(NA, k_gamma, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = sigma_i,
                         v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                         gamma_mu_prior = a_mu_prior,
                         gamma_v_i_prior = a_v_i_prior)
  alpha <- temp$alpha
  beta <- temp$beta
  Pi <- temp$Pi
  gamma <- temp$Gamma
  
  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  sigma_scale_post <- solve(tcrossprod(u) + v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_post)[,, 1], k)
  sigma <- solve(sigma_i)
  
  # Update g_i
  g_i <- sigma_i
  
  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_pi[, draw - burnin] <- Pi
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- sigma
  }
}

#Sacamos la estimación de los coeficientes de cointegración
beta <- apply(t(draws_beta) / t(draws_beta)[, 1], 2, mean) # Obtain means for every row
beta <- matrix(beta, k_w) # Transform mean vector into a matrix
beta <- round(beta, 3) # Round values
dimnames(beta) <- list(dimnames(w)[[1]], NULL) # Rename matrix dimensions

beta

# Number of non-deterministic coefficients
k_nondet <- (k_x - 0) * k

# Generate bvec object
bvec_est <- bvec(y = data$data$Y,
                 w = data$data$W,
                 x = data$data$X[, 1:4],
                 x_d = data$data$X[, -(1:4)],
                 Pi = draws_pi,
                 r = 1,
                 Gamma = draws_gamma[1:k_nondet,],
                 Sigma = draws_sigma)

summary(bvec_est)

"
La estimación de la relación de largo plazo es muy buena (-1.497 frente a -1.5 real) igual que con el método frecuentista. En cuanto a los términos de corección de equilibrio, los IC al 95% indican que son distintos de cero y toman valores muy cercanos a los reales. Todos los coeficientes autorregresivos estimados son significativos según indica el IC al 95% y se ajustan muy bien a la relación verdadera entre las variables que se han usado para simular los datos. 

En cuanto a la estimación de la matriz de varianza-covarianza, las estimaciones para $\sigma^2_{x1}$ y $\sigma^2_{x2}$ contienen al 1 en su IC al 95% y su valor estimado es muy cercano. En cuanto a las covarianzas, los intervalos de confianza incluyen al 0 al 95% por lo que recogen el verdadero valor.

Por último, vamos a estimar estos modelos para los datos macroeconómicos de LIBOR y OVERNIGHT. Después de estimar varios modelos con diferentes números de retardos, el modelo que menor AIC tiene incorpora 5 retardos y sin términos deterministas tanto en la parte de largo plazo con en la de corto plazo.
"

libor_overnight <- readxl::read_excel("libor_overnight.xlsx")
libor_overnight <- ts(libor_overnight)

vecm_liborover = VECM(libor_overnight, lag = 5, include = "none", estim = "ML")
summary(vecm_liborover)

"
Por otro lado, la estimación de manera bayesiana con 5 retardos y constante en la relación de largo plazo es la siguiente:
 " 

data <- gen_vec(libor_overnight, p = 6, r=1, iterations = 5000, burnin = 2500)
data <- add_priors(data,
                   coint = list(v_i = 0, p_tau_i = 1),
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(df = 0, scale = .0001))

# Reset random number generator for reproducibility
set.seed(100)

# Obtain data matrices
y <- t(data$data$Y)
w <- t(data$data$W)
x <- t(data$data$X)

r <- data$model$rank # Set rank

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
k_gamma <- k * k_x # Total number of non-cointegration coefficients

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta

# Priors
a_mu_prior <- data$priors$noncointegration$mu # Prior means
a_v_i_prior <- data$priors$noncointegration$v_i # Inverse of the prior covariance matrix

v_i <- data$priors$cointegration$v_i
p_tau_i <- data$priors$cointegration$p_tau_i

sigma_df_prior <- data$priors$sigma$df # Prior degrees of freedom
sigma_scale_prior <- data$priors$sigma$scale # Prior covariance matrix
sigma_df_post <- tt + sigma_df_prior # Posterior degrees of freedom

# Initial values
beta <- matrix(0, k_w, r)
beta[1:r, 1:r] <- diag(1, r)

sigma_i <- diag(1 / .0001, k)

g_i <- sigma_i

iterations <- data$model$iterations # Number of iterations of the Gibbs sampler
burnin <- data$model$burnin # Number of burn-in draws
draws <- iterations + burnin # Total number of draws

# Data containers
draws_alpha <- matrix(NA, k_alpha, iterations)
draws_beta <- matrix(NA, k_beta, iterations)
draws_pi <- matrix(NA, k * k_w, iterations)
draws_gamma <- matrix(NA, k_gamma, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = sigma_i,
                         v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                         gamma_mu_prior = a_mu_prior,
                         gamma_v_i_prior = a_v_i_prior)
  alpha <- temp$alpha
  beta <- temp$beta
  Pi <- temp$Pi
  gamma <- temp$Gamma
  
  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  sigma_scale_post <- solve(tcrossprod(u) + v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_post)[,, 1], k)
  sigma <- solve(sigma_i)
  
  # Update g_i
  g_i <- sigma_i
  
  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_pi[, draw - burnin] <- Pi
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- sigma
  }
}

#Sacamos la estimación de los coeficientes de cointegración
beta <- apply(t(draws_beta) / t(draws_beta)[, 1], 2, mean) # Obtain means for every row
beta <- matrix(beta, k_w) # Transform mean vector into a matrix
beta <- round(beta, 3) # Round values
dimnames(beta) <- list(dimnames(w)[[1]], NULL) # Rename matrix dimensions

beta

# Number of non-deterministic coefficients
k_nondet <- (k_x - 0) * k

# Generate bvec object
bvec_est_liborover <- bvec(y = data$data$Y,
                           w = data$data$W,
                           x = data$data$X[, 1:10],
                           x_d = data$data$X[, -(1:10)],
                           Pi = draws_pi,
                           r = 1,
                           Gamma = draws_gamma[1:k_nondet,],
                           Sigma = draws_sigma)

summary(bvec_est_liborover)

"
Se puede ver como las estimaciones con el método bayesiano son muy parecidas a las estimaciones frecuentistas.
"
## Análisis de la capacidad predictiva
"
Por último, se presentan los resultados de la capacidad predictiva de cada método. Para analizar la capacidad predictiva, se va a realizar un walk-forward validation. Se parte de un modelo inicial estimado con un número (normalmente pequeño) de observaciones $t$ y se predice el valor en $t+1$. Una vez estimado el valor, se calcula el error cuadrático de la estimación y esa observación pasa a a formar parte de los datos con los que se estimará otra vez el modelo. Se repite este proceso hasta haber recorrido todos los datos que se tengan. Al final se obtiene el error de predicción como la media de los errores cuadráticos de cada observación predecida. 

Como hemos visto en el apartado anterior, estimar los modelos de forma bayesiana mediante muestreadores es computacionalmente mucho más costoso que de forma frecuentista (más de 100 veces más costoso). Por esto, ha sido necesario modificar un poco el algoritmo de walk-forward validation para que sea viable. 

Tanto para los modelos BVAR como para los BVECM, la estimación del modelo se ha de hacer con 50000 iterations y 25000 burnin para asegurar la convergencia. Cada modelo que se estima tarda entre 45 segundos y 1 minuto 10 segundos. Además, la predicción de observaciones también es más lenta que para los modelos frecuentistas por lo que en general, estamos hablando que estimar un modelo y hacer una predicción tarda más o menos 1 minuto y medio. Hacer esto 1000 veces para 4 conjuntos de datos diferentes supondría 100h de cómputo. 

La estrategia que se ha decidido tomar para reducir el tiempo de cómputo ha sido la siguiente: El modelo inicial contendrá un número de observaciones iniciales tal que el número de iteraciones del walk-forward validation sea de 30 bajo a condición de que si las series son diarias se predecirá a un horizonte temporal de $t+30$ y si son cuatrimestrales de $t+4$ en vez de ser $t+1$ el horizonte temporal.

Esta decisión se ha tomado para que haya 30 estimaciones de mrse y se pueda aproximar a una distribución normal para sacar intervalos de confianza, y para que el tiempo de cómputo por cada conjunto de datos no sea superior a 1 hora (4 horas en total). Además, tiene sentido usar un horizonte temporal de $t+30$ para series diarias y $t+4$ para las cuatrimestrales por que lo que se está haciendo es estimar a un mes vista para las series diarias, y a un año vista en las series cuatrimestrales.
"

# Creamos las muestras iniciales de entrenamiento para el walk forward validation

train_varsim <- 1:round(dim(datos_var)[1]*0.1, 0)
train_usmacro <- 1:75
train_vecsim <- 1:round(dim(datos_VEC)[1]*0.1, 0)
train_liborover <- 1:1238

datosvar_train <- datos_var[train_varsim, ]
ahead_datosvar <- dim(datos_var)[1]-round(dim(datos_var)[1]*0.1, 0)

us_train <- varmat_us[train_usmacro, ]
ahead_us <- dim(varmat_us)[1]-75

datosvec_train <- datos_VEC[train_vecsim, ]
ahead_datosvec <- dim(datos_VEC)[1]-round(dim(datos_VEC)[1]*0.1, 0)

liborover_train <- libor_overnight[train_liborover, ]
ahead_liborover <- dim(libor_overnight)[1]-1238

# código para hacer el walk forward validation. Se usará n = 30 días para las diarias y n = 4 para las cuatrimestrales

############## VAR con datos simulados
# Iniciar msre

msre_x1_var <- c()
msre_x2_var <- c()

### BUCLE
for(i in 1:(round((dim(datos_var)[1] - 100) / 30, 0)-1)){
  # estimar modelo
  varfit_train <- VAR(datosvar_train, p=3, type = "none")
  
  # predecir siguientes n observaciones
  preds_varsim <- predict(varfit_train, n.ahead=30)
  
  # sacar msre
  msre_x1_var =  c(msre_x1_var, sqrt(sum((preds_varsim$fcst$y1[,1]-datos_var[(round(dim(datos_var)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_var)[1]*0.1, 0)+30*i),1])^2)/length(preds_varsim$fcst$y1[,1])))
  msre_x2_var = c(msre_x2_var, sqrt(sum((preds_varsim$fcst$y2[,1]-datos_var[(round(dim(datos_var)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_var)[1]*0.1, 0)+30*i),2])^2)/length(preds_varsim$fcst$y2[,1])))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_varsim <- 1:(round(dim(datos_var)[1]*0.1, 0) + 30*i)
  datosvar_train <- datos_var[train_varsim, ]
  
}


############### BVAR con datos simulados
# Iniciar msre

msre_x1_bvar <- c()
msre_x2_bvar <- c()

train_varsim <- 1:round(dim(datos_var)[1]*0.1, 0)
datosvar_train <- datos_var[train_varsim, ]

### BUCLE
for(i in 1:(round((dim(datos_var)[1] - 100) / 30, 0)-1)){
  # estimar modelo
  run_train <- bvar(datosvar_train, lags = 3, 
            n_draw = 50000, n_burn = 25000,
            priors = priors, mh = mh)
  
  # predecir siguientes n observaciones
  preds_bvarsim <- predict(run_train, horizon=30)
  
  # sacar msre
  est_bvarsim_x1 <- colMeans(preds_bvarsim$fcast[,,1])
  est_bvarsim_x2 <- colMeans(preds_bvarsim$fcast[,,2])
  msre_x1_bvar = c(msre_x1_bvar, sqrt(sum((est_bvarsim_x1-datos_var[(round(dim(datos_var)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_var)[1]*0.1, 0)+30*i),1])^2)/length(est_bvarsim_x1)))
  msre_x2_bvar = c(msre_x2_bvar, sqrt(sum((est_bvarsim_x2-datos_var[(round(dim(datos_var)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_var)[1]*0.1, 0)+30*i),2])^2)/length(est_bvarsim_x2)))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_varsim <- 1:(round(dim(datos_var)[1]*0.1, 0) + 30*i)
  datosvar_train <- datos_var[train_varsim, ]
  
}


# código para hacer el walk forward validation. Se usará n = 30 días para las diarias y n = 4 para las cuatrimestrales

############## VAR con USA macro data (cuatrimestral)
# Iniciar msre

msre_x1_var_us <- c()
msre_x2_var_us <- c()

### BUCLE
for(i in 1:(round((dim(varmat_us)[1] - 75) / 4, 0)-1)){
  # estimar modelo
  varfit_train_us <- VAR(us_train, p=2, type = "none")
  
  # predecir siguientes n observaciones
  preds_usvar <- predict(varfit_train_us, n.ahead=4)
  
  # sacar msre
  msre_x1_var_us = c(msre_x1_var_us, sqrt(sum((preds_usvar$fcst$Inflación[,1]-varmat_us[(75+1+(4*(i-1))):(75+4*i),1])^2)/length(preds_usvar$fcst$Inflación[,1])))
  msre_x2_var_us = c(msre_x2_var_us, sqrt(sum((preds_usvar$fcst$Desempleo[,1]-varmat_us[(75+1+(4*(i-1))):(75+4*i),2])^2)/length(preds_usvar$fcst$Desempleo[,1])))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_usmacro <- 1:(75 + 4*i)
  us_train <- varmat_us[train_usmacro, ]
  
}


############### BVAR con USA macro data (cuatrimestral)
# Iniciar msre

msre_x1_bvar_us <- c()
msre_x2_bvar_us <- c()

train_usmacro <- 1:75
us_train <- varmat_us[train_usmacro, ]
### BUCLE
for(i in 1:(round((dim(varmat_us)[1] - 75) / 4, 0)-1)){
  # estimar modelo
  run_train_us <- bvar(us_train, lags = 2, 
            n_draw = 50000, n_burn = 25000,
            priors = priors, mh = mh)
  
  # predecir siguientes n observaciones
  preds_bvarus <- predict(run_train_us, horizon=4)
  
  # sacar msre
  est_bvarsim_x1_us <- colMeans(preds_bvarus$fcast[,,1])
  est_bvarsim_x2_us <- colMeans(preds_bvarus$fcast[,,2])
  msre_x1_bvar_us = c(msre_x1_bvar_us, sqrt(sum((est_bvarsim_x1_us-varmat_us[(75+1+(4*(i-1))):(75+4*i),1])^2)/length(est_bvarsim_x1_us)))
  msre_x2_bvar_us = c(msre_x2_bvar_us, sqrt(sum((est_bvarsim_x2_us-varmat_us[(75+1+(4*(i-1))):(75+4*i),2])^2)/length(est_bvarsim_x2_us)))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_usmacro <- 1:(75 + 4*i)
  us_train <- varmat_us[train_usmacro, ]
  
}


# código para hacer el walk forward validation. Se usará n = 30 días para las diarias y n = 4 para las cuatrimestrales

############## VECM con datos simulados
# Iniciar msre

msre_x1_vecm <- c()
msre_x2_vecm <- c()

### BUCLE
for(i in 1:(round((dim(datos_VEC)[1] - 100) / 30, 0)-1)){
  # estimar modelo
  vecmfit_train <- VECM(datosvec_train, lag = 2, include = "none", estim = "ML")
  
  # predecir siguientes n observaciones
  preds_vecmsim <- predict(vecmfit_train, n.ahead=30)
  
  # sacar msre
  msre_x1_vecm =  c(msre_x1_vecm, sqrt(sum((preds_vecmsim[,1]-datos_VEC[(round(dim(datos_VEC)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_VEC)[1]*0.1, 0)+30*i),1])^2)/length(preds_vecmsim[,1])))
  msre_x2_vecm = c(msre_x2_vecm, sqrt(sum((preds_vecmsim[,2]-datos_VEC[(round(dim(datos_VEC)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_VEC)[1]*0.1, 0)+30*i),2])^2)/length(preds_vecmsim[,2])))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_vecsim <- 1:(round(dim(datos_VEC)[1]*0.1, 0) + 30*i)
  datosvec_train <- datos_VEC[train_vecsim, ]
  
}


############### BVECM con datos simulados
# Iniciar msre

msre_x1_bvecm <- c()
msre_x2_bvecm <- c()

train_vecsim <- 1:round(dim(datos_VEC)[1]*0.1, 0)
datosvec_train <- datos_VEC[train_vecsim, ]
### BUCLE
for(i in 1:(round((dim(datos_VEC)[1] - 100) / 30, 0)-1)){
  # estimar modelo
  data <- gen_vec(ts(datosvec_train), p = 3, r=1, iterations = 5000, burnin = 2500)
  data <- add_priors(data,
                   coint = list(v_i = 0, p_tau_i = 1),
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(df = 0, scale = .0001))
  
  set.seed(100)

# Obtain data matrices
y <- t(data$data$Y)
w <- t(data$data$W)
x <- t(data$data$X)

r <- data$model$rank # Set rank

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
k_gamma <- k * k_x # Total number of non-cointegration coefficients

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta

# Priors
a_mu_prior <- data$priors$noncointegration$mu # Prior means
a_v_i_prior <- data$priors$noncointegration$v_i # Inverse of the prior covariance matrix

v_i <- data$priors$cointegration$v_i
p_tau_i <- data$priors$cointegration$p_tau_i

sigma_df_prior <- data$priors$sigma$df # Prior degrees of freedom
sigma_scale_prior <- data$priors$sigma$scale # Prior covariance matrix
sigma_df_post <- tt + sigma_df_prior # Posterior degrees of freedom

# Initial values
beta <- matrix(0, k_w, r)
beta[1:r, 1:r] <- diag(1, r)

sigma_i <- diag(1 / .0001, k)

g_i <- sigma_i

iterations <- data$model$iterations # Number of iterations of the Gibbs sampler
burnin <- data$model$burnin # Number of burn-in draws
draws <- iterations + burnin # Total number of draws

# Data containers
draws_alpha <- matrix(NA, k_alpha, iterations)
draws_beta <- matrix(NA, k_beta, iterations)
draws_pi <- matrix(NA, k * k_w, iterations)
draws_gamma <- matrix(NA, k_gamma, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = sigma_i,
                           v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                           gamma_mu_prior = a_mu_prior,
                           gamma_v_i_prior = a_v_i_prior)
  alpha <- temp$alpha
  beta <- temp$beta
  Pi <- temp$Pi
  gamma <- temp$Gamma
  
  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  sigma_scale_post <- solve(tcrossprod(u) + v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_post)[,, 1], k)
  sigma <- solve(sigma_i)
  
  # Update g_i
  g_i <- sigma_i
  
  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_pi[, draw - burnin] <- Pi
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- sigma
  }
}

k_nondet <- (k_x - 0) * k

# Generate bvec object
bvec_est_sim <- bvec(y = data$data$Y,
                 w = data$data$W,
                 x = data$data$X[, 1:4],
                 x_d = data$data$X[, -(1:4)],
                 Pi = draws_pi,
                 r = 1,
                 Gamma = draws_gamma[1:k_nondet,],
                 Sigma = draws_sigma)

bvar_form_sim <- bvec_to_bvar(bvec_est_sim)


  # predecir siguientes n observaciones
  preds_bvarsim <- predict.bvar(bvar_form_sim, n.ahead=30)
  
  # sacar msre
  est_bvecsim_x1 <- preds_bvarsim$fcst$`Series 1`[, 2]
  est_bvecsim_x2 <- preds_bvarsim$fcst$`Series 2`[, 2]
  msre_x1_bvecm = c(msre_x1_bvecm, sqrt(sum((est_bvecsim_x1-datos_VEC[(round(dim(datos_VEC)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_VEC)[1]*0.1, 0)+30*i),1])^2)/length(est_bvecsim_x1)))
  msre_x2_bvecm = c(msre_x2_bvecm, sqrt(sum((est_bvecsim_x2-datos_VEC[(round(dim(datos_VEC)[1]*0.1, 0)+1+(30*(i-1))):(round(dim(datos_VEC)[1]*0.1, 0)+30*i),2])^2)/length(est_bvecsim_x2)))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_vecsim <- 1:(round(dim(datos_VEC)[1]*0.1, 0) + 30*i)
  datosvec_train <- datos_VEC[train_vecsim, ]
  
}


# código para hacer el walk forward validation. Se usará n = 30 días para las diarias y n = 4 para las cuatrimestrales

############## VECM con datos de LIBOR y OVERNIGHT
# Iniciar msre

msre_x1_vecm_liborover <- c()
msre_x2_vecm_liborover <- c()

### BUCLE
for(i in 1:(round((dim(libor_overnight)[1] - 1238) / 30, 0)-1)){
  # estimar modelo
  vecm_liborover_train = VECM(liborover_train, lag = 5, include = "none", estim = "ML", LRinclude = "both")
  
  # predecir siguientes n observaciones
  preds_liborover <- predict(vecm_liborover_train, n.ahead=30)
  
  # sacar msre
  msre_x1_vecm_liborover =  c(msre_x1_vecm_liborover, sqrt(sum((preds_liborover[,1]-libor_overnight[(1238+1+(30*(i-1))):(1238+30*i),1])^2)/length(preds_liborover[,1])))
  msre_x2_vecm_liborover = c(msre_x2_vecm_liborover, sqrt(sum((preds_liborover[,2]-libor_overnight[(1238+1+(30*(i-1))):(1238+30*i),2])^2)/length(preds_liborover[,2])))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_liborover <- 1:(1238 + 30*i)
  liborover_train <- libor_overnight[train_liborover, ]
  
}



############### BVECM con datos de LIBOR y OVERNIGHT
# Iniciar msre

msre_x1_bvecm_liborover <- c()
msre_x2_bvecm_liborover <- c()

train_liborover <- 1:1238
liborover_train <- libor_overnight[train_liborover, ]
### BUCLE
for(i in 1:(round((dim(libor_overnight)[1] - 1238) / 30, 0)-1)){
  # estimar modelo
  data <- gen_vec(ts(liborover_train), p = 6, r=1, iterations = 5000, burnin = 2500)
  data <- add_priors(data,
                   coint = list(v_i = 0, p_tau_i = 1),
                   coef = list(v_i = 0, v_i_det = 0),
                   sigma = list(df = 0, scale = .0001))
  
set.seed(100)

# Obtain data matrices
y <- t(data$data$Y)
w <- t(data$data$W)
x <- t(data$data$X)

r <- data$model$rank # Set rank

tt <- ncol(y) # Number of observations
k <- nrow(y) # Number of endogenous variables
k_w <- nrow(w) # Number of regressors in error correction term
k_x <- nrow(x) # Number of differenced regressors and unrestrictec deterministic terms
k_gamma <- k * k_x # Total number of non-cointegration coefficients

k_alpha <- k * r # Number of elements in alpha
k_beta <- k_w * r # Number of elements in beta

# Priors
a_mu_prior <- data$priors$noncointegration$mu # Prior means
a_v_i_prior <- data$priors$noncointegration$v_i # Inverse of the prior covariance matrix

v_i <- data$priors$cointegration$v_i
p_tau_i <- data$priors$cointegration$p_tau_i

sigma_df_prior <- data$priors$sigma$df # Prior degrees of freedom
sigma_scale_prior <- data$priors$sigma$scale # Prior covariance matrix
sigma_df_post <- tt + sigma_df_prior # Posterior degrees of freedom

# Initial values
beta <- matrix(0, k_w, r)
beta[1:r, 1:r] <- diag(1, r)

sigma_i <- diag(1 / .0001, k)

g_i <- sigma_i

iterations <- data$model$iterations # Number of iterations of the Gibbs sampler
burnin <- data$model$burnin # Number of burn-in draws
draws <- iterations + burnin # Total number of draws

# Data containers
draws_alpha <- matrix(NA, k_alpha, iterations)
draws_beta <- matrix(NA, k_beta, iterations)
draws_pi <- matrix(NA, k * k_w, iterations)
draws_gamma <- matrix(NA, k_gamma, iterations)
draws_sigma <- matrix(NA, k^2, iterations)

# Start Gibbs sampler
for (draw in 1:draws) {
  
  # Draw conditional mean parameters
  temp <- post_coint_kls(y = y, beta = beta, w = w, x = x, sigma_i = sigma_i,
                           v_i = v_i, p_tau_i = p_tau_i, g_i = g_i,
                           gamma_mu_prior = a_mu_prior,
                           gamma_v_i_prior = a_v_i_prior)
  alpha <- temp$alpha
  beta <- temp$beta
  Pi <- temp$Pi
  gamma <- temp$Gamma
  
  # Draw variance-covariance matrix
  u <- y - Pi %*% w - matrix(gamma, k) %*% x
  sigma_scale_post <- solve(tcrossprod(u) + v_i * alpha %*% tcrossprod(crossprod(beta, p_tau_i) %*% beta, alpha))
  sigma_i <- matrix(rWishart(1, sigma_df_post, sigma_scale_post)[,, 1], k)
  sigma <- solve(sigma_i)
  
  # Update g_i
  g_i <- sigma_i
  
  # Store draws
  if (draw > burnin) {
    draws_alpha[, draw - burnin] <- alpha
    draws_beta[, draw - burnin] <- beta
    draws_pi[, draw - burnin] <- Pi
    draws_gamma[, draw - burnin] <- gamma
    draws_sigma[, draw - burnin] <- sigma
  }
}

k_nondet <- (k_x - 0) * k

# Generate bvec object
bvec_est_liborover_train <- bvec(y = data$data$Y,
                 w = data$data$W,
                 x = data$data$X[, 1:10],
                 x_d = data$data$X[, -(1:10)],
                 Pi = draws_pi,
                 r = 1,
                 Gamma = draws_gamma[1:k_nondet,],
                 Sigma = draws_sigma)

bvar_form_liborover <- bvec_to_bvar(bvec_est_liborover_train)


  # predecir siguientes n observaciones
  preds_bvarliborover <- predict.bvar(bvar_form_liborover, n.ahead=30)
  
  # sacar msre
  est_bveclibor_x1 <- preds_bvarliborover$fcst$libor[, 2]
  est_bveclibor_x2 <- preds_bvarliborover$fcst$overnight[, 2]
  msre_x1_bvecm_liborover = c(msre_x1_bvecm_liborover, sqrt(sum((est_bveclibor_x1-libor_overnight[(1238+1+(30*(i-1))):(1238+30*i),1])^2)/length(est_bveclibor_x1)))
  msre_x2_bvecm_liborover = c(msre_x2_bvecm_liborover, sqrt(sum((est_bveclibor_x2-libor_overnight[(1238+1+(30*(i-1))):(1238+30*i),2])^2)/length(est_bveclibor_x2)))
  
  # Añadir nuevas observaciones a los datos de entrenamiento
  train_liborover <- 1:(1238 + 30*i)
  liborover_train <- libor_overnight[train_liborover, ]
  
}

