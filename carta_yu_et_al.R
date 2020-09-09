######################################################
### Carta propuesta por Yu et al 2016 (VSSI sintética)
######################################################

## Verificamos si los paquetes estan disponibles
if(!require('ggplot2')) install.packages('ggplot2')
require('ggplot2')

## Se carga el algoritmo: genetic algorithm
source('C:/Users/Carmen C/Documents/R/Proyecto-Maestria/ag.R')

# k1 = Límite de control
# k2 = Límite de advertencia
# hL = Intervalo de muestreo relajado
# hs = Intervalo de muestreo estricto
# ns = Tamaño de muestra relajado
# nL = Tamaño de muestra estricto
# m = Límite de control inferior de la CRL
# nmax = Tamaño de la muestra máximo

##########################
#k1<-7.1368
#k2<-1.9826
#hL<-6.5031
#hs<-3.3638
#ns<-12
#nL<-14
#m<-9
##########################

##Se construye la función objetivo del modelo matemático: ATS1

#Cambios de la media a evaluar
delta<-0.1
#delta<-c(0.1, 0.3, 0.5, 0.6, 0.8, 1.0, 1.5)

nmax<-14
miu<-0
des<-1

#deseado <- 200
#a<-1:600
#index <-1
#while(170>deseado || deseado>196 || (results$mejor_individuo_global[2]-results$mejor_individuo_global[1])>3){

#set.seed(a[index])
#index=index+1
set.seed(2)

funcion <- function(ns, nL, m, k1, k2, hL, hs){
  
  #Iniciación de variables
  P<-array(0, c(m+2,m+2))
  Q<-array(0, c(m+1,m+1))
  q<-matrix(0, nrow=m+1, ncol=1)
  q[1]<-1
  one<-array(1,c(m+1,1))
  I<-array(diag(m+1), c(m+1,m+1))
  pil<-array(0,c(1,m-1))
  tao<-370.4
  nmax<-14
  
  pc<-runif(1)
  
  ##Se calcula el ATS0 para construir la restricción ATS0 >= 370.4
  
  #Probabilidades con delta = 0
  p010<-1-pnorm(k1)+pnorm(-k1)
  p011<-1-pnorm(k1)+pnorm(-k1)
  p020<-pnorm(k1)-pnorm(k2)+pnorm(-k2)-pnorm(-k1)
  p021<-pnorm(k1)-pnorm(k2)+pnorm(-k2)-pnorm(-k1)
  p030<-pnorm(k2)-pnorm(-k2)
  p031<-pnorm(k2)-pnorm(-k2)
  
  #Matriz de transición P0
  P0<-array(0,c(m+2,m+2))
  P0[m+1,1]<-p020
  P0[m+1,m+1]<-p030
  P0[m+1,m+2]<-p010
  P0[m+2,m+2]<-1
  P0[1,2]<-p031
  P0[1,m+2]<-p021+p011
  for (l in 2:m){
    P0[l,m+2]<-p020+p010
    P0[l,l+1]<-p030
  }
  
  #Cálculo de ARL0
  Q0<-P0[0:m+1,0:m+1]
  I0<-diag(m+1)
  inv0<-(I0-Q0)
  ARL0<-t(q)%*%solve(inv0,tol = 1e-40)%*%one
  
  #Probabilidades de estado-estable con delta = 0
  pi00<-(1-p030)/(1-p030+p031)
  pi0l<-array(0,c(1,m-1))
  for (l in 1:m-1){
    pi0l[1,l]<-(p031*(p030^(l-1))*(1-p030))/(1-p030+p031)
  }
  pi0m<-(p031*((p030)^(m-1)))/(1-p030+p031)
  
  #Promedio esperado de Eh0
  Eh0<-pi00*hs+(sum(pi0l)+pi0m)*hL
  
  #Promedio esperado de En0
  En0<-as.integer(pi00*nL+(sum(pi0l)+pi0m)*ns)
  
  #Cálculo de ATS0
  ATS0<-ARL0*Eh0
  
  if (ATS0 >= tao){
    #Probabilidades con delta > 0
    p10<-1-pnorm(k1-(delta*(sqrt(ns))))+pnorm(-k1-(delta*(sqrt(ns))))
    p11<-1-pnorm(k1-(delta*(sqrt(nL))))+pnorm(-k1-(delta*(sqrt(nL))))
    p20<-pnorm(k1-(delta*(sqrt(ns))))-pnorm(k2-(delta*(sqrt(ns))))
    +pnorm(-k2-(delta*(sqrt(ns))))-pnorm(-k1-(delta*(sqrt(ns))))
    p21<-pnorm(k1-(delta*(sqrt(nL))))-pnorm(k2-(delta*(sqrt(nL))))
    +pnorm(-k2-(delta*(sqrt(nL))))-pnorm(-k1-(delta*(sqrt(nL))))
    p30<-pnorm(k2-(delta*(sqrt(ns))))-pnorm(-k2-(delta*(sqrt(ns))))
    p31<-pnorm(k2-(delta*(sqrt(nL))))-pnorm(-k2-(delta*(sqrt(nL))))
    
    #Matriz de probabilidades de transición P
    P[m+1,1]<-p20
    P[m+1,m+1]<-p30
    P[m+1,m+2]<-p10
    P[m+2,m+2]<-1
    P[1,2]<-p31
    P[1,m+2]<-p21+p11
    for (l in 2:m){
      P[l,m+2]<-p20+p10
      P[l,l+1]<-p30
    }
    
    #Cálculo de ARL1
    Q<-P[0:m+1,0:m+1]
    inv<-(I-Q)
    ARL1<-t(q)%*%solve(inv,tol = 1e-40)%*%one
    
    #Probabilidades de estado-estable con delta > 0
    pi0<-(1-p30)/(1-p30+p31)
    for (l in 1:m-1){
      pil[1,l]<-(p31*((p30)^(l-1))*(1-p30))/(1-p30+p31)
    }
    pim<-(p31*((p30)^(m-1)))/(1-p30+p31)
    
    #Promedio esperado de Eh[delta] y En[delta]
    Ehdelta<-pi0*hs+(sum(pil)+pim)*hL
    Endelta<-as.integer(pi0*nL+(sum(pil)+pim)*ns)
    
    #Cálculo de ATS1
    ATS1<-ARL1*Ehdelta
    
    #Promedio esperado total del tamaño de muestra En
    En<-pc*En0+(1-pc)*Endelta
    
    if (En<=nmax){
      return(round(ATS1,3) )
    }else{
      return(100000000)}
  }else{
    return(1000000000)}
}

#funcion(ns, nL, m, k1, k2, hL, hs)

##Aplicación del algoritmo genético a la función objetivo

n_poblacion <- 100
n_variables <- 7
limite_inf <- c(as.integer(1),as.integer(1),as.integer(1),0.1,0.1,0.1,0.1)
limite_sup <- c(as.integer(14),as.integer(14),as.integer(30),10,5,10,5)
n_generaciones <- 50

results <- optimizar_ga(funcion_objetivo = funcion, n_variables = n_variables, nmax, miu, des, 
                        optimizacion = "minimizar", limite_inf = limite_inf,
                        limite_sup = limite_sup, n_poblacion = n_poblacion,  
                        n_generaciones = n_generaciones, distribucion = "aleatoria",
                        verbose = 1, metodo_seleccion = "tournament")

head(results$df_resultados)
results$mejor_individuo_global
results$mejor_valor_global
#deseado<-results$mejor_valor_global
ggplot(data = results$df_resultados,
       aes(x = generacion, y = fitness)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  labs(title = "Evolución del fitness a lo largo de las generaciones") + 
  theme_bw()

#print(a[index])