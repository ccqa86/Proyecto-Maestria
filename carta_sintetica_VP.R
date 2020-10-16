#########################################################
### Carta sintética VP para procesos autocorrelacionados
#########################################################

## Verificamos si los paquetes estan disponibles
if(!require('ggplot2')) install.packages('ggplot2')
require('ggplot2')

library("writexl")

## Se carga el algoritmo: genetic algorithm
source('C:/Users/Carmen C/Documents/R/Proyecto-Maestria/ag.R')

# k1 = Límite de control relajado
# k2 = Límite de control estricto
# w1 = Límite de advertencia relajado
# w2 = Límite de advertencia estricto
# h1 = Intervalo de muestreo relajado
# h2 = Intervalo de muestreo estricto
# n1 = Tamaño de muestra relajado
# n2 = Tamaño de muestra estricto
# m = Límite de control de la CRL
# nmax = Tamaño de la muestra máximo

##########################
#k1<-3.196
#k2<-2.078
#w1<-2.385
#w2<-0.437
#h1<-14.061
#h2<-6.100
#n1<-13
#n2<-14
#m<-17
##########################

##Se construye la función objetivo del modelo matemático: ATS1
nmax<-15
miu<-0
des<-1

#deseado <- 200
#a<-1:600
#index <-1
#while(170>deseado || deseado>190 || (results$mejor_individuo_global[2]-results$mejor_individuo_global[1])>3){

#set.seed(a[index])
#index=index+1
set.seed(43)

#Niveles de psi a evaluar
psi_1<-c(0.1, 0.5)
#psi<-c(0.1, 0.5)

#Niveles de phi a evaluar
phi_1<-c(0.2, 0.4, 0.6, 0.8, 0.99)
#phi<-c(0.2, 0.4, 0.6, 0.8, 0.99)

#Cambios de la media a evaluar
delta_1<-c(0.1, 0.3, 0.5, 0.6, 0.8, 1.0)
#delta<-c(0.1, 0.3, 0.5, 0.6, 0.8, 1.0)

for (escenario in (1:6)){
 escenario_name=c("VP","VSSIWL","VSSI","VSI","VSS","X")
for (aa in 1:length(psi_1)){
  for (bb in 1:length(phi_1)){
    for (cc in 1:length(delta_1)){
      psi<-psi_1[aa]
      phi<-phi_1[bb]
      delta<-delta_1[cc]
      
      set.seed(43)

      #funcion_par <- function(psi_1=psi, phi_1=phi, delta_1=delta){
      #  return (list(psi_1, phi_1, delta_1))
      #}
      
      funcion <- function(k1, k2, w1, w2, h1, h2, n1, n2, m){
        #valores<-funcion_par()
        #psi=unlist(valores[1])
        #phi=unlist(valores[2])
        #delta=unlist(valores[3])

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
        rho_adv<-((phi^h2)*psi)/sqrt((psi+((1-psi)/n2))*(psi+((1-psi)/n1)))
        rho_cen<-((phi^h1)*psi)/sqrt((psi+((1-psi)/n1))*(psi+((1-psi)/n2)))
        rho1<-(n1*(phi^h1)*psi)/((n1*psi)+1-psi)
        rho2<-(n2*(phi^h2)*psi)/((n2*psi)+1-psi)
        
        pc<-runif(1)
        
        ##Se calcula el ATS0 para construir la restricción ATS0 >= 370.4
        
        #Probabilidades con delta = 0
        p010<-1-pnorm(k1/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)+pnorm(-k1/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)
        p011<-1-pnorm(k2/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)+pnorm(-k2/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)
        p020<-pnorm(k1/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)-pnorm(w1/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)+pnorm(-w1/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)-pnorm(-k1/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)
        p021<-pnorm(k2/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)-pnorm(w2/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)+pnorm(-w2/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)-pnorm(-k2/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)
        p030<-pnorm(w1/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)-pnorm(-w1/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)
        p031<-pnorm(w2/sqrt((1+(n2-1)*rho_adv)),mean=miu,sd=des)-pnorm(-w2/sqrt((1+(n2-1)*rho_adv)),mean=miu,sd=des)
        
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
        Eh0<-pi00*h2+(sum(pi0l)+pi0m)*h1
        
        #Promedio esperado de En0
        En0<-as.integer(pi00*n2+(sum(pi0l)+pi0m)*n1)
        
        #Cálculo de ATS0
        ATS0<-ARL0*Eh0
      
        if (ATS0 >= tao){
          #Probabilidades con delta > 0
          p10<-1-pnorm((k1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)+pnorm((-k1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)
          p11<-1-pnorm((k2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)+pnorm((-k2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)
          p20<-pnorm((k1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)-pnorm((w1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)+pnorm((-w1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)-pnorm((-k1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho_cen)),mean=miu,sd=des)
          p21<-pnorm((k2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)-pnorm((w2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)+pnorm((-w2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)-pnorm((-k2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho2)),mean=miu,sd=des)
          p30<-pnorm((w1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)-pnorm((-w1-(sqrt(n1)*delta))/sqrt((1+(n1-1)*rho1)),mean=miu,sd=des)
          p31<-pnorm((w2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho_adv)),mean=miu,sd=des)-pnorm((-w2-(sqrt(n2)*delta))/sqrt((1+(n2-1)*rho_adv)),mean=miu,sd=des)
          
          #Matriz de probabilidades de transición P
          P[m+1,1]<-p20
          P[m+1,m+1]<-p30
          P[m+1,m+2]<-p10
          P[m+2,m+2]<-1
          P[1,2]<-p31
          P[1,m+2]<-p21+p11
          for (p in 2:m){
            P[p,m+2]<-p20+p10
            P[p,p+1]<-p30
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
          Ehdelta<-pi0*h2+(sum(pil)+pim)*h1
          Endelta<-as.integer(pi0*n2+(sum(pil)+pim)*n1)
          
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

#funcion(k1, k2, w1, w2, h1, h2, n1, n2, m)

##Aplicación del algoritmo genético a la función objetivo

n_poblacion <- 300
n_variables <- 9
limite_inf <- c(miu+0.1,miu+0.1,miu+0.1,miu+0.1,miu+0.1,miu+0.1,as.integer(1),as.integer(1),as.integer(1))
limite_sup <- c(miu+4*des,(miu+4*des)*0.8,miu+4*des,(miu+4*des)*0.8,20,10,as.integer(nmax),as.integer(nmax),as.integer(30))
n_generaciones <- 75

results <- optimizar_ga(funcion_objetivo = funcion, n_variables = n_variables, nmax, miu, escenario, des, 
                        optimizacion = "minimizar", limite_inf = limite_inf,
                        limite_sup = limite_sup, n_poblacion = n_poblacion,  
                        n_generaciones = n_generaciones, distribucion = "aleatoria", metodo_seleccion = "tournament", verbose=5)

#print("-------------------------------------------------")
#head(results$df_resultados)
results$mejor_individuo_global
#results$mejor_valor_global

par(mfrow = c(1, 1))
print(ggplot(data = results$df_resultados,
       aes(x = generacion, y = fitness)) +
  geom_line(aes(group = 1)) +
  geom_point() + theme_bw() 
  + labs(title = paste("Evolución del fitness a lo largo de las generaciones para psi: ",
                       psi," phi: ",phi," delta: ", delta, " carta ", escenario_name[escenario], 
                       axis.title=element_text(size=10,face="bold"))))
    
#name=paste0("C:/Users/Carmen C/Documents/R/Proyecto-Maestria/Graficas/psi",psi*10,"phi",phi*100,"delta",delta*10,"carta",escenario_name[escenario],".jpeg")
#jpeg(filename = name)
if (aa==1 && bb==1 && cc==1 && escenario==1){
tabla=data.frame(carta=escenario_name[escenario], psi=psi, phi=phi, delta=delta, 
                 h1=unlist(results$mejor_individuo_global[5]), 
                 h2=unlist(results$mejor_individuo_global[6]), 
                 n1=unlist(results$mejor_individuo_global[7]), 
                 n2=unlist(results$mejor_individuo_global[8]), 
                 k1=unlist(results$mejor_individuo_global[1]),
                 k2=unlist(results$mejor_individuo_global[2]),
                 w1=unlist(results$mejor_individuo_global[3]), 
                 w2=unlist(results$mejor_individuo_global[4]),
                 m=unlist(results$mejor_individuo_global[9]), 
                 ATS1=unlist(results$mejor_valor_global))
}else {
  tabla=rbind(tabla,data.frame(carta=escenario_name[escenario], psi=psi, phi=phi, delta=delta, 
                               h1=unlist(results$mejor_individuo_global[5]), 
                               h2=unlist(results$mejor_individuo_global[6]), 
                               n1=unlist(results$mejor_individuo_global[7]), 
                               n2=unlist(results$mejor_individuo_global[8]), 
                               k1=unlist(results$mejor_individuo_global[1]),
                               k2=unlist(results$mejor_individuo_global[2]),
                               w1=unlist(results$mejor_individuo_global[3]), 
                               w2=unlist(results$mejor_individuo_global[4]),
                               m=unlist(results$mejor_individuo_global[9]), 
                               ATS1=unlist(results$mejor_valor_global)))
}
    }
  }
}
 #print(tabla)
 #write_xlsx(tabla,paste0("C:/Users/Carmen C/Documents/R/Proyecto-Maestria/Tablas/carta",escenario_name[escenario],".xlsx"))
}
write_xlsx(tabla,paste0("C:/Users/Carmen C/Documents/R/Proyecto-Maestria/Tablas/resultados.xlsx"))
#dev.off()
#deseado<-results$mejor_valor_global
#}
#print(a[index])