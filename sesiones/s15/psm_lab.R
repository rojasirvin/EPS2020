# Espacio de trabajo ----
rm(list = ls())
options(scipen=999) # Prevenir notación científica

setwd("C:/Users/Irvin/Documents/EPS2020/sesiones/s15")

library(pacman)
library(janitor)
library(sandwich)
library(AER)
library(estimatr)
library(clubSandwich)
library(data.table)
library(MatchIt)
library(Matching)

p_load(tidyverse, foreign, reshape2, psych, qwraps2, forcats, readxl, 
       broom, lmtest, margins, plm, rdrobust, multiwayvcov,
       wesanderson, sandwich, stargazer,
       readstata13, pscore, optmatch, kdensity, MatchIt, bootstrap, matlib, dplyr)

data.smoke <-read_csv("./cattaneo_smoking.csv")   %>% 
  clean_names() %>% 
  mutate(smoker=ifelse(mbsmoke=="smoker",1,0))




#Una primera aproximación

binaria <- "smoker"
variables1 <- c("mage", "medu", "nprenatal")
f1 <- as.formula(paste(binaria,
                 paste(variables1,
                       collapse ="+"),
                 sep= " ~ "))
print(f1)


#La fórmula está en f1
ps1 <- glm(formula=f1,
            family=binomial,
            data=data.smoke)

X  <- ps1$fitted
Y  <- data.smoke$bweight
Tr  <- data.smoke$smoker

rr  <- Matchby(Y=Y, Tr=Tr, X=X, M=1, by="")
summary(rr)


MatchBalance(f1,
             data = data.smoke,
             match.out = rr,
             ks = TRUE, #pedir la prueba de KS de distribuciones iguales
             nboots=500,
             digits=5,
             paired=TRUE,
             print.level=2) # 0=nada de output, 1=resumen, 2=todo








# Esta funciona

variables2 <- c("mmarried", "fbaby", "medu", "nprenatal", "foreign", "mhisp", "fage")

f2 <- as.formula(paste(binaria,
                       paste(variables2,
                             collapse ="+"),
                       sep= " ~ "))

# Estimamos el PS

ps2 <- glm(formula=f2,
           family=binomial,
           data=data.smoke)

X  <- ps2$fitted
Y  <- data.smoke$bweight
Tr  <- data.smoke$smoker

rr  <- Matchby(Y=Y, Tr=Tr, X=X, M=1, by="")
summary(rr)


MatchBalance(f2,
             data = data.smoke,
             match.out = rr,
             ks = TRUE, #pedir la prueba de KS de distribuciones iguales
             nboots=500,
             digits=5,
             paired=TRUE,
             print.level=2) # 0=nada de output, 1=resumen, 2=todo
















m.out <- matchit(smoker ~ mmarried + mage + mage_sq + fbaby + medu + nprenatal + foreign + mhisp + fage,
                 method = "nearest",
                 discard = "hull.control", #quedarnos con soporte común
                 data = data.smoke)











#1a. Wald 

mean_cliente<-data.morocco %>%
  group_by(treatment) %>% 
  summarize(p_cliente=mean(client, na.rm=F)) %>% 
  ungroup()


mean_gasto<-data.morocco %>%
  group_by(treatment) %>% 
  summarize(m_gasto=mean(expense_total, na.rm=F)) %>% 
  ungroup()

#Neceistamos la diferencia de gastos y de probabilidad de ser cliente
dif_gasto <- mean_gasto[2,2]-mean_gasto[1,2]
dif_cliente <- mean_cliente[2,2]-mean_cliente[1,2]

Wald <- as.numeric(dif_gasto / dif_cliente)
Wald


#1b Wald es el estimador de VI cuando el instrumento es dicotómico
Wald_vi <- ivreg(expense_total ~ client  | treatment,
                data=data.morocco)

#Notemos que obtenemos directamente el error estándar
summary(Wald_vi)



#1c Forma reducida

res_fr<- lm(expense_total ~ treatment +
              members_resid_bl + nadults_resid_bl +
              head_age_bl + act_livestock_bl + act_business_bl +
              borrowed_total_bl + members_resid_d_bl +
              nadults_resid_d_bl + head_age_d_bl + act_livestock_d_bl +
              act_business_d_bl + borrowed_total_d_bl +
              ccm_resp_activ + other_resp_activ + ccm_resp_activ_d + 
              other_resp_activ_d + factor(paire),
            na.action=na.omit,
            data=data.morocco)

coef_test(res_fr, vcov = "CR1S", 
          cluster = data.morocco$demi_paire, test = "naive-t")[1:2,]







#1d primera etapa

res_fs<- lm(client ~ treatment +
              members_resid_bl + nadults_resid_bl +
                head_age_bl + act_livestock_bl + act_business_bl +
                borrowed_total_bl + members_resid_d_bl +
                nadults_resid_d_bl + head_age_d_bl + act_livestock_d_bl +
                act_business_d_bl + borrowed_total_d_bl +
                ccm_resp_activ + other_resp_activ + ccm_resp_activ_d + 
                other_resp_activ_d + factor(paire),
              na.action=na.omit,
              data=data.morocco)

coef_test(res_fs, vcov = "CR1S", 
          cluster = data.morocco$demi_paire, test = "naive-t")[1:2,]






#1e Columna 3, Panel A de la Tabla 9 en Crepon et al.

res_mco <- lm(expense_total ~ client +
                members_resid_bl + nadults_resid_bl +
                head_age_bl + act_livestock_bl + act_business_bl +
                borrowed_total_bl + members_resid_d_bl +
                nadults_resid_d_bl + head_age_d_bl + act_livestock_d_bl +
                act_business_d_bl + borrowed_total_d_bl +
                ccm_resp_activ + other_resp_activ + ccm_resp_activ_d + 
                other_resp_activ_d + factor(paire),
              na.action=na.omit,
              data=filter(data.morocco,treatment==1))

coef_test(res_mco, vcov = "CR1S", 
          cluster = filter(data.morocco,treatment==1)$demi_paire, test = "naive-t")[1:2,]




#1g Columna 3, Panel B, Tabla 9. LATE:

res_iv <- ivreg(expense_total ~ client + members_resid_bl + nadults_resid_bl
     + head_age_bl + act_livestock_bl + act_business_bl 
     + borrowed_total_bl + members_resid_d_bl + nadults_resid_d_bl
     + head_age_d_bl + act_livestock_d_bl + act_business_d_bl 
     + borrowed_total_d_bl + ccm_resp_activ + other_resp_activ 
     + ccm_resp_activ_d  + other_resp_activ_d + factor(paire) |
       treatment +  members_resid_bl + nadults_resid_bl
     + head_age_bl + act_livestock_bl + act_business_bl 
     + borrowed_total_bl + members_resid_d_bl + nadults_resid_d_bl
     + head_age_d_bl + act_livestock_d_bl + act_business_d_bl 
     + borrowed_total_d_bl + ccm_resp_activ + other_resp_activ 
     + ccm_resp_activ_d  + other_resp_activ_d + factor(paire),
     data=data.morocco)

summary(res_iv)$coefficients[1:2,]

#Errores agrupados a nivel pareja (paire)
coef_test(res_iv, vcov = "CR1S", 
          cluster = data.morocco$demi_paire, test = "naive-t")[1:2,]






#Pregunta 2

data.morocco<- read.csv("crepon_morocco_analysis.csv")%>%
  select(treatment,client,expense_total )

obs <- nrow(data.morocco)

#Solo para ejemplificar, la siguiente línea pide el remuestreo, es decir, obtener una muestra de tamaño N de la muestra original, con reemplazo
data.b <-data.morocco[sample(nrow(data.morroco),obs, replace = TRUE),]


#2a
#Simplemente repetimos el proceso anterior B veces
set.seed(820)
B=50
#Inicializamos el vector donde guardaremos los W estimados
Wrep50_1 <- data.frame(W=matrix(ncol = 1, nrow = B))

for (i in 1:B)
{
  data.b <-data.morocco[sample(nrow(data.morroco),obs, replace = TRUE),]
  
  #Literalmente pegamos el cálculo de un W, hecho en la Pregunta 1a
  
  mean_cliente<-data.b %>%
    group_by(treatment) %>% 
    summarize(p_cliente=mean(client, na.rm=F)) %>% 
    ungroup()
  
  mean_gasto<-data.b %>%
    group_by(treatment) %>% 
    summarize(m_gasto=mean(expense_total, na.rm=F)) %>% 
    ungroup()
  
  dif_gasto <- mean_gasto[2,2]-mean_gasto[1,2]
  dif_cliente <- mean_cliente[2,2]-mean_cliente[1,2]
  
  #Y lo guardamos en la posición adecuada
  
  Wrep50_1[i,1] <- as.numeric(dif_gasto / dif_cliente)
}

#El error estimado es simplemente la desviación estándar de los B estadísticos estimados
sd(Wrep50_1$W)




#2b Simplemente cambiamos la semilla
set.seed(920)
B=50
#Inicializamos el vector donde guardaremos los W estimados
Wrep50_2 <- data.frame(W=matrix(ncol = 1, nrow = B))

for (i in 1:B)
{
  data.b <-data.morocco[sample(nrow(data.morroco),obs, replace = TRUE),]
  
  mean_cliente<-data.b %>%
    group_by(treatment) %>% 
    summarize(p_cliente=mean(client, na.rm=F)) %>% 
    ungroup()
  
  mean_gasto<-data.b %>%
    group_by(treatment) %>% 
    summarize(m_gasto=mean(expense_total, na.rm=F)) %>% 
    ungroup()
  
  dif_gasto <- mean_gasto[2,2]-mean_gasto[1,2]
  dif_cliente <- mean_cliente[2,2]-mean_cliente[1,2]
  
  Wrep50_2[i,1] <- as.numeric(dif_gasto / dif_cliente)
}

sd(Wrep50_2$W)



#2c regresamos a la primera semilla pero ahora B=1000
set.seed(820)
B=1000

Wrep1000 <- data.frame(W=matrix(ncol = 1, nrow = B))

for (i in 1:B)
{
  data.b <-data.morocco[sample(nrow(data.morroco),obs, replace = TRUE),]
  
  mean_cliente<-data.b %>%
    group_by(treatment) %>% 
    summarize(p_cliente=mean(client, na.rm=F)) %>% 
    ungroup()
  
  mean_gasto<-data.b %>%
    group_by(treatment) %>% 
    summarize(m_gasto=mean(expense_total, na.rm=F)) %>% 
    ungroup()
  
  dif_gasto <- mean_gasto[2,2]-mean_gasto[1,2]
  dif_cliente <- mean_cliente[2,2]-mean_cliente[1,2]
  
  Wrep1000[i,1] <- as.numeric(dif_gasto / dif_cliente)
}

#El error estimado es simplemente la desviación estándar de los B estadísticos estimados
sd(Wrep1000$W)



#Pregunta 4

#4a 
data.angrist<-read_csv("./STAR_public_use.csv",
                 locale = locale(encoding = "latin1"))   %>% 
  clean_names()

data.angrist<-data.angrist %>% 
  filter(noshow==0) %>% 
  mutate(credits_earned2=ifelse(is.na(prob_year1),NA,credits_earned2)) 
  
  
reg_credits_earned2 <-lm(credits_earned2 ~ ssp + sfp+ sfsp+
             factor(sex)+
             factor(mtongue)+
             factor(hsgroup)+
             factor(numcourses_nov1)+
             factor(lastmin)+
             factor(mom_edn)+
             factor(dad_edn),
           data.angrist)
coeftest(reg_credits_earned2, vcov = vcovHC(reg_credits_earned2, "HC1"))[1:4,]



#4b z-score para credits_earned2

credits_earned2_stats <- data.angrist %>% 
    filter(control==1) %>% 
    summarize(media=mean(credits_earned2,na.rm=T),desvest=sd(credits_earned2,na.rm=T))
  
data.angrist <- data.angrist %>% 
    mutate(credits_earned2_sd=(credits_earned2-credits_earned2_stats$media)/credits_earned2_stats$desvest)

#Media 0
data.angrist %>%
  filter(control==1) %>% 
  summarize(media=mean(credits_earned2_sd,na.rm=T))

#Desv. Std. 1
data.angrist %>%
  filter(control==1) %>% 
  summarize(desvest=sd(credits_earned2_sd,na.rm=T)) 




#4c Regresión con variable estandarizada
reg_credits_earned2_sd<-lm(credits_earned2_sd ~ ssp + sfp+ sfsp+
             factor(sex)+
             factor(mtongue)+
             factor(hsgroup)+
             factor(numcourses_nov1)+
             factor(lastmin)+
             factor(mom_edn)+
             factor(dad_edn),
           data.angrist)

coeftest(reg_credits_earned2_sd, vcov = vcovHC(reg_credits_earned2_sd, "HC1"))[1:4,]

#Noten que los estadísticos t son exactamente iguales
#La inferencia no cambia, solo la interpretación
coeftest(reg_credits_earned2, vcov = vcovHC(reg_credits_earned2, "HC1"))[1:4,]


#4d
gpa_year2_stats <- data.angrist %>% 
  filter(control==1) %>% 
  summarize(media=mean(gpa_year2,na.rm=T),desvest=sd(gpa_year2,na.rm=T))

data.angrist <- data.angrist %>% 
  mutate(gpa_year2_sd=(gpa_year2-gpa_year2_stats$media)/gpa_year2_stats$desvest)

#Agregamos las dos variables
data.angrist <- data.angrist %>% 
  mutate(promedio_vars=rowMeans(select(.,credits_earned2_sd,gpa_year2_sd)))

promedio_vars_stats <- data.angrist %>% 
  filter(control==1) %>% 
  summarize(media=mean(promedio_vars,na.rm=T),desvest=sd(promedio_vars,na.rm=T))

data.angrist <- data.angrist %>% 
  mutate(promedio_vars_sd=(promedio_vars-promedio_vars_stats$media)/promedio_vars_stats$desvest)


#Media 0
data.angrist %>%
  filter(control==1) %>% 
  summarize(media=mean(promedio_vars_sd,na.rm=T))

#Desv. Std. 1
data.angrist %>%
  filter(control==1) %>% 
  summarize(desvest=sd(promedio_vars_sd,na.rm=T)) 


#4e Regresión con el índice
reg_promedio_vars_sd<-lm(promedio_vars_sd ~ ssp + sfp+ sfsp+
                             factor(sex)+
                             factor(mtongue)+
                             factor(hsgroup)+
                             factor(numcourses_nov1)+
                             factor(lastmin)+
                             factor(mom_edn)+
                             factor(dad_edn),
                           data.angrist)

coeftest(reg_promedio_vars_sd, vcov = vcovHC(reg_credits_earned2_sd, "HC1"))[1:4,]






#Pregunta 5
data.pvalues<-read_csv("./pvalues.csv",
                       locale = locale(encoding = "latin1"))  

alpha <- 0.05


#5a Corrección con familias originales
#Ordenar de menor a mayor por familia
data.fam.or <- data.pvalues
setorder(data.fam.or,familia,p) #usé setorder de la libreria data.table

#Número de posición
data.fam.or <- data.fam.or %>% 
  group_by(familia) %>% 
  mutate(posicion=seq(along.with = familia)) %>% 
  mutate(numerohipotesis=max(posicion)) %>% 
  ungroup()

#Reglas (1 = rechazar H0, 0 = no rechazar)
data.fam.or <- data.fam.or %>% 
  mutate(regla_sincorr=ifelse(p<.05,1,0)) %>% 
  mutate(alpha_bonferroni=alpha/numerohipotesis) %>% 
  mutate(corrector_bh=posicion/numerohipotesis) %>% 
  mutate(q=corrector_bh*alpha) %>% 
  mutate(regla_bonferroni=ifelse(p<alpha_bonferroni,1,0)) %>% 
  mutate(regla_bh=ifelse(p<q,1,0))

setorder(data.fam.or,hipotesis)


data.fam.or %>% 
  select(hipotesis,regla_sincorr,regla_bonferroni,regla_bh)




#5b Ahora familias es la familia corregida
data.fam.corr <- data.pvalues
setorder(data.fam.corr,familia_corregida,p) #usé setorder de la libreria data.table

#Número de posición
data.fam.corr <- data.fam.corr %>% 
  group_by(familia_corregida) %>% 
  mutate(posicion=seq(along.with = familia_corregida)) %>% 
  mutate(numerohipotesis=max(posicion)) %>% 
  ungroup()

#Reglas (1 = rechazar H0, 0 = no rechazar)
data.fam.corr <- data.fam.corr %>% 
  mutate(regla_sincorr=ifelse(p<.05,1,0)) %>% 
  mutate(alpha_bonferroni=alpha/numerohipotesis) %>% 
  mutate(corrector_bh=posicion/numerohipotesis) %>% 
  mutate(q=corrector_bh*alpha) %>% 
  mutate(regla_bonferroni=ifelse(p<alpha_bonferroni,1,0)) %>% 
  mutate(regla_bh=ifelse(p<q,1,0))

setorder(data.fam.corr,hipotesis)

data.fam.corr %>% 
  select(hipotesis,regla_sincorr,regla_bonferroni,regla_bh)



#5c Sin considerar las familias
data.nofam <- data.pvalues
setorder(data.nofam,p)

#Es como si solo hubiera una sola familia
#Hago esto para usar mi mismo código
data.nofam <- data.nofam %>% 
  mutate(granfamilia=1)

#Número de posición
data.nofam <- data.nofam %>% 
  group_by(granfamilia) %>% 
  mutate(posicion=seq(along.with = granfamilia)) %>% 
  mutate(numerohipotesis=max(posicion)) %>% 
  ungroup()

#Reglas (1 = rechazar H0, 0 = no rechazar)
data.nofam <- data.nofam %>% 
  mutate(regla_sincorr=ifelse(p<.05,1,0)) %>% 
  mutate(alpha_bonferroni=alpha/numerohipotesis) %>% 
  mutate(corrector_bh=posicion/numerohipotesis) %>% 
  mutate(q=corrector_bh*alpha) %>% 
  mutate(regla_bonferroni=ifelse(p<alpha_bonferroni,1,0)) %>% 
  mutate(regla_bh=ifelse(p<q,1,0))

setorder(data.nofam,hipotesis)

data.nofam %>% 
  select(hipotesis,regla_sincorr,regla_bonferroni,regla_bh)

  