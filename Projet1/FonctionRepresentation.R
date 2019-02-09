
# Packages necessaires ----------------------------------------------------


library(rgl) #package pour representer la fonction
library(GoFKernel) # pour utiliser la fonction inverse
#library(pracma)

# Donnees de la fonction a integrer ---------------------------------------


a = 2.5
b = 2.6
c = 4
d = 4.2
# la fonction est a integrer entre [a,b] et [c,d]
N = 5000
resultat = 1.881348



fonction <- function(X,Y) {
  r = (cos(X*Y-2*Y**2 )**2 )*(Y*X**4 + 2*(Y**2))
  return(r)
}


persp3d(fonction,xlim = c(a,b), ylim = c(c,d),col = "lightblue")



# Methode de rejet --------------------------------------------------------


K = d * b**4 + 2 * d**2 # est un majorant de f dans [a,b] * [c,d] on trouve K = 227.20992

Xnaif = runif(n =N,min = a,max = b )
Ynaif = runif(n = N,min = c,max = d)
Vnaif = runif(n = N,min = 0,max = K)
conNaif = (Vnaif < fonction(Xnaif,Ynaif))*1
monteCarloNaif = (b-a)*(d-c)* K * cumsum(conNaif)/ (1:N)
plot(monteCarloNaif, type = 'lines', col = 'blue',ylim = c(1.7,2))
abline(h = resultat)
# Fin de la methose de rejet






# Methode avec les probas uniformes ---------------------------------------


Xunif = runif(n =N,min = a,max = b )
Yunif = runif(n = N,min = c,max = d)
moneteCarloUnif = (b-a)* (d-c)* cumsum(fonction(X = Xunif,Y=Yunif)) / (1:N) 
lines(moneteCarloUnif,col = 'red')
# Fin de la methode avec les probas uniformes








# Methode avec premiere probabilite personnalise --------------------------


# D'apres le graphe de notre fonction, elle est presque constante selon l'axe des X et 
# parabolique en suivant l'axe des Y. Cela nous amene a considerer une densite de la forme
# f(x)g(y) ou f(x) est la probabilite uniforme sur [a,b] et 
# g des la forme : (y-4.1)^2 + epsilon. En effet, la probabilite ne doit pas s'annuler sur le pave considere
epsilon = 0.3
# on calcule la constante Const pour que const*f(x)*g(y) soit une probabilite sur [a,b] * [c,d]
Const = 1/((2/3)*0.1**3 + 0.2*epsilon)
pPer <- function(X,Y){
  return((1/(b-a)) * Const*((Y-4.1)**2+epsilon))
}

Xper = runif(n = N,min = a,max = b)
Uper = runif(n = N)
foRepartitionY = function(y) {
  return( (Const/3)*((y-4.1)**3+0.1**3)+(Const)*epsilon*(y-c)) 
}

#persp3d(pPer,xlim = c(a,b), ylim = c(c,d),col = "lightblue")

# on utilise la fonction inverse de R qui calcule numeriquement l'inverse d'une fonction donnee
# comme argument avec son domaine de definition.v En effet pour simuler les points Y on a besoin d'inverser la 
# fonction de repartition des Y
fInver = inverse(foRepartitionY,lower = c,upper = d)
Yper = c()
for (i in 1:N) {
  Yper[i] = fInver(Uper[i])
}

montePerso = fonction(Xper,Yper)/pPer(Xper,Yper)
monteCarloPerso = cumsum(montePerso) / (1:N)

lines(monteCarloPerso,col = 'green')

# Fin de la methode avec probabilite personnalise





# Deuxieme Methode de proba personnalise ----------------------------------



# Methode avec minimum en moyenne. Afin d'ameliorer la methode precedente, on a remarque que la fonction n'atteint pas 
# son minimum au point 4.1. le minimum bouge legerement en fonction de X. On a decide alors de faire la moyenne entre le
# point ou f atteint son minimum pour les valeurs extremes de x, a savoir a et b.
# de sorte que la densite s'ecrit maintenant : (constMoyen/(b-a))* ((y-ymoyen)^2 + epsilon)
# constMoyen est un scalaire multiplicatif pour normaliser et obtenir une densite sur [a,b]*[c,d]
# le reste des calculs est similaire a la mehode precedente
Y = seq(c,d,length = 100)

min1 = which.min(fonction(X = a,Y = Y))
ymin1 = Y[min1]

min2 = which.min(fonction(X = b,Y = Y))
ymin2 = Y[min2]

ymoyen = mean(c(ymin1,ymin2))

epsilonMoyen = 0.3
constMoyen = ( ( (d-ymoyen)**3 + (ymoyen - c)**3 )/3 + epsilonMoyen * (d-c) )**(-1)

pPerMoy <- function(X,Y){
  return(10 * constMoyen*((Y-ymoyen)**2+epsilonMoyen))
}

foRepartitionYmoyen = function(y) {
  return( (constMoyen/3)*((y-ymoyen)**3+(ymoyen-c)**3)+(constMoyen)*epsilonMoyen*(y-c)) 
}

fInverMoyen = inverse(foRepartitionYmoyen,lower = c,upper = d)

XperMoyen = runif(n = N,min = a,max = b)
UperMoyen = runif(n = N)

YperMoyen = c()
for (i in 1:N) {
  YperMoyen[i] = fInverMoyen(UperMoyen[i])
}

montePersoMoyen = fonction(XperMoyen,YperMoyen)/pPerMoy(XperMoyen,YperMoyen)
monteCarloPersoMoyen = cumsum(montePersoMoyen) / (1:N)
lines(monteCarloPersoMoyen, col = 'darkorchid')
# fin de la methode




# Methode avec copule de Franck avec des theta opposes et proba uniforme  --------


thetFranck = -1
cFranckOppose <- function(u,v) {
  r <- thetFranck * exp(-thetFranck * (u+v)) * (1 - exp(-thetFranck))
  r <- r / (exp(-thetFranck)-1 + (exp(-thetFranck*u)-1)*(exp(-thetFranck*v) - 1))**2
  s <- -thetFranck * exp(thetFranck * (u+v)) * (1 - exp(thetFranck))
  s <- s / (exp(thetFranck)-1 + (exp(thetFranck*u)-1)*(exp(thetFranck*v) - 1))**2
  return((r+s)/2)
}
#persp3d(cFranckOppose,xlim = c(0,1), ylim = c(0,1),col = "lightblue")

# simuler le couple de variables aleatoire (U,V)(notation cours page 17)
# pour cela on doit calculer la fonction de repartition conditionelle
# fonction de repartition conditionnelle
repartitionConditionelle <- function(v) {
  # thetFranck et u sont des parametres
  numerateur1 = (exp(-thetFranck * v)-1)*exp(-thetFranck * u)
  denominateur1 = exp(-thetFranck)-1 + (exp(-thetFranck * v)-1)*(exp(-thetFranck * u)-1)
  somme1 = numerateur1/denominateur1
  numerateur2 = (exp(thetFranck * v)-1)*exp(thetFranck * u)
  denominateur2 = exp(thetFranck)-1 + (exp(thetFranck * v)-1)*(exp(thetFranck * u)-1)
  somme2 = numerateur2/denominateur2
  return((somme1+somme2)/2)
}
U = runif(N)
Z = runif(N)
V <- c()
for (i in 1:N) {
  u = U[i] # mise a jour du point u pour chaque V[i]
  fInverFranck = inverse(repartitionConditionelle,lower = 0,upper = 1)
  V[i] = fInverFranck(Z[i])
}



#Simuler les variables aleatoires (X,Y) (notation page 17)
X  = qunif(U,min =a,max = b)
Y = qunif(V,min =c,max = d)


# Calculer la densite (formule page 11)

densiteFranckOppose <- function(x,y) {
  # (x,y) dans [2.5,2.6] x [4,4.2]
  fx = 1/(b-a)
  gy = 1/(d-c)
  Fx = (x-a)/(b-a)
  Gy = (y-c)/(d-c)
  return(fx*gy*cFranckOppose(Fx,Gy))
}
monteFranckOppose = fonction(X,Y)/densiteFranckOppose(X,Y)
monteFranckOpposeMoyen = cumsum(monteFranckOppose) / (1:N)
lines(monteFranckOpposeMoyen, col = 'aquamarine4')

# Fin de la methode de copule des Franck avec probas uniformes











# Methode copule de Franck avec des theta opposes et une marginale uniforme pour X et parabolique pour Y--------



thetFranck6 = -1
cFranckOppose6 <- function(u,v) {
  r <- thetFranck6 * exp(-thetFranck6 * (u+v)) * (1 - exp(-thetFranck6))
  r <- r / (exp(-thetFranck6)-1 + (exp(-thetFranck6*u)-1)*(exp(-thetFranck6*v) - 1))**2
  s <- -thetFranck6 * exp(thetFranck6 * (u+v)) * (1 - exp(thetFranck6))
  s <- s / (exp(thetFranck6)-1 + (exp(thetFranck6*u)-1)*(exp(thetFranck6*v) - 1))**2
  return((r+s)/2)
}


# simuler le couple de variables aleatoire (U,V)(notation cours page 17)

# fonction de repartition conditionnelle
repartitionConditionelle6 <- function(v) {
  # thetFranck6 et u sont des parametres
  numerateur1 = (exp(-thetFranck6 * v)-1)*exp(-thetFranck6 * u)
  denominateur1 = exp(-thetFranck6)-1 + (exp(-thetFranck6 * v)-1)*(exp(-thetFranck6 * u)-1)
  somme1 = numerateur1/denominateur1
  numerateur2 = (exp(thetFranck6 * v)-1)*exp(thetFranck6 * u)
  denominateur2 = exp(thetFranck6)-1 + (exp(thetFranck6 * v)-1)*(exp(thetFranck6 * u)-1)
  somme2 = numerateur2/denominateur2
  return((somme1+somme2)/2)
}
U = runif(N)
Z = runif(N)
V <- c()
for (i in 1:N) {
  u = U[i]
  fInverFranck = inverse(repartitionConditionelle6,lower = 0,upper = 1)
  V[i] = fInverFranck(Z[i])
}

#Simuler les variables aleatoires (X,Y) (notation page 17)
X  = qunif(U,min =a,max = b)
for (i in 1:N) {
  Y[i] = fInverMoyen(V[i]) # fonction d'inversion de la quatrieme methode
}

# Calculer la densite (formule page 11)

densiteFranckOppose6 <- function(x,y) {
  # (x,y) dans [2.5,2.6] x [4,4.2]
  fx = 1/(b-a)
  gy = constMoyen*((y-ymoyen)**2+epsilonMoyen)
  Fx = (x-a)/(b-a)
  Gy = foRepartitionYmoyen(y)
  return(fx*gy*cFranckOppose(Fx,Gy))
}

#persp3d(densiteFranckOppose6,xlim = c(2.5,2.6), ylim = c(4,4.2),col = "lightblue")


monteFranckOppose6 = fonction(X,Y)/densiteFranckOppose6(X,Y)
monteFranckOpposeMoyen = cumsum(monteFranckOppose) / (1:N)
lines(monteFranckOpposeMoyen, col = 'chocolate2')





legend("bottomright", legend=c("rejet", "proba uniforme","proba personalise",
                               "proba personalise2","copule marginale uniforme",
                               "copule marginale uniforme et parabolique"),
       col=c("blue", "red",'green','darkorchid','aquamarine4','chocolate2'),
       lty=c(1,2,3,4,5,6))




# Tracer les histogrammes et calculer les ecarts types --------------------




histo = 200
rejet <- c()
proba_unif <- c()
methode_perso <- c()
methode_perso2 <- c()
copule_uniforme <- c()
copule_parabol <- c()
for (i in 1:histo) {
  # premiere methode
  Xnaif = runif(n =N,min = a,max = b )
  Ynaif = runif(n = N,min = c,max = d)
  Vnaif = runif(n = N,min = 0,max = K)
  conNaif = (Vnaif < fonction(Xnaif,Ynaif))*1
  rejet[i] = (b-a)*(d-c)* K * sum(conNaif)/ N
  # deuxieme methode
  Xunif = runif(n =N,min = a,max = b )
  Yunif = runif(n = N,min = c,max = d)
  moneteCarloUnif = (b-a)* (d-c)* sum(fonction(X = Xunif,Y=Yunif)) / N
  proba_unif[i] = moneteCarloUnif
  
  # troisieme methode
  Xper = runif(n = N,min = a,max = b)
  Uper = runif(n = N)
  Yper = c()
  for (j in 1:N) {
    Yper[j] = fInver(Uper[j])
  }
  montePerso = fonction(Xper,Yper)/pPer(Xper,Yper)
  methode_perso[i]  = sum(montePerso) / N
  
  # quatrieme methode
  XperMoyen = runif(n = N,min = a,max = b)
  UperMoyen = runif(n = N)
  YperMoyen = c()
  for (j in 1:N) {
    YperMoyen[j] = fInverMoyen(UperMoyen[j])
  }
  montePersoMoyen = fonction(XperMoyen,YperMoyen)/pPerMoy(XperMoyen,YperMoyen)
  methode_perso2[i] = sum(montePersoMoyen) / N
  
  
  # cinquieme methode
  U = runif(N)
  Z = runif(N)
  V <- c()
  for (j in 1:N) {
    u = U[j]
    fInverFranck = inverse(repartitionConditionelle,lower = 0,upper = 1)
    V[j] = fInverFranck(Z[j])
  }
  X  = qunif(U,min =a,max = b)
  Y = qunif(V,min =c,max = d)
  monteFranckOppose = fonction(X,Y)/densiteFranckOppose(X,Y)
  copule_uniforme[i] = sum(monteFranckOppose) / (N)
  
  # sixieme methode
  U = runif(N)
  Z = runif(N)
  V <- c()
  for (j in 1:N) {
    u = U[j]
    fInverFranck = inverse(repartitionConditionelle,lower = 0,upper = 1)
    V[j] = fInverFranck(Z[j])
  }
  X  = qunif(U,min =a,max = b)
  for (j in 1:N) {
    Y[j] = fInverMoyen(V[j]) # fonction d'inversion de la quatrieme methode
  }
  monteFranckOppose6 = fonction(X,Y)/densiteFranckOppose6(X,Y)
  copule_parabol[i] = sum(monteFranckOppose6) / (N)
  
}

par(mfrow=c(1,6))
hist(rejet)
hist(proba_unif)
hist(methode_perso)
hist(methode_perso2)
hist(copule_uniforme)
hist(copule_parabol)


sd(rejet)
sd(proba_unif)
sd(methode_perso)
sd(methode_perso2)
sd(copule_uniforme)
sd(copule_parabol)
# Commentaire : la methode de rejet est la moins performante. Les autres methodes ont des performances assez comparables
# leurs ecarts types se situent aux environs de 0.02.
