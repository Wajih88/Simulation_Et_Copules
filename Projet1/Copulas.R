# Plotting some copulas
library(rgl)



# Copule de Ali-Mikhail-Haq (-1<= theta <1 )
# pour thetAli > 0 : concentration autour des valeurs (0,0) et (1,1)
# pour thetAlo < 0 : concentration autour des valeurs (0,1) et (1,0)

thetAli = -0.5
AliMikhail <- function(u,v) {
  nume = - ((-u*v + u + v -1)*thetAli**2 - thetAli * (u*v + u + v -2) - 1)
  return(nume / (( 1 - thetAli * (1-u) * (1-v) )**3) )
}
persp3d(AliMikhail,xlim = c(0,1), ylim = c(0,1),col = "lightblue")


# Copule de clayton (theta > 0) 
## Concentration autour des petites valeurs

thetClayton = 10
ClaytonDensity <- function(u,v) {
  return(  (thetClayton + 1) * ( (u*v)**(-thetClayton-1) ) * (u**(-thetClayton) + v**(-thetClayton) - 1 )**(-(2*thetClayton+1)/thetClayton)    )
}
persp3d(ClaytonDensity,xlim = c(0,1), ylim = c(0,1),col = "lightblue")


# Copule de Gumbel (theta >= 1)
# Concentration autour des grandes valeurs
thetaGumbel = 13
GumbelDensity <- function (u,v) {
  Ulo = (-log(u))**thetaGumbel
  Vlo = (-log(v))**thetaGumbel
  Ulo1 = (-log(u))**(thetaGumbel-1)
  Vlo1 = (-log(v))**(thetaGumbel-1)
  Copule = exp(-(Ulo + Vlo)**(1/thetaGumbel))
  loca = (Ulo + Vlo)**(1/thetaGumbel)
  loca2 = (Ulo + Vlo)**(1/thetaGumbel - 2)
  return( Ulo1 * Vlo1 * Copule *(thetaGumbel + loca - 1)  *loca2 /(u*v))
}
persp3d(GumbelDensity,xlim = c(0,1), ylim = c(0,1), col = "lightblue")

# Copule de Franck (theta != 0)
## pour theta > 0 concentration autour de (0,0) et (1,1)
## pour theta <0 concentration autour de (0,1) et (1,0)
thetFranck = -2

cFranck <- function(u,v) {
  r <- thetFranck * exp(-thetFranck * (u+v)) * (1 - exp(-thetFranck))
  r <- r / (exp(-thetFranck)-1 + (exp(-thetFranck*u)-1)*(exp(-thetFranck*v) - 1))**2
  return(r)
}

persp3d(cFranck,xlim = c(0,1), ylim = c(0,1),col = "lightblue")


# Copule de Joe (theta >= 1)
# concentration autour des grandes valeurs
thetaJoe = 3

JoeDensity <- function(u,v) {
  facteur1 = (1-u)**thetaJoe + (1-v)**thetaJoe - ((1-u)**thetaJoe)*((1-v)**thetaJoe)
  locu = (1-u)**(thetaJoe-1)
  locv = (1-v)**(thetaJoe-1)
  somme1 = thetaJoe * locu * locv * (facteur1**(1/thetaJoe-1))
  facteur2 = locu * ((1-(1-v))**thetaJoe)
  facteur3 = (1/thetaJoe-1)*locv*(-thetaJoe + thetaJoe*((1-u)**thetaJoe))*(facteur1**(1/thetaJoe-2))
  somme2 = facteur2*facteur3
  return(somme1 + somme2)
}

persp3d(JoeDensity,xlim = c(0,1), ylim = c(0,1),col = "lightblue")


# Methode avec moyenne de deux copules de Franck avec thetas opposes
# il y a concentration autour des points (0,0) (1,0) (0,1) (1,1)
thetFranck = -5
cFranckOppose <- function(u,v) {
  r <- thetFranck * exp(-thetFranck * (u+v)) * (1 - exp(-thetFranck))
  r <- r / (exp(-thetFranck)-1 + (exp(-thetFranck*u)-1)*(exp(-thetFranck*v) - 1))**2
  s <- -thetFranck * exp(thetFranck * (u+v)) * (1 - exp(thetFranck))
  s <- s / (exp(thetFranck)-1 + (exp(thetFranck*u)-1)*(exp(thetFranck*v) - 1))**2
  return((r+s)/2)
}
persp3d(cFranckOppose,xlim = c(0,1), ylim = c(0,1),col = "lightblue")
