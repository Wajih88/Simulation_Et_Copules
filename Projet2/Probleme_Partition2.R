# Initialisation des parametres -------------------------------------------

N = 2000 # nombre maximal que peut atteindre une valeur 
d = 1000 # nombre de valeurs  dans le vecteur




# Fonction Elementaires ---------------------------------------------------
fctH = function(vectValeurs,vectCoordones) {
  loc1 = sum(vectValeurs*vectCoordones)
  loc2 = sum(vectValeurs * (1-vectCoordones))
  return(abs(loc1 - loc2))
}


# Premiere Transition -----------------------------------------------------

transition1 = function(vectValeurs,vectCoord,betta) {
  # on tire une valeur au hasard et on change son appartenance 
  coord = sample(1:d,1)
  vectCoordNew = vectCoord
  vectCoordNew[coord] = 1-vectCoordNew[coord]
  Delta = fctH(vectValeurs,vectCoordNew)-fctH(vectValeurs,vectCoord)
  Qxy = 1/length(which(vectCoord == vectCoord[coord]))
  Qyx = 1/(length(which(vectCoord == 1-vectCoord[coord]))+1)
  alpha = min(exp(-betta * Delta)*Qxy/Qyx,1)
  U = runif(1)
  if(U < alpha) {return(vectCoordNew)}
  else { return(vectCoord)}
}


# Metropolis premiere transition ------------------------------------------
# La transition semble bien marcher pour des vecteurs pris aleatoirement
vectValeurs = sample(1:N,d,replace = TRUE)
vectCoord = sample(0:1,d,replace = TRUE)
nmax1 = 100
betta0 = 0.05
vectCoordBest = vectCoord
Hbest = fctH(vectValeurs,vectCoord)
valeursH = c()
for (i in 1:nmax1) {
  betta = betta0 * sqrt(i)
  vectCoord = transition1(vectValeurs,vectCoord,betta)
  valeursH[i] = fctH(vectValeurs,vectCoord)
  if(valeursH[i]<Hbest){
    Hbest = valeursH[i]
    vectCoordBest = vectCoord
  }
}
plot(valeursH)
print(Hbest)



# Deuxieme Transition -----------------------------------------------------

transition2 = function(vectValeurs,vectCoord,betta) {
  bool = (sum(vectValeurs * vectCoord) > sum(vectValeurs* (1-vectCoord)))*1
  # bool = 1 si la somme associee aux 1 du vecteur coordonnes est plus grande
  # que la somme associee au 0 du vecteur coordonnes
  vectCoordNew = vectCoord
  coord = sample(which(vectCoord==bool),1)
  vectCoordNew[coord] = 1 - vectCoordNew[coord]
  Delta = fctH(vectValeurs,vectCoordNew)-fctH(vectValeurs,vectCoord)
  
  Qxy = 1/(length(which(vectCoord == vectCoord[coord])))
  # la somme des elements dans l'autre ensemble que coord + la valeur de celui-ci
  terme1 = sum(vectValeurs* (1-bool)) + vectValeurs[coord] 
  # la somme des elements dans le mm ensemble que coord - la valeur de celui-ci
  terme2 = sum(vectValeurs * bool) - vectValeurs[coord]
  if(terme1 < terme2) {
    # Dans ce cas la transition de y a x est impossible et la transition de x a y est toujours accepte
    Qyx = 0
  }
  else{
    Qyx = 1/(length(which(vectCoord == 1-vectCoord[coord]))+1)
  }
  alpha = min(exp(-betta * Delta)*Qxy/Qyx,1)  
  U = runif(1)
  if(U < alpha) {return(vectCoordNew)}
  else { return(vectCoord)}
}


# Metropolis deuxieme transition ------------------------------------------

vectValeurs = sample(1:N,d,replace = TRUE)
vectCoord = sample(0:1,d,replace = TRUE)
nmax2 = 1000
betta0 = 0.2
vectCoordBest = vectCoord
Hbest = fctH(vectValeurs,vectCoord)
valeursH = c()
for (i in 1:nmax2) {
  betta = betta0 * sqrt(i)
  vectCoord = transition2(vectValeurs,vectCoord,betta)
  valeursH[i] = fctH(vectValeurs,vectCoord)
  if(valeursH[i]<Hbest){
    Hbest = valeursH[i]
    vectCoordBest = vectCoord
  }
}
plot(valeursH)
print(Hbest)





# Test avec valeurs particulieres avec premiere transition ----------------
# La premiere transition n'arrive pas a trouver le minimum globale
vectValeurs = 1:N
vectValeurs[N+1] = N*(N+1)/2
vectCoord = sample(0:1,N+1,replace = TRUE)
nmax2 = 1000
betta0 = 0.01
vectCoordBest = vectCoord
Hbest = fctH(vectValeurs,vectCoord)
valeursH = c()
for (i in 1:nmax2) {
  betta = betta0 * sqrt(i)
  vectCoord = transition1(vectValeurs,vectCoord,betta)
  valeursH[i] = fctH(vectValeurs,vectCoord)
  if(valeursH[i]<Hbest){
    Hbest = valeursH[i]
    vectCoordBest = vectCoord
  }
}
plot(valeursH)
print(Hbest)




# Test avec valeurs particulieres avec deuxieme transition -----------------------------------------
# la deuxieme transition arrive a atteindre le minimum global de ce vecteur extreme
vectValeurs = 1:N
vectValeurs[N+1] = N*(N+1)/2
vectCoord = sample(0:1,N+1,replace = TRUE)
nmax2 = 100
betta0 = 0.2
vectCoordBest = vectCoord
Hbest = fctH(vectValeurs,vectCoord)
valeursH = c()
for (i in 1:nmax2) {
  betta = betta0 * sqrt(i)
  vectCoord = transition2(vectValeurs,vectCoord,betta)
  valeursH[i] = fctH(vectValeurs,vectCoord)
  if(valeursH[i]<Hbest){
    Hbest = valeursH[i]
    vectCoordBest = vectCoord
  }
}
plot(valeursH)
print(Hbest)