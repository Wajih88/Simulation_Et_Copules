
# Initialisation des parametres -------------------------------------------

N = 500 # nombre maximal que peut atteindre une valeur dans les deux vecteurs
d = 50 # nombre initial de valeurs  dans les deux vecteurs





# fonctions utiles --------------------------------------------------------

Hvaleur = function(vecteur1,vecteur2) { 
  loc1 = sum(vecteur1)
  loc2 = sum(vecteur2)
  return(abs(loc2-loc1))
}

transition1 = function(vecteur1,vecteur2,betta) {
  loc1 = sum(vecteur1)+sum(vecteur2)
  loc2 = Hvaleur(vecteur1,vecteur2)
  if ((sum(vecteur1) == sum(vecteur2)) | ( loc1 %% 2 == 1 & loc2 == 1 )) 
  {
    print('les deux vecteurs sont equilibrees')
    newList <- list("vecteur1" = vecteur1, "vecteur2" = vecteur2)
    return(newList)
  }
  
  
  else if(sum(vecteur2)>sum(vecteur1)) {
    # dans ce cas on va prendre un element de vecteur2 au hasard et on l'ajoute a vecteur1
    m = length(vecteur2)
    coord = sample(1:m,1)
    vectnew1 = vecteur1
    n = length(vecteur1)
    vectnew1[n+1] = vecteur2[coord]
    vectnew2 = vecteur2[-coord]

  }
  else if (sum(vecteur1)>sum(vecteur2)){
    # dans ce cas on va prendre un element de vecteur1 au hasard et on l'ajoute a vecteur2
    n = length(vecteur1)
    coord = sample(1:n,1)
    vectnew1 = vecteur1[-coord]
    m = length(vecteur2)
    vectnew2 = vecteur2
    vectnew2[m+1] = vecteur1[coord]
  }
  Delta = Hvaleur(vectnew1,vectnew2) - Hvaleur(vecteur1,vecteur2)
  alpha = min(exp(-betta * Delta),1)
  U = runif(1)
  if (U<alpha) {
    newList <- list("vecteur1" = vectnew1, "vecteur2" = vectnew2)
    return(newList)
  }
  else {
    newList <- list("vecteur1" = vecteur1, "vecteur2" = vecteur2)
    return(newList)
  }
}

# Algorithme de Metropolis ------------------------------------------------


# Initialisation des vecteurs 
vecteur1 = sample(1:N,d)
vecteur2 = sample(1:N,d)

nmax = 300
valeursH = c()
betta0 = 0.007
maliste = transition1(vecteur1,vecteur2,betta0)
vecteur1best = maliste$vecteur1
vecteur2best = maliste$vecteur2
Hbest = Hvaleur(maliste$vecteur1,maliste$vecteur2)
for (i in 1:nmax) {
  betta = betta0 * sqrt(i)
  maliste = transition1(maliste$vecteur1,maliste$vecteur2,betta)
  valeursH[i] = Hvaleur(maliste$vecteur1,maliste$vecteur2)
  if(valeursH[i]<Hbest){
    Hbest = valeursH[i]
    vecteur1best = maliste$vecteur1
    vecteur2best = maliste$vecteur2
  }
}
plot(valeursH)
print(Hbest)
