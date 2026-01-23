####################################
# import des librairies necessaires
####################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

#####################
# paramètres globaux
#####################

# paramètres de l'option
S0 = 100  
K = 100   
T = 1.0   
r = 0.05  
sigma = 0.2  


# paramètres pour le schéma explicite 
M = 800
N = 180000  
z = norm.ppf(0.99)
l = np.abs((r - 0.5 * sigma**2)) * T + (z * sigma * np.sqrt(T)) 
delta_min = np.log(S0) - l
delta_max = np.log(S0) + l
delta = (delta_max - delta_min) / M 
h = T / N  
delta_values = np.linspace(delta_min, delta_max, M+1) 

###########################
# définition des fonctions
###########################

def stabilite_explicite_BS(sigma, h, delta):
    if h < (delta**2)/ (sigma**2):
        return True
    else:
        return False


def call_bs(S,K,T,r,sigma):
    d1 = (np.log(S/K) + (r+ 0.5* sigma**2)*T)/ (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    call_price = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    return call_price


def schema_explicite_BS_log(S0, K, T, r, sigma, M, N,h,delta,delta_values,delta_max):
    S_values = np.exp(delta_values)
    u_price = np.maximum(S_values - K, 0.0)
    A = (1 -r*h - (sigma**2 * h)/(delta**2))
    B = ((h*sigma**2)/(2*delta**2) + h*(r-0.5*sigma**2)/(2*delta))
    C = ((h*sigma**2)/(2*delta**2) - h*(r-0.5*sigma**2)/(2*delta))
    u_new = np.zeros(M+1)

    for i in range(N,0,-1):
        new_t = (i-1) * h
        # mise à jour des prix
        # [1:-1] de l'indice 1 à l'indice M-1, [2:] de l'indice 2 à M, [:-2] de l'indice 0 à M-2
        u_new[1:-1] = (u_price[1:-1]*A + u_price[2:]*B + u_price[:-2]*C)

        # conditions aux bornes (Dirichlet)
        u_new[0] = 0.0
        u_new[M] = np.exp(delta_max) - K * np.exp(-r *( T - new_t))

        # mise à jour des prix
        u_price[:] = u_new[:]
   
    price = u_price[M//2]  
    
    return price
       

def main():
    ################################################################################
    # Calcul des prix pour différentes valeurs de M et vérification de la stabilité
    ################################################################################

    M_values = np.array([50,100,800,1000,1200,1400,1600,1800,2000])
    price_explicite = np.zeros(len(M_values))
    erreur_explicite = np.zeros(len(M_values))
    stable = np.zeros(len(M_values), dtype=bool)
    bs_price = call_bs(S0,K,T,r,sigma)


    for i, m in enumerate(M_values) :
        delta_m = (delta_max - delta_min) / m
        delta_values_m = np.linspace (delta_min, delta_max, m+1)
        stable[i] = stabilite_explicite_BS(sigma, h, delta_m)

        if stable[i] == False :
            print(f"Stabilité non respectée pour M={m}")
            price_explicite[i] = np.nan
            erreur_explicite[i] = np.nan
        else: 
            print(f"Stabilité respectée pour M={m}")
            price_explicite[i] = schema_explicite_BS_log(S0, K, T, r, sigma, m, N,h,delta_m,delta_values_m,delta_max)
            erreur_explicite[i] = np.abs(price_explicite[i] - bs_price)


    ##########################
    # affichage des résultats 
    ##########################

    plt.figure(figsize=(12,6))
    plt.plot(M_values, price_explicite, 'o-', label='Prix schéma explicite')
    plt.axhline(y=bs_price, color='red', linestyle='--', label='Prix Black-Scholes fermé')
    plt.xlabel("points d'espace M") 
    plt.ylabel("Prix")
    plt.title("Convergence du schéma explicite vers le prix analytique")
    plt.legend()
    plt.grid()
    plt.gca().ticklabel_format(style='plain', axis='y', useOffset=False)
    plt.tight_layout()
    plt.show()


    plt.figure(figsize=(12,6))
    plt.plot(M_values, erreur_explicite, 'o-')
    plt.xlabel("points d'espace M")
    plt.ylabel("Erreur")
    plt.title("évolution de l'erreur du schéma explicite en fonction de M")
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True, which='both')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()