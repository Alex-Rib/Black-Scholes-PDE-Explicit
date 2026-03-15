# Schéma Explicite pour l'EDP de Black-Scholes

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Finance](https://img.shields.io/badge/Finance-Derivatives-green)
![Status](https://img.shields.io/badge/Status-Educational-orange)

## 📊 Description

Implémentation d'un **schéma aux différences finies (Schéma Explicite)** pour résoudre l'équation aux dérivées partielles (EDP) de Black-Scholes et pricer un Call Européen.


## 🎯 Objectifs

- Résoudre numériquement l'EDP de Black-Scholes par différences finies
- Étudier la convergence vers la solution analytique (formule de Black-Scholes)
- Analyser la condition de stabilité du schéma explicite
- Étudier l'ordre de convergence de l'erreur en fonction du pas d'espace $\delta$.

## 📐 Modèle Mathématique

### EDP de Black-Scholes transformée

En posant  $x = \ln(S)$, le prix de l'option $u(t,x)$ satisfait l'EDP suivante :

$$\frac{\partial u}{\partial t} + \left(r - \frac{1}{2}\sigma^2\right)\frac{\partial u}{\partial x} + \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2} - ru = 0$$

### Conditions aux bornes (Dirichlet) et Terminale

- **Condition Terminale** :
  $$u(T, x) = \max(e^x - K, 0)$$

- **Borne inférieur ($x \to x_{min}$)** : $u = 0$
- **Borne supérieur ($x \to x_{max}$)** : $u = e^x - K e^{-r(T-t)}$

## 🔧 Méthode Numérique

### Discrétisation

Nous notons les pas de discrétisation selon les notations du cours :
- **Pas de temps** : $h = T/N$
- **Pas d'espace** : $\delta = (x_{max} - x_{min})/M$

### Schéma Explicite 

La valeur $u_i^{n-1}$ (au temps précédent) est calculée explicitement à partir des valeurs au temps $n$ :

$$u_i^{n-1} = A \cdot u_i^n + B \cdot u_{i+1}^n + C \cdot u_{i-1}^n$$

Les coefficients $A, B, C$ sont donnés par :

$$A = 1 - rh - \frac{\sigma^2 h}{\delta^2}$$

$$B = \frac{h\sigma^2}{2\delta^2} + \frac{h(r - 0.5\sigma^2)}{2\delta}$$

$$C = \frac{h\sigma^2}{2\delta^2} - \frac{h(r - 0.5\sigma^2)}{2\delta}$$

### Condition de stabilité 

Pour garantir la stabilité du schéma explicite, la condition suivante doit être respectée :

$$h < \frac{\delta^2}{\sigma^2}$$

Le script vérifie automatiquement cette condition : si elle n'est pas respectée pour un $M$ donné, le calcul est ignoré.

## 📊 Paramètres numériques du Modèle


- **M** : Nombre de points d'espace (testé de 50 à 2000)
- **N** : Nombre de pas de temps (180 000)
- **Intervalle** : Domaine du log-prix à 99% de confiance

## 📈 Résultats

Le script génère les visualisations suivantes :

1. **Convergence du Prix** : Comparaison entre le prix obtenue par schéma et le prix analytique (Black-Scholes) pour différentes finesses de grille ($M$).
2. **Analyse d'Erreur** : Graphique log-log montrant la décroissance de l'erreur absolue en fonction de $M$.


## 📚 Prérequis

```bash
pip install numpy matplotlib scipy
```

## 🔜 Perspectives

Ce projet a pour vocation d'être étendu aux méthodes plus avancées :
- **Crank-Nicolson** : Pour améliorer la convergence temporelle (ordre 2 en temps contre ordre 1 pour l'explicite).


## 👨‍💻 Auteur

Alexandre R. - Université Paris Cité
