*This is a school project I did last year, on tensors and singular value decomposition.
I had to answer a series of questions, and the following text references those questions to help my professor better understand my mistakes.
It is in french.* 

##  Structure du projet

Chaque partie est divisée en trois fichiers, deux fichiers cpp (un main et un fichier contenant les differences fonctions) et un header.
Il faut effectuer le building sur (NUMERO_PARTIE)_main.cpp et compiler pour avoir les resultats des différentes questions. 

## Diverses problemes

### Partie 3 

[Dernière question] 

La matrice B ne donne pas la décomposition valide qui se trouve dans les slides de validation. Toutefois, les résultats sont bons pour les matrices A et B. Je ne sais pas vraiment où se trouve le problème sachant que le code marche pour les deux autres matrices.

### Partie 4

[Dernière question]

La fonction pmod() ne donne pas les résultats demandés. Il est possible que cela soit une erreur d'indices. Certains coefficients sont bons, mais pas au bon endroit. Il existe plusieurs 'cout' dans la fonction pour suivre l'évolution des indices i et j qui définissent les élements des tenseurs. Ces log ont été mis en commentaires dans le fichier '4_tenseur.cpp', mais peuvent être utilisés pour trouver l'erreur si elle est liée aux indices. 

### Partie 5

Par manque de temps, je n'ai pas pu écrire le code lié à cette partie. 




