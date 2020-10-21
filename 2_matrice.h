#ifndef MATRICE_H
#define MATRICE_H

#include "1_vecteur.h"
#include "1_vecteur.cpp"
#include <cmath>
#include <iostream>

using uint = unsigned int;
class Matrice
{
public:
  void affiche();
  Matrice(uint l, uint c);
  Matrice(Vecteur d);
  Matrice(Vecteur *vecteurs, uint taille);
  ~Matrice();                    //destructeur
  Matrice(const Matrice &other); //constructeur de recopie
  Matrice operator=(const Matrice &other);
  Vecteur &operator[](uint column); //accès à une colonne de la matrice (un vecteur)
  Matrice operator+(const Matrice &other);
  Matrice operator-(const Matrice &other);
  Matrice operator*(const Matrice &other); //produit de deux matrices
  Vecteur mvprod(Vecteur vec);
  Matrice transpose();
  Matrice submat(uint il, uint jl, uint ic, uint jc);
  uint get_nb_lines();
  uint get_nb_columns();

private:
  Vecteur **mat = nullptr;
  uint dims[2];
};

float norm(Matrice mat);
Matrice operator*(float a, Matrice mat);
Matrice outer(Vecteur u, Vecteur v);
void setPrecision(Matrice &A); //arrondir à 0 avec une precision 10^-6

#endif
