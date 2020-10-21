#ifndef TENSEUR_H
#define TENSEUR_H

#include "1_vecteur.h"
#include "2_matrice.h"
#include <cmath>
#include <iostream>

using uint = unsigned int;

class Tenseur
{
public:
  Tenseur(uint *tab, uint taille);
  Tenseur(Vecteur vecteur, uint *tab, uint taille);
  Tenseur(Matrice A, uint *tab, uint taille, uint k);
  ~Tenseur();                    // destructeur
  Tenseur(const Tenseur &other); // constructeur de recopie
  Tenseur operator=(const Tenseur &other);
  float &operator[](uint index);
  Tenseur operator+(const Tenseur &other);
  Tenseur operator-(const Tenseur &other);
  Matrice mode(uint k);
  void affiche();
  uint get_nbelts();
  uint get_ordre();
  uint *get_tab_dims();

private:
  uint ordre;     // ordre d
  uint *tab_dims; // [n_1,...,n_d]
  uint nbelts;
  Vecteur *tenseur_vec;
};

int modulo(int a, int b);
Tenseur pmod(Tenseur tenseur, Matrice matrice, uint k);
uint phi(uint ordre, uint *dims, uint *indices);            //fonction phi
void invphi(uint ordre, uint *dims, uint i, uint *indices); //fonction reciproque de phi

#endif
