#ifndef MSVD_H
#define MSVD_H

#include "1_vecteur.h"
#include "2_matrice.h"
#include <cmath>
#include <iostream>

using uint = unsigned int;

void cppidcolgivens(float x, float z, float &c, float &s);
Vecteur householder(Vecteur x, float &b);
Matrice reductridiag(Matrice &D);
void qrsym(Matrice &A, Matrice &Q);
Matrice qrpivot(Matrice &A, Matrice &Q);
void svd(Matrice A, Matrice decompositions[3]);
bool isEqualZero(float f);                //verifie si la variable arrondi est quasi-nulle
bool isDiagonal(Matrice &A);              //verifie si la matrice est diagonale
float get_maximum(float *array, int dim); //retourne l'element max d'un array

#endif
