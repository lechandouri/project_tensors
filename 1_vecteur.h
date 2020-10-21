#ifndef VECTEUR_H // inclure fichier une seule fois
#define VECTEUR_H

#include <cmath>
#include <iostream>

// rajouter friend?

using uint = unsigned int; // raccourci

class Vecteur
{
public:
  void affiche();
  Vecteur(uint dimension);
  Vecteur(float *coefficients, uint dimension);
  ~Vecteur();                              //destructeur
  Vecteur(const Vecteur &other);           //constructeur de recopie
  Vecteur operator=(const Vecteur &other); //operateur d'affectation
  Vecteur operator+(const Vecteur &other);
  Vecteur operator-(const Vecteur &other);
  float &operator[](uint line);
  Vecteur subvec(uint i, uint j); //extraction de Vecteur
  uint get_dim();

private:
  //membres de la classe en privée
  float *tab;
  uint dim;
};

float dot(Vecteur vec1, Vecteur vec2);
float norm(Vecteur vec);
Vecteur operator*(float a, Vecteur vec);

#endif
