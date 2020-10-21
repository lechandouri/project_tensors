#ifndef VECTEUR_CPP // inclure fichier une seule fois
#define VECTEUR_CPP

#include "1_vecteur.h"
using namespace std;

void Vecteur::affiche()
{
  cout << "dim: " << dim;
  cout << " - vec: [ ";

  for (uint i = 0; i < dim; i++)
  {
    cout << tab[i] << ' ';
  }

  cout << ']' << endl;
}

Vecteur::Vecteur(uint dimension) : dim(dimension)
{
  tab = new float[dimension];
  for (uint i = 0; i < dimension; i++) // pour forcer le zéro exact
  {
    tab[i] = 0.;
  }
}

Vecteur::Vecteur(float *coefficients, uint dimension) : dim(dimension)
{
  tab = new float[dimension];

  for (uint i = 0; i < dimension; i++)
  {
    tab[i] = coefficients[i]; // initialisation de chaque element du tableau
  }
}

Vecteur::~Vecteur() { delete[] tab; }

Vecteur::Vecteur(const Vecteur &other) : dim(other.dim)
{
  tab = new float[other.dim];

  for (uint i = 0; i < other.dim; i++)
  {
    tab[i] = other.tab[i];
  }
}

Vecteur Vecteur::operator=(const Vecteur &other)
{
  delete[] tab;
  dim = other.dim;
  tab = new float[dim];

  for (uint i = 0; i < dim; i++)
  {
    tab[i] = other.tab[i];
  }

  return *this;
}

Vecteur Vecteur::operator+(const Vecteur &other)
{
  if (other.dim != dim)
    throw runtime_error("Les vecteurs n'ont pas la meme taille");

  float coefficients[dim];

  for (uint i = 0; i < dim; i++)
  {
    coefficients[i] = tab[i] + other.tab[i];
  }

  return Vecteur(coefficients, dim);
}

Vecteur Vecteur::operator-(const Vecteur &other)
{
  if (other.dim != dim)
    throw runtime_error("Les vecteurs n'ont pas la meme taille");

  float coefficients[dim];

  for (uint i = 0; i < dim; i++)
  {
    coefficients[i] = tab[i] - other.tab[i];
  }

  return Vecteur(coefficients, dim);
}

float &Vecteur::operator[](uint line)
{
  if (line >= dim)
    throw runtime_error("L'indice excede la taille du vecteur");

  return tab[line];
}

Vecteur Vecteur::subvec(uint i, uint j)
{
  if (i > j)
    throw runtime_error(
        "L'indice initial doit etre plus petit ou egal a l'indice final");

  if (j >= dim)
    throw runtime_error("L'indice final excede la taille du vecteur");

  float coefficients[j - i + 1];

  for (uint k = i; k <= j; k++)
  {
    coefficients[k - i] = tab[k];
  }

  return Vecteur(coefficients, j - i + 1);
}

uint Vecteur::get_dim() { return dim; }

float dot(Vecteur vec1, Vecteur vec2)
{
  if (vec1.get_dim() != vec2.get_dim())
    throw runtime_error("Les vecteurs n'ont pas la meme taille");

  float s = 0.;

  for (uint i = 0; i < vec1.get_dim(); i++)
  {
    s = s + vec1[i] * vec2[i];
  }

  return s;
}

float norm(Vecteur vec)
{
  float result = 0.;

  for (uint i = 0; i < vec.get_dim(); i++)
  {
    result = result + pow(vec[i], 2);
  }

  return sqrt(result);
}

Vecteur operator*(float a, Vecteur vec)
{
  float coefficients[vec.get_dim()];

  for (uint i = 0; i < vec.get_dim(); i++)
  {
    coefficients[i] = a * vec[i];
  }

  return Vecteur(coefficients, vec.get_dim());
}

#endif
