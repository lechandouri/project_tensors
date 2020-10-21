#include "2_matrice.h"

using namespace std;

void Matrice::affiche()
{
  cout << "dims: (" << dims[0] << ", " << dims[1] << ") - mat:" << endl;

  for (uint l = 0; l < dims[0]; l++)
  {
    for (uint c = 0; c < dims[1]; c++)
    {
      cout << (*mat[c])[l] << ' ';
    }

    cout << endl;
  }
}

Matrice::Matrice(uint l, uint c) : dims{l, c}
{
  mat = new Vecteur
      *[c]; // 'c' correspond au nombres d'elements (nombre de colonnes)

  for (uint j = 0; j < c; j++)
  {
    mat[j] = new Vecteur(l);
  }
}

Matrice::Matrice(Vecteur d) //matrice diagonale
{
  uint dim = d.get_dim();
  dims[1] = dims[0] = dim;
  mat = new Vecteur *[dim];

  for (uint k = 0; k < dim; k++)
  {
    mat[k] = new Vecteur(dim);
    (*mat[k])[k] = d[k];
  }
}

Matrice::Matrice(Vecteur *vecteurs, uint dimension)
    : dims{vecteurs[0].get_dim(), dimension}
{
  mat = new Vecteur *[dimension];

  for (uint j = 0; j < dimension; j++)
  {
    mat[j] = new Vecteur(vecteurs[j]);
  }
}

Matrice::~Matrice() { delete[] mat; }

Matrice::Matrice(const Matrice &other) : dims{other.dims[0], other.dims[1]}
{
  mat = new Vecteur *[dims[1]];

  for (uint j = 0; j < dims[1]; j++)
  {
    mat[j] = new Vecteur(*other.mat[j]);
  }
}

Matrice Matrice::operator=(const Matrice &other)
{
  delete[] mat;
  dims[0] = other.dims[0];
  dims[1] = other.dims[1];
  mat = new Vecteur *[dims[1]];

  for (uint j = 0; j < dims[1]; j++)
  {
    mat[j] = new Vecteur(*other.mat[j]);
  }

  return *this;
}

Vecteur &Matrice::operator[](uint column)
{
  if (column >= dims[1])
    throw runtime_error("L'indice excede le nombre de colonnes de la matrice");

  return *mat[column];
}

Matrice Matrice::operator+(const Matrice &other)
{
  if (dims[0] != other.dims[0] || dims[1] != other.dims[1])
    throw runtime_error("Les matrices n'ont pas la meme dimension");

  Matrice result(dims[0], dims[1]);

  for (uint j = 0; j < dims[1]; j++)
  {
    result.mat[j] = new Vecteur(*mat[j] + *other.mat[j]);
  }

  return result;
}

Matrice Matrice::operator-(const Matrice &other)
{
  if (dims[0] != other.dims[0] || dims[1] != other.dims[1])
    throw runtime_error("Les matrices n'ont pas la meme dimension");

  Matrice result(dims[0], dims[1]);

  for (uint j = 0; j < dims[1]; j++)
  {
    result.mat[j] = new Vecteur(*mat[j] - *other.mat[j]);
  }

  return result;
}

Matrice Matrice::operator*(const Matrice &other)
{
  if (dims[1] != other.dims[0])
    throw runtime_error(
        "Les matrices ne peuvent pas etre multipliees entre elles");

  Matrice result(dims[0], other.dims[1]);
  // dims[0] : nombre de lignes, dims[1] : nombre de colonnes
  float coeff_line[dims[1]]; // elements de la ligne i

  for (uint i = 0; i < dims[0]; i++)
  {
    for (uint j = 0; j < dims[1]; j++)
    {
      coeff_line[j] = (*mat[j])[i];
    }

    // vecteur de taille dims[1] correspondant à la ligne i de la matrice d
    //  bas
    Vecteur line(coeff_line, dims[1]);
    for (uint j = 0; j < other.dims[1]; j++)
    {
      result[j][i] =
          dot(line, *other.mat[j]); // element (i,j) du produit des matrices
    }
  }

  return result;
}

Vecteur Matrice::mvprod(Vecteur vec)
{
  if (dims[1] != vec.get_dim())
    throw runtime_error(
        "La matrice et le vecteur ne peuvent pas etre multiplies entre eux");

  float result[dims[1]];
  float coeff_line[dims[1]];

  for (uint i = 0; i < dims[0]; i++)
  {
    for (uint j = 0; j < dims[1]; j++)
    {
      coeff_line[j] = (*mat[j])[i];
    }

    result[i] = dot(Vecteur(coeff_line, dims[1]), vec);
  }

  return Vecteur(result, dims[1]);
}

Matrice Matrice::transpose()
{
  Matrice result(dims[1], dims[0]);

  for (uint i = 0; i < dims[1]; i++)
  {
    for (uint j = 0; j < dims[0]; j++)
    {
      result[j][i] = (*mat[i])[j];
    }
  }

  return result;
}

Matrice Matrice::submat(uint il, uint jl, uint ic, uint jc)
{
  if (il > jl || ic > jc)
    throw runtime_error("Les indices initiaux doivent etre plus petits ou "
                        "egaux aux indices finaux");

  if (jl >= dims[0] || jc >= dims[1])
    throw runtime_error("Les dimensions de la sous-matrice sont trop grandes");

  Matrice result(jl - il + 1, jc - ic + 1);

  for (uint l = il; l <= jl; l++)
  {
    for (uint c = ic; c <= jc; c++)
    {
      result[c - ic][l - il] = (*mat[c])[l];
    }
  }

  return result;
}

uint Matrice::get_nb_lines() { return dims[0]; }

uint Matrice::get_nb_columns() { return dims[1]; }

float norm(Matrice mat)
{
  uint i_max = mat.get_nb_lines();
  uint j_max = mat.get_nb_columns();
  float result = 0;

  for (uint i = 0; i < i_max; i++)
  {
    for (uint j = 0; j < j_max; j++)
    {
      result = result + abs(mat[j][i]);
    }
  }

  return sqrt(result);
}

Matrice operator*(float a, Matrice mat)
{
  uint i_max = mat.get_nb_lines();
  uint j_max = mat.get_nb_columns();
  Matrice result(i_max, j_max);

  for (uint i = 0; i < i_max; i++)
  {
    for (uint j = 0; j < j_max; j++)
    {
      result[j][i] = a * mat[j][i];
    }
  }

  return result;
}

Matrice outer(Vecteur u, Vecteur v)
{
  if (u.get_dim() != v.get_dim())
    throw runtime_error("Les vecteurs n'ont pas la meme taille");

  uint dim = u.get_dim();
  Matrice result(dim, dim);

  for (uint i = 0; i < dim; i++)
  {
    for (uint j = 0; j < dim; j++)
    {
      result[j][i] = u[i] * v[j];
    }
  }

  return result;
}

void setPrecision(Matrice &A) //pour arrondir à 0
{
  for (int t = 0; t < A.get_nb_lines(); t++)
  {
    for (int s = 0; s < A.get_nb_columns(); s++)
    {
      if (abs(A[s][t]) <= pow(10, -6))
      {
        A[s][t] = 0;
      }
    }
  }
}
