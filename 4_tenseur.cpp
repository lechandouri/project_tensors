#include "4_tenseur.h"
#include "1_vecteur.cpp"
#include "1_vecteur.h"
#include "2_matrice.cpp"
#include "2_matrice.h"

using namespace std;

Tenseur::Tenseur(uint *tab, uint taille) : ordre(taille)
{
  tab_dims = new uint[taille];
  nbelts = 1;
  for (uint i = 0; i < taille; i++)
  {
    tab_dims[i] = tab[i];
    nbelts = nbelts * tab[i];
  }
  tenseur_vec = new Vecteur(nbelts);
}

Tenseur::Tenseur(Vecteur vecteur, uint *tab, uint taille) : ordre(taille)
{
  tab_dims = new uint[taille];
  nbelts = 1;
  for (uint i = 0; i < taille; i++)
  {
    tab_dims[i] = tab[i];
    nbelts = nbelts * tab[i];
  }
  tenseur_vec = new Vecteur(vecteur);
}

Tenseur::Tenseur(Matrice A, uint *tab, uint taille, uint k) : ordre(taille)
{
  tab_dims = new uint[taille];
  nbelts = 1;

  for (uint j = 0; j < taille; j++)
  {
    tab_dims[j] = tab[j];
    nbelts = nbelts * tab[j];
  }

  tenseur_vec = new Vecteur(nbelts);
  uint indices[ordre];
  uint column_k;
  uint line_k;
  uint indices_j[ordre - 1];
  uint tab_dims_j[ordre - 1];

  for (int t = 0; t < nbelts; t++)
  {
    invphi(ordre, tab_dims, t, indices);
    line_k = indices[k];

    // copie indices et tab_dims et supprimer l'element indice par k
    for (int s = 0; s <= ordre - 1; s++)
    {
      if (s < k)
      {
        indices_j[s] = indices[s];
        tab_dims_j[s] = tab_dims[s];
      }
      else
      {
        indices_j[s] = indices[s + 1];
        tab_dims_j[s] = tab_dims[s + 1];
      }
    }

    column_k = phi(ordre - 1, tab_dims_j, indices_j); // position colonne j_k
    (*tenseur_vec)[t] = A[column_k][line_k];
  }
}

Tenseur::~Tenseur() { delete[] tab_dims; }

Tenseur::Tenseur(const Tenseur &other) : ordre(other.ordre)
{
  tab_dims = new uint[other.ordre];
  for (uint i = 0; i < other.ordre; i++)
  {
    tab_dims[i] = other.tab_dims[i];
  }
  nbelts = other.nbelts;
  tenseur_vec = other.tenseur_vec;
}

Tenseur Tenseur::operator=(const Tenseur &other)
{
  delete[] tab_dims;
  ordre = other.ordre;
  tab_dims = new uint[other.ordre];
  for (uint i = 0; i < other.ordre; i++)
  {
    tab_dims[i] = other.tab_dims[i];
  }
  nbelts = other.nbelts;
  tenseur_vec = other.tenseur_vec;
  return *this;
}

float &Tenseur::operator[](uint index)
{
  if (index >= nbelts)
    throw runtime_error("L'indice excede la taille du vecteur");

  return (*tenseur_vec)[index];
}

Tenseur Tenseur::operator+(const Tenseur &other)
{
  for (uint i = 0; i < other.ordre; i++)
  {
    if (other.tab_dims[i] != tab_dims[i])
    {
      throw runtime_error("Les tenseurs n'ont pas les memes dimensions");
    }
  }
  return Tenseur((*tenseur_vec + *other.tenseur_vec), other.tab_dims,
                 other.ordre);
}

Tenseur Tenseur::operator-(const Tenseur &other)
{
  for (uint i = 0; i < other.ordre; i++)
  {
    if (other.tab_dims[i] != tab_dims[i])
    {
      throw runtime_error("Les tenseurs n'ont pas les memes dimensions");
    }
  }

  return Tenseur((*tenseur_vec - *other.tenseur_vec), other.tab_dims,
                 other.ordre);
}

void Tenseur::affiche()
{
  cout << "taille: " << nbelts;
  cout << " - tenseur vec: [ ";

  for (uint i = 0; i < nbelts; i++)
  {
    cout << (*tenseur_vec)[i] << ' ';
  }

  cout << ']' << endl;
}

uint Tenseur::get_nbelts() { return nbelts; }
uint Tenseur::get_ordre() { return ordre; }
uint *Tenseur::get_tab_dims() { return tab_dims; }

uint phi(uint ordre, uint *dims, uint *indices)
{
  int phi = indices[ordre - 1];
  uint n_ind = 1;

  for (int i = ordre - 2; i >= 0; i--)
  {
    n_ind = n_ind * dims[i + 1];
    phi = phi + n_ind * indices[i];
  }

  return phi;
}

int modulo(int a, int b) { return a % b; }

void invphi(uint ordre, uint *dims, uint i, uint *indices)
{
  int f[ordre];
  f[ordre - 1] = i;
  indices[ordre - 1] = modulo(f[ordre - 1], dims[ordre - 1]);

  for (int t = ordre - 2; t >= 0; t--)
  {
    f[t] = ((f[t + 1] - indices[t + 1] - 1) / dims[t + 1]);
    indices[t] = modulo(f[t], dims[t]);
  }
}

Matrice Tenseur::mode(uint k)
{
  Matrice T_mode(tab_dims[k], nbelts / tab_dims[k]);
  uint indices[ordre]; // (i_1,...,i_d) pour i fixe
  float valeur_i;
  uint line_k;
  uint column_k;
  uint indices_j[ordre - 1];
  uint tab_dims_j[ordre - 1];

  for (uint i = 0; i < nbelts; i++)
  {
    invphi(ordre, tab_dims, i, indices);
    valeur_i = (*tenseur_vec)[i];
    line_k = indices[k]; // position ligne i_k

    // copie indices et tab_dims et supprime l'element indice par k
    for (int t = 0; t <= ordre - 1; t++)
    {
      if (t < k)
      {
        indices_j[t] = indices[t];
        tab_dims_j[t] = tab_dims[t];
      }
      else
      {
        indices_j[t] = indices[t + 1];
        tab_dims_j[t] = tab_dims[t + 1];
      }
    }

    column_k = phi(ordre - 1, tab_dims_j, indices_j); // position colonne j_k
    T_mode[column_k][line_k] = valeur_i;
  }

  return T_mode;
}

Tenseur pmod(Tenseur tenseur, Matrice matrice, uint k)
{
  uint nbelts = tenseur.get_nbelts();
  int m = matrice.get_nb_lines();
  int n = matrice.get_nb_columns();
  //cout << "n_k: " << n << endl;
  uint ordre = tenseur.get_ordre();
  uint *tab_dims = tenseur.get_tab_dims();
  uint tab_dims_pmod[ordre];
  uint indices_pmod[ordre];
  uint indices[ordre];
  //cout << "m_k: " << m << endl;

  for (int i = 0; i < ordre; i++)
  {
    tab_dims_pmod[i] = (i == k) ? m : tab_dims[i];
  }

  int nbelts_pmod = (nbelts / tab_dims[k]) * m;
  Vecteur tenseur_vec_pmod(nbelts_pmod);

  for (int index = 0; index < nbelts_pmod; index++)
  {
    //cout << "INDEX: " << index << endl;
    invphi(ordre, tab_dims_pmod, index, indices_pmod);
    tenseur_vec_pmod[index] = 0;
    uint i = indices_pmod[k];
    //cout << "I: " << i;

    //cout << " indices: " << endl;
    for (int j = 0; j < ordre; j++)
    {
      indices[j] = indices_pmod[j];
      //cout << indices[j] << ' ';
    }
    //cout << endl;

    for (int j = 0; j < tab_dims[k]; j++)
    {
      indices[k] = j;
      int jj = phi(ordre, tab_dims, indices);
      //cout << " J: " << jj;
      tenseur_vec_pmod[index] = tenseur_vec_pmod[index] + matrice[j][i] * tenseur[phi(ordre, tab_dims, indices)];
    }
    //cout << endl;
  }

  return Tenseur(tenseur_vec_pmod, tab_dims_pmod, ordre);
}
