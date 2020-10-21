#include "1_vecteur.h"
#include "1_vecteur.cpp"
using namespace std;

int main()
{
  float coeff_u[] = {1, 1, 1};
  float coeff_v[] = {3, 4, 0, 0};

  /* 1 */
  cout << "Creation et affichage du vecteur U " << endl;
  Vecteur u(coeff_u, 3);
  u.affiche();
  cout << "Creation et affichage du vecteur V " << endl;

  Vecteur v(coeff_v, 4);
  v.affiche();

  /* 2 */
  cout << "Copie de U dans t" << endl;
  Vecteur t(u);

  /* 3 */
  cout << "Modification et affichage de U" << endl;
  u[2] = 0;
  u.affiche();
  t.affiche();

  /* 4 */
  cout << "Produit scalaire (V)^TV et norme de V" << endl;
  cout << dot(v, v) << endl;
  cout << norm(v) << endl;

  /* 5 */
  cout << "Affichage de (1/norm(V))V" << endl;
  ((1 / norm(v)) * v).affiche();

  /* 6 */
  cout << "Affichage du vecteur V" << endl;
  Vecteur w = v.subvec(1, 3);
  v.affiche();
  cout << "Affichage du vecteur W " << endl;
  w.affiche();

  /* 7 */
  cout << "Affichage de U+W" << endl;
  (u + w).affiche();
  cout << "Affichage de U-W" << endl;
  (u - w).affiche();
}