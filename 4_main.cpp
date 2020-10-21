#include "4_tenseur.cpp"
#include "4_tenseur.h"
using namespace std;

int main()
{
  uint dims_T[] = {2, 2, 2};
  uint taille_T = 3;

  /* 1 */
  cout << "Creation du tenseur T" << endl;
  Tenseur T(dims_T, taille_T);
  T.affiche();

  uint dims_U[] = {2, 2, 2};
  uint taille_U = 3;

  /* 2 */

  cout << "Creation du tenseur U" << endl;
  Tenseur U(dims_U, taille_U);
  for (uint i = 0; i < T.get_nbelts(); i++)
  {
    U[i] = 1;
  }
  U.affiche();

  /* 3 */
  cout << "Creation du tenseur V=U+T" << endl;
  Tenseur V(U + T);
  V.affiche();

  cout << "Creation du tenseur W=U-T" << endl;
  Tenseur W(U - T);
  W.affiche();

  /* 4 */
  cout << "Affichage de U" << endl;
  U.affiche();
  cout << "Modification de U" << endl;
  for (uint i = 0; i < U.get_nbelts(); i++)
  {
    U[i] = U[i] * (-1);
  }
  U.affiche();
  cout << "Affichage de W" << endl;
  W.affiche();

  /* 5 */
  cout << "Mode 1 de T" << endl;
  (T.mode(0)).affiche();

  /* 6 */
  float mode2_c1[] = {1, 4};
  float mode2_c2[] = {3, 1. / 3};
  float mode2_c3[] = {0, 3. / 2};
  float mode2_c4[] = {-1, 2};
  Vecteur mode2_columns[] = {Vecteur(mode2_c1, 2), Vecteur(mode2_c2, 2),
                             Vecteur(mode2_c3, 2), Vecteur(mode2_c4, 2)};
  Matrice mode2(mode2_columns, 4);
  T = Tenseur(mode2, dims_T, 3, 1);
  cout << "Modification de T" << endl;
  T.affiche();
  cout << "Mode 2 de T"<< endl;
  (T.mode(1)).affiche();

  /* 7 */
  float A_c1[] = {3, 0, 0};
  float A_c2[] = {-1, 6, -3};
  Vecteur A_columns[] = {Vecteur(A_c1, 3), Vecteur(A_c2, 3)};
  Matrice A(A_columns, 2);
  cout << " --- Matrice A --- " << endl;
  A.affiche();

  cout << "Produit modal de T et A : " << endl;
  Tenseur S(pmod(T, A, 2));
  S.affiche();
  cout << "Mode 3 de S" << endl;
  (S.mode(2)).affiche();

  /* 8 */
  Tenseur R(S + S);
  cout << "Affichage du tenseur R  =  S + S" << endl;
  R.affiche();
}
