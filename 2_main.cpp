#include "2_matrice.h"
#include "2_matrice.cpp"
using namespace std;

int main()
{
  float A_c1[] = {1, 1, 0};
  float A_c2[] = {-0.5, 2, -1};
  float A_c3[] = {0, -1, 1};
  float B_c1[] = {-2, 0};
  float B_c2[] = {3, 1};
  float E_diag[] = {3, 2, 1};

  /* 1 */
  cout << "Creation et affichage de A et B" << endl;
  Vecteur A_columns[] = {Vecteur(A_c1, 3), Vecteur(A_c2, 3), Vecteur(A_c3, 3)};
  Vecteur B_columns[] = {Vecteur(B_c1, 2), Vecteur(B_c2, 2)};
  Matrice A(A_columns, 3);
  Matrice B(B_columns, 2);
  A.affiche();
  B.affiche();

  /* 2 */
  cout << "Affichage de B" << endl;
  Matrice C(B);
  B[1][0] = 0; //B[1,2]
  B.affiche();
  cout << "Creation et affichage de C" << endl;

  C.affiche();

  /* 3 */
  cout << "Creation et affichage de D" << endl;
  Matrice D = A.submat(0, 2, 0, 1); //A[1:3][1:2]
  D.affiche();

  /* 4 */
  cout << "Creation et affichage de E" << endl;
  Matrice E(Vecteur(E_diag, 3));
  E.affiche();

  /* 5 */
  cout << "Affichage de B+C" << endl;
  (B + C).affiche();
  cout << "Creation et affichage de C-B" << endl;
  (C - B).affiche();
  cout << "Creation et affichage de D*C" << endl;
  (D * C).affiche();

  /* 6 */
  cout << "Norme de Frobenius pour C" << endl;
  cout << norm(C) << endl;

  /* 7 */
  cout << "Affichage de 0.5*(B+(B)^T)" << endl;
  Matrice F = 0.5 * (B + B.transpose());
  F.affiche();
}
