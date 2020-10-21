#include "1_vecteur.h"
#include "2_matrice.cpp"
#include "2_matrice.h"
#include "3_matrice_svd.cpp"
#include "3_matrice_svd.h"

using namespace std;

int main()
{
  float z_1 = 1;
  float z_2 = 2;
  float c;
  float s;

  /* 1 */
  cout << "Resultat de 'givens' pour z_1 = 2 et z_2 = 1" << endl;
  cppidcolgivens(z_1, z_2, c, s);
  cout << "Valeur de s : " << s << endl;
  cout << "Valeur de c : " << c << endl;

  /* 2 */
  float beta_x;
  float beta_y;
  float beta_z;
  float coeff_x[] = {-1., 0};
  float coeff_y[] = {1. / sqrt(2), 1. / sqrt(2)};
  float coeff_z[] = {-4.};
  Vecteur x(coeff_x, 2);
  Vecteur y(coeff_y, 2);
  Vecteur z(coeff_z, 1);

  cout << "Resultat de 'householder' pour x" << endl;
  householder(x, beta_x).affiche();
  cout << "Valeur de beta : " << beta_x << endl;

  cout << "Resultat de 'householder' pour y" << endl;
  householder(y, beta_y).affiche();
  cout << "Valeur de beta : " << beta_y << endl;

  cout << "Resultat de 'householder' pour z" << endl;
  householder(z, beta_z).affiche();
  cout << "Valeur de beta : " << beta_z << endl;

  /* 3 */

  float M_c1[] = {10, -6};
  float M_c2[] = {-6, 10};
  Vecteur M_columns[] = {Vecteur(M_c1, 2), Vecteur(M_c2, 2)};
  Matrice M(M_columns, 2);
  Matrice Q(2, 2);

  cout << " --- Matrice M --- " << endl;
  M.affiche();
  cout << "Resultat pour QR Symetrique appliquee a M " << endl;
  qrsym(M, Q);
  setPrecision(Q);
  Q.affiche();

  /* 
  
  float N_c1[] = {10, -6, 0};
  float N_c2[] = {-6, 10, 0};
  float N_c3[] = {0, 0, 1};
  Vecteur N_columns[] = {Vecteur(N_c1, 3), Vecteur(N_c2, 3), Vecteur(N_c3, 3)};
  Matrice N(N_columns, 3);
  Matrice R(3, 3);

  cout << " --- Matrice N (bonus pour la validation) --- " << endl;
  N.affiche();
  cout << "Resultat pour QR Symetrique appliquee a N " << endl;
  qrsym(N, R);
  setPrecision(R);
  R.affiche();

  */

  /* 4 */

  float O_c1[] = {1, 1. / 2, 1. / 3};
  float O_c2[] = {1. / 2, 1. / 3, 1. / 4};
  float O_c3[] = {1. / 3, 1. / 4, 1. / 5};
  Vecteur O_columns[] = {Vecteur(O_c1, 3), Vecteur(O_c2, 3), Vecteur(O_c3, 3)};
  Matrice O(O_columns, 3);
  Matrice O_copy(O);
  Matrice S(3, 3);

  cout << " --- Matrice O (M dans le sujet) --- " << endl;
  O.affiche();
  cout << "Resultat pour QR Symetrique avec Pivot appliquee a O " << endl;
  Matrice PI = qrpivot(O, S);
  cout << "Matrice PI" << endl;
  setPrecision(PI);
  PI.affiche();
  cout << "Matrice Q" << endl;
  setPrecision(S);
  S.affiche();
  cout << "Matrice R" << endl;
  Matrice R_1 = S.transpose() * O_copy * PI;
  setPrecision(R_1);
  R_1.affiche();

  /* 5 */
  float A_c1[] = {1, 0};
  float A_c2[] = {0, -1};
  Vecteur A_columns[] = {Vecteur(A_c1, 2), Vecteur(A_c2, 2)};
  Matrice A(A_columns, 2);
  cout << " --- Matrice A --- " << endl;
  A.affiche();
  uint m_A = A.get_nb_lines();
  uint n_A = A.get_nb_columns();
  Matrice U_A(m_A, m_A);
  Matrice SIGMA_A(m_A, n_A);
  Matrice V_A(n_A, n_A);
  Matrice decompositions_A[3] = {U_A, SIGMA_A, V_A};
  svd(A, decompositions_A);
  cout << "Decomposition en valeurs singulieres de la matrice A" << endl;
  cout << "Matrice U" << endl;
  decompositions_A[0].affiche();
  cout << "Matrice SIGMA" << endl;
  decompositions_A[1].affiche();
  cout << "Matrice V" << endl;
  decompositions_A[2].affiche();

  float B_c1[] = {2. * sqrt(2), -sqrt(2)};
  float B_c2[] = {-2. * sqrt(2), -sqrt(2)};
  Vecteur B_columns[] = {Vecteur(B_c1, 2), Vecteur(B_c2, 2)};
  Matrice B(B_columns, 2);
  cout << " --- Matrice B --- " << endl;
  B.affiche();
  uint m_B = B.get_nb_lines();
  uint n_B = B.get_nb_columns();
  Matrice U_B(m_B, m_B);
  Matrice SIGMA_B(m_B, n_B);
  Matrice V_B(n_B, n_B);
  Matrice decompositions_B[3] = {U_B, SIGMA_B, V_B};
  svd(B, decompositions_B);
  cout << "Decomposition en valeurs singulieres de la matrice B" << endl;
  cout << "Matrice U" << endl;
  setPrecision(decompositions_B[0]);
  decompositions_B[0].affiche();
  cout << "Matrice SIGMA" << endl;
  setPrecision(decompositions_B[1]);
  decompositions_B[1].affiche();
  cout << "Matrice V" << endl;
  setPrecision(decompositions_B[2]);
  decompositions_B[2].affiche();

  float C_c1[] = {1. / 2, sqrt(3) / 2};
  float C_c2[] = {3 * sqrt(3) / 2, -3. / 2};
  float C_c3[] = {0, 0};
  Vecteur C_columns[] = {Vecteur(C_c1, 2), Vecteur(C_c2, 2), Vecteur(C_c3, 2)};
  Matrice C(C_columns, 3);
  setPrecision(C);
  cout << " --- Matrice C --- " << endl;
  C.affiche();
  uint m_C = C.get_nb_lines();
  uint n_C = C.get_nb_columns();
  Matrice U_C(m_C, m_C);
  Matrice SIGMA_C(m_C, n_C);
  Matrice V_C(n_C, n_C);
  Matrice decompositions_C[3] = {U_C, SIGMA_C, V_C};
  svd(C, decompositions_C);
  cout << "Decomposition en valeurs singulieres de la matrice C" << endl;
  cout << "Matrice U" << endl;
  setPrecision(decompositions_C[0]);

  decompositions_C[0].affiche();
  cout << "Matrice SIGMA" << endl;
  setPrecision(decompositions_C[1]);
  decompositions_C[1].affiche();
  cout << "Matrice V" << endl;
  setPrecision(decompositions_C[2]);
  decompositions_C[2].affiche();
}
