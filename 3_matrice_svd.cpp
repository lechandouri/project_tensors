#include "3_matrice_svd.h"
#include "1_vecteur.h"
#include "2_matrice.h"

using namespace std;

void cppidcolgivens(float x, float z, float &c, float &s)
{
  if (z == 0)
  {
    c = 1;
    s = 0;
    return;
  }

  float tau;

  if (abs(z) > abs(x))
  {
    tau = -x / z;
    s = 1 / (sqrt(1 + pow(tau, 2)));
    c = s * tau;
  }
  else
  {
    tau = -z / x;
    c = 1 / (sqrt(1 + pow(tau, 2)));
    s = c * tau;
  }
}

Vecteur householder(Vecteur x, float &beta)
{
  float sigma;

  if (x.get_dim() == 1) // cas particulier d'un vecteur de dimension 1
  {
    sigma = 0;
  }
  else
  {
    Vecteur x2n = x.subvec(1, x.get_dim() - 1);
    sigma = dot(x2n, x2n);
  }

  Vecteur v(x);
  v[0] = 1;

  if (sigma == 0)
  {
    if (x[0] >= 0)
      beta = 0;
    else
      beta = 2;
    return v;
  }

  float mu = sqrt(pow(x[0], 2) + sigma);

  if (x[0] <= 0)
  {
    v[0] = x[0] - mu;
  }
  else
  {
    v[0] = -sigma / (x[0] + mu);
  }

  beta = 2 * pow(v[0], 2) / (sigma + pow(v[0], 2));

  return (1 / v[0]) * v;
}

Matrice reductridiag(Matrice &D)
{
  uint l = D.get_nb_columns() - 1;
  float d = (D[l - 1][l - 1] - D[l][l]) / 2;
  int signe = d >= 0 ? 1 : -1;
  float mu = D[l][l] - pow(D[l - 1][l], 2) /
                           (d + signe * sqrt(pow(d, 2) + pow(D[l - 1][l], 2)));
  float x = D[0][0] - mu;
  float z = D[0][1];
  float c;
  float s;
  float tau_1;
  float tau_2;
  Matrice Z(l + 1, l + 1);

  for (uint i = 0; i <= l; i++)
  {
    Z[i][i] = 1;
  }

  for (int k = 0; k <= l - 1; k++)
  {
    cppidcolgivens(x, z, c, s);

    for (int j = 0; j <= l; j++)
    {
      tau_1 = D[k][j];
      tau_2 = D[k + 1][j];
      D[k][j] = c * tau_1 - s * tau_2;
      D[k + 1][j] = s * tau_1 + c * tau_2;
      tau_1 = Z[k][j];
      tau_2 = Z[k + 1][j];
      Z[k][j] = c * tau_1 - s * tau_2;
      Z[k + 1][j] = s * tau_1 + c * tau_2;
    }

    for (int j = 0; j <= l; j++)
    {
      tau_1 = D[j][k];
      tau_2 = D[j][k + 1];
      D[j][k] = c * tau_1 - s * tau_2;
      D[j][k + 1] = s * tau_1 + c * tau_2;
    }

    if (k < l - 1)
    {
      x = D[k][k + 1];
      z = D[k][k + 2];
    }
  }
  return Z;
}

void qrsym(Matrice &A, Matrice &Q)
{
  int n = A.get_nb_columns() - 1; // n est la derniere colonne de A
  float beta;

  for (int i = 0; i <= n; i++)
  {
    for (int j = 0; j <= n; j++)
    {
      Q[j][i] = (i == j) ? 1 : 0;
    }
  }

  for (int k = 0; k <= n - 2; k++)
  {
    Vecteur v = householder(A[k].subvec(k + 1, n), beta);
    Vecteur p = beta * (A.submat(k + 1, n, k + 1, n).mvprod(v));
    Vecteur w = p - (beta / 2) * (dot(p, v) * v);
    A[k][k + 1] = norm(A[k].subvec(k + 1, n));
    A[k + 1][k] = A[k][k + 1];

    Matrice Q_tmp = beta * (outer(v, v) * Q.submat(k + 1, n, k + 1, n));

    for (int i = k + 1; i <= n; i++)
    {
      for (int j = k + 1; j <= n; j++)
      {
        A[j][i] = A[j][i] - v[i - (k + 1)] * w[j - (k + 1)] -
                  w[i - (k + 1)] * v[j - (k + 1)];
        Q[j][i] = Q[j][i] - Q_tmp[j - (k + 1)][i - (k + 1)];
      }
    }
  }

  Matrice T(n + 1, n + 1);
  for (int j = n; j >= 0; j--)
  {
    T[j][j] = A[j][j];
    if (j > 0)
    {
      T[j][j - 1] = A[j - 1][j];
      T[j - 1][j] = T[j][j - 1];
    }
  }

  while (!isDiagonal(T))
  {
    for (int i = 0; i <= n - 1; i++)
    {
      if (abs(T[i + 1][i]) + abs(T[i][i + 1]) <=
          pow(10, -9) * (abs(T[i][i]) + abs(T[i + 1][i + 1])))
      {
        T[i + 1][i] = 0;
        T[i][i + 1] = 0;
      }
    }

    int p;
    int q;

    for (p = 0; p <= n; p++)
    {
      if (p == 0)
      {
        if (!isEqualZero(T[p][p + 1]))
          break;
      }
      else if (p < n)
      {
        if (!(isEqualZero(T[p][p + 1]) && isEqualZero(T[p][p - 1])))
        {
          break;
        }
      }
    }

    for (q = 0; q <= n; q++)
    {
      if (q == 0)
      {
        if (!isEqualZero(T[n - q][(n - q) - 1]))
          break;
      }
      else if (q < n)
      {
        if (!(isEqualZero(T[n - q][n - q + 1]) &&
              isEqualZero(T[n - q][n - q - 1])))
        {
          {
            break;
          }
        }
      }
    }

    Matrice T_2 = T.submat(p, n - q, p, n - q);

    if (p + q < n + 1)
    {
      Matrice Z(reductridiag(T_2));
      Matrice T_hat(T);

      Matrice Q_hat(Q);
      for (uint i = p; i <= n - q; i++)
      {
        for (uint j = p; j <= n - q; j++)
        {
          T_hat[j][i] = T_2[j - p][i - p];
          Q_hat[j][i] = Z[j - p][i - p];
        }
      }

      T = (1. / 2) * (T_hat + T_hat.transpose());
      Q = Q * Q_hat;
    }
  }
  Matrice Q_Q_T(Q * (Q.transpose()));
  // remplacer valeurs quasi-nulles par zero
  setPrecision(Q_Q_T);
}

Matrice qrpivot(Matrice &A, Matrice &Q)
{
  int m = A.get_nb_lines() - 1; // indice de la derniere ligne de A
  int n = A.get_nb_columns() - 1;
  Matrice PI(n + 1, n + 1);
  float c[n + 1];
  float beta;

  for (int j = 0; j <= n; j++)
  {
    PI[j][j] = 1; // PI Matrice identite
    c[j] = dot(A[j], A[j]);
  }

  int k;
  int r = -1;
  float tau = get_maximum(c, n + 1);

  while (tau > 0 && r < n)
  {
    r = r + 1;

    for (k = r; k <= n; k++)
    {
      if (isEqualZero(c[k] - tau))
      {
        break;
      }
    }

    // permutation des vecteurs extraits de A
    float permut_vec_1;
    for (int i = 0; i <= m; i++)
    {
      permut_vec_1 = A[r][i];
      A[r][i] = A[k][i];
      A[k][i] = permut_vec_1;
    }

    // permutation des reels de c
    float permut_sca = c[r];
    c[r] = c[k];
    c[k] = permut_sca;

    // permutation des vecteurs extraits de PI
    float permut_vec_2;
    for (int i = 0; i <= n; i++)
    {
      permut_vec_2 = PI[r][i];
      PI[r][i] = PI[k][i];
      PI[k][i] = permut_vec_2;
    }

    Vecteur v = householder(A[r].subvec(r, m), beta);
    Matrice beta_vvt_A(beta * outer(v, v) * (A.submat(r, m, r, n)));
    for (int i = r; i <= m; i++)
    {
      for (int j = r; j <= n; j++)
      {
        A[j][i] = A[j][i] - beta_vvt_A[j - r][i - r];
      }
    }

    for (int i = r + 1; i <= m; i++)
    {
      A[r][i] = v[i - r];
    }

    for (int i = r + 1; i <= n; i++)
    {
      c[i] = c[i] - pow(A[i][r], 2);
    }

    if (r < n)
    {
      float c_bis[n - r];
      for (int i = r + 1; i <= n; i++)
      {
        c_bis[i - r - 1] = c[i]; // {c_r+1, ..., c_n}
      }
      tau = get_maximum(c_bis, n - r);
    }
    else
    {
      tau = 0;
    }
  }
  // Q matrice identite
  for (int i = 0; i <= m; i++)
  {
    for (int j = 0; j <= m; j++)
    {
      Q[j][i] = (i == j) ? 1. : 0;
    }
  }

  Vecteur v_(m + 1);
  for (int j = n; j >= 0; j--)
  {
    v_[j] = 1;
    if (j + 1 <= m)
    {
      for (int l = j + 1; l <= m; l++)
      {
        v_[l] = A[j][l];
      }

      beta = 2.0 / (1 + pow(norm(A[j].subvec(j + 1, m)), 2));
    }

    else // cas particulier où m=n et donc j+1=m+1>m
    {
      beta = 2;
    }

    Matrice beta_vvt_Q(beta * outer(v_.subvec(j, m), v_.subvec(j, m)) *
                       (Q.submat(j, m, j, m)));

    for (int s = j; s <= m; s++)
    {
      for (int t = j; t <= m; t++)
      {
        Q[t][s] = Q[t][s] - beta_vvt_Q[t - j][s - j];
      }
    }
  }
  return PI;
}

void svd(Matrice A, Matrice decompositions[3])
{
  uint m = A.get_nb_lines() - 1;
  uint n = A.get_nb_columns() - 1;

  if (m >= n)
  {
    Matrice Q_1(n + 1, n + 1);
    Matrice Q_2(m + 1, m + 1);
    Matrice A_T_A(A.transpose() * A);
    qrsym(A_T_A, Q_1);
    Matrice A_Q_1(A * Q_1);
    Matrice A_Q_1_copy(A_Q_1);
    Matrice PI = qrpivot(A_Q_1, Q_2);
    Matrice R(Q_2.transpose() * A_Q_1_copy * PI);

    for (int j = 0; j <= n; j++)
    {
      if (R[j][j] < 0)
      {
        for (int i = 0; i <= m; i++)
        {
          Q_2[j][i] = -Q_2[j][i];
        }
      }
    }

    R = Q_2.transpose() * A_Q_1_copy * PI;

    decompositions[0] = Q_2;      // U
    decompositions[1] = R;        // SIGMA
    decompositions[2] = Q_1 * PI; // V
  }
  else
  {
    Matrice Q_1(m + 1, m + 1);
    Matrice Q_2(n + 1, n + 1);
    Matrice A_A_T(A * A.transpose());
    qrsym(A_A_T, Q_1);
    Matrice A_T_Q_1(A.transpose() * Q_1);
    Matrice A_T_Q_1_copy(A_T_Q_1);
    Matrice PI = qrpivot(A_T_Q_1, Q_2);
    Matrice R(Q_2.transpose() * A_T_Q_1_copy * PI);

    for (int i = 0; i <= m; i++)
    {
      if (R[i][i] < 0)
      {
        for (int j = 0; j <= n; j++)
        {
          Q_2[i][j] = -Q_2[i][j];
        }
      }
    }

    R = Q_2.transpose() * A_T_Q_1_copy * PI;

    decompositions[0] = Q_1 * PI;      // U
    decompositions[1] = R.transpose(); // SIGMA
    decompositions[2] = Q_2;           // V
  }
}

bool isEqualZero(float f) { return abs(f) < pow(10, -6); }

bool isDiagonal(Matrice &A)
{
  for (uint i = 0; i < A.get_nb_lines(); i++)
    for (uint j = 0; j < A.get_nb_columns(); j++)
      if (i != j && !isEqualZero(A[j][i]))
      {
        return false;
      }
  return true;
}

float get_maximum(float *liste, int dim)
{
  float maxi = liste[0];
  for (int i = 0; i < dim; i++)
  {
    if (liste[i] > maxi)
    {
      maxi = liste[i];
    }
  }
  return maxi;
}
