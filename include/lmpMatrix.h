#ifndef _MMmatrix_H_
#define _MMmatrix_H_

#include <assert.h>
#include <errno.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifndef DIM
#define DIM 3
#endif

/*** when use template we put the clarification
 * and defination together in the header file  ***/
using std::vector;
namespace MMatrix {
template <class T>
class Mmatrix {
 public:
  Mmatrix(){};
  ~Mmatrix(){};

  // memory allocate
  static T* allocVec(const int length);
  static T** allocMtx(const int rows, const int colms);

  inline static void freeVec(T* vector);
  inline static void freeMtx(const int rows, T** matrix);

  // initialization
  inline static void initVec(const int n, T* vector);
  inline static void initMtx(const int rows, const int colms, T** matrix);
  inline static void initMtxI(const int rows, const int colms, T** matrix);
  inline static void initMtx33(T[DIM][DIM]);
  inline static void initMtxI33(T[DIM][DIM]);

  // print
  inline static void printVec(const int n, T* vector);
  inline static void printVec(const vector<T> vector);
  inline static void printMtx(const int rows, const int colms, T** matrix);
  inline static void printMtx3x3(T mat[3][3]);
  inline static void printMtx2x2(T mat[2][2]);
  inline static void printMtx3x3(const char filename[], const char head[],
                                 T a[3][3]);
  inline static void printVec6(const char filename[], const char head[],
                               T a[6]);

  // vector-vector
  inline static void scalarDotVec(T a, vector<T> v, vector<T>& b);
  inline static double vecInnProd(const vector<T>& a, const vector<T>& b);
  inline static double vecInnProd33(const vector<T>& a, const vector<T>& b);
  inline static double vecInnProd33(T a[3], T b[3]);
  inline static double calNorm(const int dim, std::vector<T> vect);
  inline static void vecOutProd(vector<T> a, vector<T> b, T** M);
  inline static void crossProd(vector<T> a, vector<T> b, vector<T>& res);
  inline static void crossProd33(const T a[3], const T b[3], T res[3]);
  inline static void crossProd33(const vector<T>& a, const vector<T>& b,
                                 vector<T>& res);

  // vector-matrix
  inline static void mtxDotVec(T** M, vector<T> x, vector<T>& b);

  // matrix-matrix
  inline static void mtxAddMtx(int row, int colm, T** A, T** B, T** C);
  inline static void mtxMultmtx(const int rowA, const int clnA, const int clnB,
                                T** A, T** B, T** C);

  // fixed sized matrix-matrix
  inline static void mtx33Multmtx33(T a[3][3], T b[3][3], T c[3][3]);
  inline static void mtx22Multmtx22(T a[2][2], T b[2][2], T c[2][2]);
  inline static void mtx33copymtx33(T a[3][3], T b[3][3]);
  inline static void mtx33Addmtx33(T a[3][3], T b[3][3], T c[3][3]);
  inline static void mtx22Addmtx22(T a[2][2], T b[2][2], T c[2][2]);
  inline static int mtx33Inverse(T a[3][3], T b[3][3]);

  inline static void transPoseMtx33(T a[3][3], T b[3][3]);
  inline static void dim2outprod(T av[2], T bv[2], T mtx[2][2]);
};

/********************************************************
 * matrix multiplication C = A x B
 * *******************************************************/
template <class T>
inline void Mmatrix<T>::mtxMultmtx(const int rowA, const int clnA,
                                   const int clnB, T** A, T** B, T** C) {
  for (int i = 0; i < rowA; i++) {
    for (int j = 0; j < clnB; j++) {
      C[i][j] = 0.0;
      for (int k = 0; k < clnA; k++) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

/********************************************************
 * perform outer product btw two vectors <a|b>
 * *******************************************************/
template <class T>
inline double Mmatrix<T>::vecInnProd(const vector<T>& a, const vector<T>& b) {
  const int n = a.size();
  double res = 0;
  for (int i = 0; i < n; i++) res += a[i] * b[i];
  return res;
}

template <class T>
inline double Mmatrix<T>::vecInnProd33(const vector<T>& a, const vector<T>& b) {
  double res = 0;
  for (int i = 0; i < 3; i++) res += a[i] * b[i];
  return res;
}

template <class T>
inline double Mmatrix<T>::vecInnProd33(T a[3], T b[3]) {
  double res = 0;
  for (int i = 0; i < 3; i++) res += a[i] * b[i];
  return res;
}

/********************************************************
 * perform cross product btw two vectors cross(a, b, c)
 * *******************************************************/
template <class T>
inline void Mmatrix<T>::crossProd(vector<T> a, vector<T> b, vector<T>& res) {
  if (a.size() == 3) {
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
    return;
  }
}

template <class T>
inline void Mmatrix<T>::crossProd33(const T a[3], const T b[3], T res[3]) {
  res[0] = a[1] * b[2] - a[2] * b[1];
  res[1] = a[2] * b[0] - a[0] * b[2];
  res[2] = a[0] * b[1] - a[1] * b[0];
}

template <class T>
inline void Mmatrix<T>::crossProd33(const vector<T>& a, const vector<T>& b,
                                    vector<T>& res) {
  res[0] = a[1] * b[2] - a[2] * b[1];
  res[1] = a[2] * b[0] - a[0] * b[2];
  res[2] = a[0] * b[1] - a[1] * b[0];
}

/********************************************************
 * perform outer product btw two vectors a x b = M
 * *******************************************************/
template <class T>
inline void Mmatrix<T>::vecOutProd(vector<T> a, vector<T> b, T** M) {
  const int N = int(a.size());
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      M[i][j] = a[i] * b[j];
    }
  }
}

/********************************************************
 * perform matrix vector calculation b = Ax
 * *******************************************************/
template <class T>
inline void Mmatrix<T>::mtxDotVec(T** M, vector<T> x, vector<T>& b) {
  const int N = int(x.size());
  for (int i = 0; i < N; i++) {
    b[i] = 0;
    for (int j = 0; j < N; j++) {
      b[i] += M[i][j] * x[j];
    }
  }
  return;
}

/********************************************************
 * perform scalar times vector a x = b
 * *******************************************************/
template <class T>
inline void Mmatrix<T>::scalarDotVec(T a, vector<T> v, vector<T>& b) {
  const int N = int(v.size());
  for (int i = 0; i < N; i++) b[i] = a * v[i];
}

/********************************************************
 * matrix add matrix  A + B = C
 * *******************************************************/
template <class T>
inline void Mmatrix<T>::mtxAddMtx(int row, int colm, T** A, T** B, T** C) {
  for (int i = 0; i < row; i++)
    for (int j = 0; j < colm; j++) C[i][j] = A[i][j] + B[i][j];
}

/********************************************************
 * calculate the norm of vector
 * *******************************************************/
template <class T>
inline double Mmatrix<T>::calNorm(const int dim, std::vector<T> vect) {
  T sum = 0;
  for (int i = 0; i < dim; i++) sum += vect[i] * vect[i];
  return sum;
}

template <class T>
T* Mmatrix<T>::allocVec(const int length) {
  T* vector;
  vector = (T*)calloc(1, length * sizeof(T));

  if (vector != (T*)NULL) {
    return vector;
  } else {
    fprintf(stderr, "vector allocation : %d\n", errno);
    exit(1);
  }
}

template <class T>
void Mmatrix<T>::initVec(const int n, T* vector) {
  printf("The size n is %d\n", n);
  for (int i = 0; i < n; i++) vector[i] = 0.0;
}

template <class T>
void Mmatrix<T>::freeVec(T* vector) {
  assert(vector != (T*)NULL);
  if (vector != (T*)NULL) {
    free(vector);
  }
  return;
}

template <class T>
void Mmatrix<T>::printVec(const int n, T* vector) {
  printf("The size n is %d\n", n);
  for (int i = 0; i < n; i++) {
    std::cout << vector[i] << std::endl;
  }
  return;
}

template <class T>
void Mmatrix<T>::printVec(std::vector<T> vector) {
  const int n = vector.size();
  printf("The size n is %d\n", n);
  for (int i = 0; i < n; i++) {
    std::cout << vector[i] << std::endl;
  }
  return;
}

template <class T>
T** Mmatrix<T>::allocMtx(const int rows, const int colms) {
  T** matrix;
  matrix = new T*[rows];
  for (int i = 0; i < rows; i++) {
    matrix[i] = new T[colms];
    assert(matrix[i] != (T*)NULL);
  }
  if (matrix != (T**)NULL) {
    return matrix;
  } else {
    fprintf(stderr, "matrix allocation : %d\n", errno);
    abort();
  }
}

template <class T>
void Mmatrix<T>::initMtx(const int rows, const int colms, T** matrix) {
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < colms; j++) matrix[i][j] = 0.0;
}

template <class T>
void Mmatrix<T>::initMtxI(const int rows, const int colms, T** matrix) {
  initMtx(rows, colms, matrix);
  for (int i = 0; i < colms; i++) matrix[i][i] = 1.0;
}

template <class T>
void Mmatrix<T>::freeMtx(const int rows, T** matrix) {
  if (matrix != (T**)NULL) {
    for (int i = 0; i < rows; i++) {
      delete[] matrix[i];
    }
    delete[] matrix;
  }
}

template <class T>
void Mmatrix<T>::printMtx(const int rows, const int colms, T** matrix) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < colms; j++) {
      // std::cout<<i<<" "<<j<<" "<<matrix[i][j]<<"\t";
      fprintf(stdout, "M[%d][%d] %.4f", i, j, matrix[i][j]);
    }
    std::cout << std::endl;
  }
}

/*-------------------------------------------------------------------------
 *      Function:     Matrix33Mult33
 *      Description:  Multiplies two 3X3 matrices
 *
 *      Arguments:
 *          a    3 X 3 array containing components of the first matrix
 *          b    3 X 3 array containing components of the second matrix
 *          c    3 X 3 array in which to return to the caller the
 *               resulting matrix after multiplying matrix <a> by
 *               matrix <b>.
 *------------------------------------------------------------------------*/

template <class T>
void Mmatrix<T>::mtx33Multmtx33(T a[3][3], T b[3][3], T c[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      c[i][j] = 0.0;
      for (int k = 0; k < 3; k++) c[i][j] += a[i][k] * b[k][j];
    }
}

template <class T>
void Mmatrix<T>::mtx22Multmtx22(T a[2][2], T b[2][2], T c[2][2]) {
  const int dim = 2;
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++) {
      c[i][j] = 0.0;
      for (int k = 0; k < dim; k++) c[i][j] += a[i][k] * b[k][j];
    }
}

/*-------------------------------------------------------------------------
 *      Function:     Matrix33Mult33
 *      Description:  copy 3x3 matrix b to 3x3 matrix a
 *
 *      Arguments:
 *               3 x 3 matrix a
 *               3 x 3 matrix b
 *------------------------------------------------------------------------*/
template <class T>
void Mmatrix<T>::mtx33copymtx33(T a[3][3], T b[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) a[i][j] = b[i][j];
}

template <class T>
void Mmatrix<T>::initMtx33(T mtx[DIM][DIM]) {
  for (int i = 0; i < DIM; i++)
    for (int j = 0; j < DIM; j++) mtx[i][j] = 0.0;
}

template <class T>
void Mmatrix<T>::initMtxI33(T mtx[DIM][DIM]) {
  for (int i = 0; i < DIM; i++)
    for (int j = 0; j < DIM; j++) mtx[i][j] = 0.0;

  mtx[0][0] = 1.0;
  mtx[1][1] = 1.0;
  mtx[2][2] = 1.0;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Transpose
 *      Description:  Transpose rows and columns of a 3 X 3 matrix
 *
 *      Arguments:
 *          mat    3 X 3 source matrix to be transposed.
 *          trans  3 X 3 matrix in which to return to the caller the
 *                 transpose of matrix <mat>.
 *
 *------------------------------------------------------------------------*/
template <class T>
void Mmatrix<T>::transPoseMtx33(T mat[3][3], T trans[3][3]) {
  int m, n;
  for (m = 0; m < 3; m++)
    for (n = 0; n < 3; n++) trans[n][m] = mat[m][n];
}

template <class T>
void Mmatrix<T>::mtx33Addmtx33(T a[3][3], T b[3][3], T c[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) c[i][j] = a[i][j] + b[i][j];
}

template <class T>
void Mmatrix<T>::mtx22Addmtx22(T a[2][2], T b[2][2], T c[2][2]) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) c[i][j] = a[i][j] + b[i][j];
}

template <class T>
void Mmatrix<T>::printMtx3x3(T a[3][3]) {
  for (int i = 0; i < DIM; i++) {
    printf("[ %.5f , %.5f  ,%.5f ]\n", a[i][0], a[i][1], a[i][2]);
  }
}

template <class T>
void Mmatrix<T>::printMtx2x2(T a[2][2]) {
  for (int i = 0; i < 2; i++) printf("[ %.5f , %.5f]\n", a[i][0], a[i][1]);
}

template <class T>
void Mmatrix<T>::printMtx3x3(const char filename[], const char head[],
                             T a[3][3]) {
  FILE* ptFile;
  if ((ptFile = fopen(filename, "a")) != (NULL)) {
    fprintf(ptFile, "%s\n", head);
    for (int i = 0; i < DIM; i++) {
      fprintf(ptFile, "[ %.5f , %.5f  ,%.5f ]\n", a[i][0], a[i][1], a[i][2]);
    }
    fflush(ptFile);
    fclose(ptFile);
  } else {
    fprintf(stderr, "can not open %s\n", filename);
    exit(1);
  }
}

template <class T>
void Mmatrix<T>::printVec6(const char filename[], const char head[], T a[6]) {
  FILE* ptFile;
  if ((ptFile = fopen(filename, "a")) != (NULL)) {
    fprintf(ptFile, "%s\n", head);
    fprintf(ptFile,
            "[ %.5f , %.5f  ,%.5f ]\n"
            "[ %.5f , %.5f  ,%.5f ]\n",
            a[0], a[1], a[2], a[3], a[4], a[5]);
    fflush(ptFile);
    fclose(ptFile);
  } else {
    fprintf(stderr, "can not open %s\n", filename);
    exit(1);
  }
}

/*
 *      (When the numbers involved were on the order of 10e-12 or
 *      smaller)
 *
 *      Returns:  1 on success
 *                0 if the matrix was not invertible
 */
template <class T>
int Mmatrix<T>::mtx33Inverse(T a[3][3], T b[3][3]) {
  int i, j, k;
  T p, fmax, fval, eps = 1.0e-20;
  T tmp[3][3];

  /*
   *      Initialize the inverse to the identity matrix and create
   *      a temporary copy of the source matrix.
   */
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      b[i][j] = (T)(i == j);
      tmp[i][j] = a[i][j];
    }
  }

  for (i = 0; i < 3; i++) {
    fmax = fabs(tmp[i][i]);
    for (j = i + 1; j < 3; j++) {
      if (fabs(tmp[j][i]) > fmax) {
        fmax = fabs(tmp[j][i]);
        for (k = 0; k < 3; k++) {
          p = tmp[i][k];
          tmp[i][k] = tmp[j][k];
          tmp[j][k] = p;
          p = b[i][k];
          b[i][k] = b[j][k];
          b[j][k] = p;
        }
      }
    }
    /*
     *          If can't do the inversion, return 0
     */
    if (fmax < eps) {
      printf("Matrix33Invert: fmax < eps, cannot invert!\n");
      return (0);
    }

    fval = 1.0 / tmp[i][i];

    for (j = 0; j < 3; j++) {
      tmp[i][j] *= fval;
      b[i][j] *= fval;
    }

    for (k = 0; k < 3; k++) {
      if (k != i) {
        fval = tmp[k][i];
        for (j = 0; j < 3; j++) {
          tmp[k][j] -= fval * tmp[i][j];
          b[k][j] -= fval * b[i][j];
        }
      }
    }
  }
  return (1);
}

/*** calculate dim2 vector out product, the result is mtx[2][2] ***/
template <class T>
void Mmatrix<T>::dim2outprod(T va[2], T vb[2], T mtx[2][2]) {
  mtx[0][0] = va[0] * vb[0];
  mtx[0][1] = va[0] * vb[1];
  mtx[1][0] = va[1] * vb[0];
  mtx[1][1] = va[1] * vb[1];
}

}  // namespace MMatrix
#endif /**myMatrix_H_ **/
