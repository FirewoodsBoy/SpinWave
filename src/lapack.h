#if !defined(LAPACK_H)
#define LAPACK_H
#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif
#include <iostream>
#include <complex>
#include <memory>
#include <vector>

typedef long int integer;
typedef unsigned long int ulong;
typedef std::complex<double> Complex;
const double EPS = 1e-9;
namespace lapack
{
    extern "C"
    {
        int zheev_(char *jobz, char *uplo, integer *n, Complex *a, integer *lda, double *w, Complex *work, integer *lwork, double *rwork, integer *info);
        int zgemv_(char *trans, integer *m, integer *n, Complex *alpha, Complex *a, integer *lda, Complex *x, integer *incx, Complex *beta, Complex *y, integer *incy);
        int zgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, Complex *alpha, Complex *a, integer *lda, Complex *b, integer *ldb, Complex *beta, Complex *c__, integer *ldc);
        int dsyev_(char *jobz, char *uplo, integer *n, double *a, integer *lda, double *w, double *work, integer *lwork, integer *info);
        int dgemv_(char *trans, integer *m, integer *n, double *alpha, double *a, integer *lda, double *x, integer *incx, double *beta, double *y, integer *incy);
        int dgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, double *alpha, double *a, integer *lda, double *b, integer *ldb, double *beta, double *c__, integer *ldc);
        int zpotrf_(char *uplo, integer *n, Complex *a, integer *lda, integer *info);
        int dpotrf_(char *uplo, integer *n, double *a, integer *lda, integer *info);
        int ztrtri_(char *uplo, char* diag, integer* n, Complex *a, integer* lda, integer *info);
        int dtrtri_(char *uplo, char* diag, integer* n, double *a, integer* lda, integer *info);
        int zpstrf_(char *uplo, integer *n, Complex *a, integer *lda, integer *piv, integer* rank, double* tol, double *work, integer* info);
    }
} // namespace lapack
enum class Trans
{
    N,
    T,
    C
};
char transToChar(Trans t);
int eigenForHermitian(Complex *H, size_t N, double *eigen, bool ifArray = 0, bool upDn = 1);
void complexMv(Complex *A, Complex *v1, Complex *v2, size_t N, Trans trans = Trans::N, Complex alpha = {1, 0}, Complex beta = {0, 0});
void complexMult(Complex *A, Complex *B, Complex *C, size_t N, Trans trans1 = Trans::N, Trans trans2 = Trans::N, Complex alpha = {1, 0}, Complex beta = {0, 0});
int eigenForSym(double *H, size_t N, double *eigen, bool ifArray = 0, bool upDn = 1);
void realMv(double *A, double *v1, double *v2, size_t N, Trans trans = Trans::N, double alpha = 1, double beta = 0);
void realMult(double *A, double *B, double *C, size_t N, Trans trans1 = Trans::N, Trans trans2 = Trans::N);
int cholesky(Complex *H, integer N, bool upLo = true, bool ifZero = true);
int cholesky(double *H, integer N, bool upLo = true, bool ifZero = true);
int triInv(Complex *K, integer N, bool upLo = true);
int triInv(double *K, integer N, bool upLo = true);
std::unique_ptr<double[]> eig(double *H, size_t N, bool ifArray = 0);
std::unique_ptr<double[]> eig(Complex *H, size_t N, bool ifArray = 0);
std::unique_ptr<Complex[]> matMult(Complex *A, Complex *B, size_t N, Trans trans1 = Trans::N, Trans trans2 = Trans::N);
std::unique_ptr<double[]> matMult(double *A, double *B, size_t N, Trans trans1 = Trans::N, Trans trans2 = Trans::N);
/***********************Auxiliary functions*************************/
// std::ostream &operator<<(std::ostream &os, const Complex &a);
// {
//     os << a.real() << '+' << a.imag() << 'i';
//     return os;
// }
template <typename T>
// This function print a matrix to ostream;
// For lapack, MAT[i+1] plus a row with mat[i], thus MATij(row i, col j) should be MAT[i+j*M];
// When trans default, print matrix can directly be used by matlab.
void matPrint(const T *MAT, size_t M, size_t N = 0, std::ostream &os = std::cout, char split = ' ', bool trans = 0)
{
    N = (N == 0) ? M : N;
    if (trans)
    {
        for (size_t i = 0; i < M; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                os << MAT[i * N + j] << split;
            }
            os << '\n';
        }
    }
    else
    {
        for (size_t i = 0; i < M; i++)
        {

            for (size_t j = 0; j < N; j++)
            {
                os << MAT[i + j * M] << split;
            }
            os << '\n';
        }
    }
}
template <typename T>
void arrPrint(const T *arr, size_t l, std::ostream &os = std::cout, char split = ' ')
{
    for (int i = 0; i < l; i++)
    {
        os << arr[i] << split;
    }
}
template <typename T>
void chop(T *arr, size_t l)
{
    for (int i = 0; i < l; i++)
    {
        if (std::abs(arr[i]) < EPS)
        {
            arr[i] = 0;
        }
    }
}
void chop(Complex *arr, size_t l);
template <typename T>
std::unique_ptr<T[]> copy(const T *arr, size_t l)
{
    if (!arr || l <= 0)
    {
        return nullptr;
    }
    auto copyOfArr = std::unique_ptr<T[]>(new T[l]);
    std::copy(arr, arr + l, copyOfArr.get());
    return copyOfArr;
}
#endif // LAPACK_H