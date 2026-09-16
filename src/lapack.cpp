#include "lapack.h"
#include <stdexcept>

char transToChar(Trans t)
{
    switch (t)
    {
    case Trans::N:
        return 'N';
    case Trans::T:
        return 'T';
    case Trans::C:
        return 'C';
    default:
        return 'N';
    }
}

int eigenForHermitian(Complex *H, size_t N, double *eigen, bool ifArray, bool upDn)
{
    /**************************************
     * up=true down=false
     * when ifArray=true, H becomes eigenvectors.
     * Note that H is always changed.
     * ************************************/
    integer n = N;
    integer info = 0;
    char jobz = ifArray ? 'V' : 'N';
    char uplo = upDn ? 'U' : 'L';

    integer lwork_query = -1;
    Complex work_query;
    double *rwork = new double[3 * n - 2];
    std::unique_ptr<double[]> rwork_holder(rwork);

    lapack::zheev_(&jobz, &uplo, &n, H, &n, eigen, &work_query, &lwork_query, rwork, &info);
    if (info != 0)
    {
        throw std::runtime_error("eigenForHermitian failed with info = " + std::to_string(info));
    }
    integer lwork_optimal = static_cast<integer>(work_query.real());
    if (lwork_optimal < 1)
        lwork_optimal = 1;
    auto work = std::unique_ptr<Complex[]>(new Complex[lwork_optimal]);
    lapack::zheev_(&jobz, &uplo, &n, H, &n, eigen, work.get(), &lwork_optimal, rwork, &info);
    return info;
}

void complexMv(Complex *A, Complex *v1, Complex *v2, size_t N, Trans trans, Complex alpha, Complex beta)
{
    /*************************************
     * if trans=Trans::N, A=A;
     * if trans=Trans::T, A=tansverse(A);
     * if trans=Trans::C, A=h.c.(A); (for real condition, T=C)
     * trans1 for A, trans2 for B;
     * ***********************************/
    char t = transToChar(trans);
    integer n = N;
    integer incx = 1, incy = 1;
    lapack::zgemv_(&t, &n, &n, &alpha, A, &n, v1, &incx, &beta, v2, &incy);
}

void complexMult(Complex *A, Complex *B, Complex *C, size_t N, Trans trans1, Trans trans2, Complex alpha, Complex beta)
{
    /*************************************
     * if trans=Trans::N, A=A;
     * if trans=Trans::T, A=tansverse(A);
     * if trans=Trans::C, A=h.c.(A); (for real condition, T=C)
     * trans1 for A, trans2 for B;
     * ***********************************/
    char t1 = transToChar(trans1);
    char t2 = transToChar(trans2);
    integer n = N;
    lapack::zgemm_(&t1, &t2, &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
}

int eigenForSym(double *H, size_t N, double *eigen, bool ifArray, bool upDn)
{
    /**************************************
     * up=true down=false
     * when ifArray=true, H becomes eigenvectors.
     * Note that H is always changed.
     * ************************************/
    integer n = N;
    integer info = 0;
    char jobz = ifArray ? 'V' : 'N';
    char uplo = upDn ? 'U' : 'L';

    // 查询最优工作空间大小
    integer lwork_query = -1;
    double work_query;
    lapack::dsyev_(&jobz, &uplo, &n, H, &n, eigen, &work_query, &lwork_query, &info);
    if (info != 0)
    {
        throw std::runtime_error("eigenForHermitian failed with info = " + std::to_string(info));
    }
    integer lwork_optimal = static_cast<integer>(work_query);
    if (lwork_optimal < 1)
        lwork_optimal = 1;
    auto work = std::unique_ptr<double[]>(new double[lwork_optimal]);
    lapack::dsyev_(&jobz, &uplo, &n, H, &n, eigen, work.get(), &lwork_optimal, &info);
    return info;
}

void realMv(double *A, double *v1, double *v2, size_t N, Trans trans, double alpha, double beta)
{
    /*************************************
     * if trans=Trans::N, A=A;
     * if trans=Trans::T, A=tansverse(A);
     * if trans=Trans::C, A=h.c.(A); (for real condition, T=C)
     * trans1 for A, trans2 for B;
     * ***********************************/
    char t = transToChar(trans);
    integer n = N;
    integer incx = 1, incy = 1;
    lapack::dgemv_(&t, &n, &n, &alpha, A, &n, v1, &incx, &beta, v2, &incy);
}

void realMult(double *A, double *B, double *C, size_t N, Trans trans1, Trans trans2)
{
    /*************************************
     * if trans=Trans::N, A=A;
     * if trans=Trans::T, A=tansverse(A);
     * if trans=Trans::C, A=h.c.(A); (for real condition, T=C)
     * trans1 for A, trans2 for B;
     * ***********************************/
    char t1 = transToChar(trans1);
    char t2 = transToChar(trans2);
    integer n = N;
    double alpha = 1;
    double beta = 0;
    lapack::dgemm_(&t1, &t2, &n, &n, &n, &alpha, A, &n, B, &n, &beta, C, &n);
}

int cholesky(Complex *H, integer N, bool upLo, bool ifZero)
{
    char uplo = upLo ? 'U' : 'L';

    if (ifZero)
    {
        if (upLo)
        {
            for (int i = 0; i < N - 1; i++)
            {
                std::fill(H + i * N + i + 1, H + i * N + N, 0.0);
            }
        }
        else
        {
            for (int i = 1; i < N; i++)
            {
                std::fill(H + i * N, H + i * N + i, 0.0);
            }
        }
    }
    integer info;
    lapack::zpotrf_(&uplo, &N, H, &N, &info);
    return info;
}

int cholesky(double *H, integer N, bool upLo, bool ifZero)
{
    char uplo = upLo ? 'U' : 'L';

    if (ifZero)
    {
        if (upLo)
        {
            for (int i = 0; i < N - 1; i++)
            {
                std::fill(H + i * N + i + 1, H + i * N + N, 0.0);
            }
        }
        else
        {
            for (int i = 1; i < N; i++)
            {
                std::fill(H + i * N, H + i * N + i, 0.0);
            }
        }
    }
    integer info;
    lapack::dpotrf_(&uplo, &N, H, &N, &info);
    return info;
}

int triInv(Complex *K, integer N, bool upLo)
{
    char uplo = upLo ? 'U' : 'L';
    char diag = 'N';
    integer info;
    lapack::ztrtri_(&uplo, &diag, &N, K, &N, &info);
    return info;
}

int triInv(double *K, integer N, bool upLo)
{
    char uplo = upLo ? 'U' : 'L';
    char diag = 'N';
    integer info;
    lapack::dtrtri_(&uplo, &diag, &N, K, &N, &info);
    return info;
}

std::unique_ptr<double[]> eig(double *H, size_t N, bool ifArray)
{
    auto eigens = std::unique_ptr<double[]>(new double[N]);
    eigenForSym(H, N, eigens.get(), ifArray);
    return eigens;
}

std::unique_ptr<double[]> eig(Complex *H, size_t N, bool ifArray)
{
    auto eigens = std::unique_ptr<double[]>(new double[N]);
    eigenForHermitian(H, N, eigens.get(), ifArray);
    return eigens;
}

std::unique_ptr<Complex[]> matMult(Complex *A, Complex *B, size_t N, Trans trans1, Trans trans2)
{
    std::unique_ptr<Complex[]> C = std::unique_ptr<Complex[]>(new Complex[N * N]);
    complexMult(A, B, C.get(), N, trans1, trans2);
    return C;
}

std::unique_ptr<double[]> matMult(double *A, double *B, size_t N, Trans trans1, Trans trans2)
{
    std::unique_ptr<double[]> C = std::unique_ptr<double[]>(new double[N * N]);
    realMult(A, B, C.get(), N, trans1, trans2);
    return C;
}

void chop(Complex *arr, size_t l)
{
    for (size_t i = 0; i < l; i++)
    {
        if (std::abs(arr[i].real()) < EPS)
        {
            arr[i].real(0);
        }
        if (std::abs(arr[i].imag()) < EPS)
        {
            arr[i].imag(0);
        }
    }
}
