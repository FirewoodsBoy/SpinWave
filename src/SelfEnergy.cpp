#include "SelfEnergy.h"

std::vector<Complex> oneLoop(Vec2<size_t> q, size_t type, const Bosons &B3, const SpinWave &sw,
                             Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum)
{
    double delta = (omegaRange[1] - omegaRange[0]) / omegaNum;
    std::vector<Complex> sigma(omegaNum + 1);
    std::vector<Complex> dos(omegaNum + 1);
    size_t Nm = sw.getNm();
    size_t N = sw.getN();
    size_t q1 = q[0], q2 = q[1];
    double N2 = std::pow(N, 2);
    if (B3.getOrder(3).size())
    {
        for (size_t t1 = 0; t1 < Nm; t1++)
        {
            for (size_t t2 = t1; t2 < Nm; t2++)
            {
                Bosons B1f = B3.grep({{type, 1}, {t1, 0}, {t2, 0}});
                // Bosons B2f = B3.grep({{t1, 1}, {t2, 1}, {type, 0}});
                std::vector<boson> backward1 = {{type, 1}, {t1, 1}, {t2, 1}};
                std::sort(backward1.begin(), backward1.end(), [](const boson &b1, const boson &b2)
                          { return b1.t < b2.t; });
                Bosons B1b = B3.grep(backward1);
                // std::vector<boson> backward2 = {{type, 0}, {t1, 0}, {t2, 0}};
                // std::sort(backward2.begin(), backward2.end(), [](const boson &b1, const boson &b2)
                //           { return b1.t < b2.t; });
                // Bosons B2b = B3.grep(backward2);
                size_t T{}, I, J;
                for (size_t i = 0; i < 3; i++)
                {
                    if (backward1[i].t == type)
                    {
                        T = i;
                        break;
                    }
                }
                I = (T + 1) % 3;
                J = (T + 2) % 3;

                std::vector<std::vector<Complex>> sigmas;
                sigmas.reserve(threadNum);
                // std::vector<std::vector<Complex>> doss;
                // doss.reserve(threadNum);
                for (size_t i = 0; i < threadNum; i++)
                {
                    sigmas.emplace_back(omegaNum + 1);
                    // doss.emplace_back(omegaNum + 1);
                }

#pragma omp parallel num_threads(threadNum)
                {
                    size_t myRank = omp_get_thread_num();
                    size_t kq1, kq2;
                    double omega;
                    Complex V1, V2, Vf, Vb;
                    for (size_t k1 = myRank; k1 < N; k1 += threadNum)
                    {
                        for (size_t k2 = 0; k2 < N; k2++)
                        {
                            kq1 = (2 * N + N / 2 - q1 - k1) % N;
                            kq2 = (2 * N + N / 2 - q2 - k2) % N;
                            std::vector<Vec2<size_t>> ks = {{q1, q2}, {k1, k2}, {kq1, kq2}};
                            // forward scattering
                            if (t1 == t2)
                            {
                                V1 = vertex(B1f, ks);
                                std::swap(ks[1], ks[2]);
                                V1 += vertex(B1f, ks);
                                V2 = std::conj(V1);
                                Vf = V1 * V2 / 2.;
                            }
                            else
                            {
                                V1 = vertex(B1f, ks);
                                V2 = std::conj(V1);
                                Vf = V1 * V2;
                            }
                            // backward scattering
                            if (backward1[T].t == backward1[J].t)
                            {
                                V1 = 0;
                                std::vector<size_t> perm = {0, 1, 2};
                                do
                                {
                                    ks[perm[0]] = {q1, q2};
                                    ks[perm[1]] = {k1, k2};
                                    ks[perm[2]] = {kq1, kq2};
                                    V1 += vertex(B1b, ks);
                                } while (std::next_permutation(perm.begin(), perm.end()));
                                V2 = std::conj(V1);
                                Vb = V1 * V2 / 2.;
                            }
                            else
                            {
                                ks[T] = {q1, q2};
                                ks[I] = {k1, k2};
                                ks[J] = {kq1, kq2};
                                V1 = vertex(B1b, ks);
                                V2 = std::conj(V1);
                                Vb = V1 * V2;
                                if (backward1[I].t == backward1[J].t)
                                {
                                    std::swap(ks[I], ks[J]);
                                    V1 += vertex(B1b, ks);
                                    V2 = std::conj(V1);
                                    Vb = V1 * V2 / 2.;
                                }
                                else if (backward1[I].t == backward1[T].t)
                                {
                                    std::swap(ks[I], ks[T]);
                                    V1 += vertex(B1b, ks);
                                    V2 = std::conj(V1);
                                    Vb = V1 * V2;
                                }
                            }
                            Complex cf = sw.getEk(t1, k1, k2) + sw.getEk(t2, kq1, kq2) - eps;
                            Complex cb = sw.getEk(backward1[I].t, N - k1, N - k2) + sw.getEk(backward1[J].t, N - kq1, N - kq2) - eps;
                            for (size_t o = 0; o < omegaNum + 1; o++)
                            {
                                omega = o * delta + omegaRange[0];
                                // Complex d = 1. / (omega - cf);
                                // doss[myRank][o] += d;
                                sigmas[myRank][o] += Vf / (omega - cf);
                                sigmas[myRank][o] += Vb / (-omega - cb);
                            }
                        }
                    }
                }
                // parallel finish

                for (size_t i = 0; i < threadNum; i++)
                {
                    for (size_t o = 0; o < omegaNum + 1; o++)
                    {
                        sigma[o] += sigmas[i][o] / N2;
                        // dos[o] += doss[i][o] / N2;
                    }
                }
            }
        }
    }
    return sigma;
}
std::vector<Complex> twoLoop(Vec2<size_t> q, size_t type, const Bosons &B4, const SpinWave &sw,
                             Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum)
{
    double delta = (omegaRange[1] - omegaRange[0]) / omegaNum;
    std::vector<Complex> sigma(omegaNum + 1);
    std::vector<Complex> dos(omegaNum + 1);
    size_t Nm = sw.getNm();
    size_t N = sw.getN();
    size_t q1 = q[0], q2 = q[1];
    double N4 = std::pow(N, 4);
    double omega{};
    // Zero loop scattering (flavor change process)
    // Forward scattering
    for (size_t t1 = 0; t1 < Nm; t1++)
    {
        if (t1 != type)
        {
            Bosons Bf = B4.grep({{type, true}, {t1, false}});
            Complex V1 = vertex(Bf, {{q1, q2}, {N - q1, N - q2}});
            Complex V = V1 * std::conj(V1);
            for (size_t o = 0; o < omegaNum; o++)
            {
                omega = o * delta + omegaRange[0];
                sigma[o] += V / (omega - sw.getEk(t1, q1, q2) + eps);
            }
        }
    }
    // Backward scattering
    for (size_t t1 = 0; t1 < Nm; t1++)
    {
        if (t1 != type)
        {
            Complex V1{}, V_{};
            if (t1 < type)
            {
                Bosons Bb = B4.grep({{t1, true}, {type, true}});
                V1 = vertex(Bb, {{N - q1, N - q2}, {q1, q2}});
                V_ = V1 * std::conj(V1);
            }
            else
            {
                Bosons Bb = B4.grep({{type, true}, {t1, true}});
                V1 = vertex(Bb, {{q1, q2}, {N - q1, N - q2}});
                V_ = V1 * std::conj(V1);
            }
            for (size_t o = 0; o < omegaNum; o++)
            {
                omega = o * delta + omegaRange[0];
                sigma[o] += V_ / (-omega - sw.getEk(t1, q1, q2) + eps);
            }
        }
    }
    if (B4.getOrder(4).size() > 0)
    {
        for (size_t t1 = 0; t1 < Nm; t1++)
        {
            for (size_t t2 = t1; t2 < Nm; t2++)
            {
                for (size_t t3 = t2; t3 < Nm; t3++)
                {
                    // Forward scattering
                    Bosons Bf = B4.grep({{type, true}, {t1, false}, {t2, false}, {t3, false}});
                    Bosons Bf1;
                    if (t1 == t3)
                    {
                        Bf1 = Bf.swap({1, 2, 3});
                        Bf1 *= sqrt(1. / 6.);
                    }
                    else if (t1 == t2)
                    {
                        Bf1 = Bf.swap({1, 2});
                        Bf1 *= sqrt(0.5);
                    }
                    else if (t2 == t3)
                    {
                        Bf1 = Bf.swap({2, 3});
                        Bf1 *= sqrt(0.5);
                    }
                    else
                    {
                        Bf1 = Bf;
                    }
                    Bf.clear();
                    // std::cout<<Bf<<std::endl;

                    std::vector<boson> backward = {{type, true}, {t1, true}, {t2, true}, {t3, true}};
                    std::sort(backward.begin(), backward.end(), [](const boson &b1, const boson &b2)
                              { return b1.t < b2.t; });
                    Bosons Bb = B4.grep(backward);
                    size_t T{}, I, J, K;
                    for (size_t i = 0; i < 4; i++)
                    {
                        if (backward[i].t == type)
                        {
                            T = i;
                            break;
                        }
                    }
                    I = (T + 1) % 4;
                    J = (T + 2) % 4;
                    K = (T + 3) % 4;
                    Bosons Bb1;
                    if (t1 == t3 && t1 == type)
                    {
                        Bb1 = Bb.swap({0,1,2,3});
                        Bb1 *= sqrt(1. / 6.);
                    }else if (backward[T].t == backward[J].t){
                        Bb1 = Bb.swap({T, I, J});
                        Bb1 *= sqrt(0.5);
                    }else if (backward[T].t == backward[I].t){
                        Bb1 = Bb.swap({T, I});
                        if(backward[J].t == backward[K].t){
                            Bb1 += Bb1.swap(J, K);
                            Bb1 *= sqrt(0.5);
                        }
                    }else if ((backward[I].t == backward[K].t)){
                        Bb1 = Bb.swap({I, J, K});
                        Bb1 *= sqrt(1. / 6.);
                    }else if (backward[I].t == backward[J].t){
                        Bb1 = Bb.swap({I, J});
                        Bb1 *= sqrt(0.5);
                    }else if (backward[J].t == backward[K].t){
                        Bb1 = Bb.swap({J, K});
                        Bb1 *= sqrt(0.5);
                    }else{
                        Bb1 = Bb;
                    }
                    Bb.clear();
                    std::vector<std::vector<Complex>> sigmas;
                    sigmas.reserve(threadNum);
                    for (size_t i = 0; i < threadNum; i++)
                    {
                        sigmas.emplace_back(omegaNum + 1);
                    }
#pragma omp parallel num_threads(threadNum)
                    {
                        size_t myRank = omp_get_thread_num();
                        size_t kq1, kq2;
                        double omega;
                        Complex V1, V2, Vf, Vb;
                        for (size_t k1 = myRank; k1 < N; k1 += threadNum)
                        {
                            for (size_t k2 = 0; k2 < N; k2++)
                            {
                                std::vector<Vec2<size_t>> ks = {{q1, q2}, {k1, k2}, {}, {}};
                                auto Bfp = Bf1.reduce(ks, 2, 3);
                                ks[T] = {q1, q2};
                                ks[I] = {k1, k2};
                                auto Bbp = Bb1.reduce(ks, J, K);
                                for (size_t p1 = 0; p1 < N; p1++)
                                {
                                    for (size_t p2 = 0; p2 < N; p2++)
                                    {
                                        kq1 = (-q1 - k1 - p1 + 3 * N) % N;
                                        kq2 = (-q2 - k2 - p2 + 3 * N) % N;
                                        V1 = Bfp[p1 * N + p2];
                                        V2 = std::conj(V1);
                                        Vf = V1 * V2;
                                        // backward scattering
                                        V1 = Bbp[p1 * N + p2];
                                        V2 = std::conj(V1);
                                        Vb = V1 * V2;
                                        // Vf = 1.;
                                        // Vb = 1.;
                                        Complex cf = sw.getEk(t1, k1, k2) + sw.getEk(t2, p1, p2) + sw.getEk(t3, kq1, kq2) - eps;
                                        Complex cb = sw.getEk(backward[I].t, N - k1, N - k2) + sw.getEk(backward[J].t, N - p1, N - p2) + sw.getEk(backward[K].t, N - kq1, N - kq2) - eps;
                                        // std::cout << Vf << ' ' << Vb <<std::endl;
                                        for (size_t o = 0; o < omegaNum + 1; o++)
                                        {
                                            omega = o * delta + omegaRange[0];
                                            // Complex d = 1. / (omega - cf);
                                            // doss[myRank][o] += d;
                                            sigmas[myRank][o] += Vf / (omega - cf);
                                            sigmas[myRank][o] += Vb / (-omega - cb);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // parallel finish

                    for (size_t i = 0; i < threadNum; i++)
                    {
                        for (size_t o = 0; o < omegaNum + 1; o++)
                        {
                            sigma[o] += sigmas[i][o] / N4;
                            // dos[o] += doss[i][o] / N2;
                        }
                    }
                }
            }
        }
    }
    return sigma;
}
std::vector<double> twoMagnonDos(Vec2<size_t> q, const SpinWave &sw,
                                 Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum)
{
    std::vector<double> dos(omegaNum + 1);
    double delta = (omegaRange[1] - omegaRange[0]) / omegaNum;
    size_t N = sw.getN(), Nm = sw.getNm();
    double N2 = std::pow(N, 2);
    size_t q1 = N - q[0], q2 = N - q[1];

    for (size_t t1 = 0; t1 < Nm; t1++)
    {
        for (size_t t2 = t1; t2 < Nm; t2++)
        {
            std::vector<std::vector<double>> doss;
            doss.reserve(threadNum);
            // std::vector<std::vector<Complex>> doss;
            // doss.reserve(threadNum);
            for (size_t i = 0; i < threadNum; i++)
            {
                doss.emplace_back(omegaNum + 1);
                // doss.emplace_back(omegaNum + 1);
            }
            double sym = 0.5;
            if (t1 != t2)
            {
                sym = 1;
            }

#pragma omp parallel num_threads(threadNum)
            {
                size_t myRank = omp_get_thread_num();
                size_t kq1, kq2;
                double omega;
                for (size_t k1 = myRank; k1 < N; k1 += threadNum)
                {
                    for (size_t k2 = 0; k2 < N; k2++)
                    {
                        kq1 = (2 * N + N / 2 - q1 - k1) % N;
                        kq2 = (2 * N + N / 2 - q2 - k2) % N;
                        Complex cf = (sw.getEk(t1, k1, k2) + sw.getEk(t2, kq1, kq2) - eps);
                        for (size_t o = 0; o < omegaNum + 1; o++)
                        {
                            omega = o * delta + omegaRange[0];
                            doss[myRank][o] += -(sym / (omega - cf)).imag() / M_PI;
                            // doss[myRank][o] += sym/(omegaRange[1] - omegaRange[0]);
                        }
                    }
                }
            }
            // parallel finish
            for (size_t i = 0; i < threadNum; i++)
            {
                for (size_t o = 0; o < omegaNum + 1; o++)
                {
                    dos[o] += doss[i][o] / N2;
                }
            }
        }
    }
    return dos;
}
std::vector<double> threeMagnonDos(Vec2<size_t> q, const SpinWave &sw,
                                   Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum)
{
    std::vector<double> dos(omegaNum + 1);
    double delta = (omegaRange[1] - omegaRange[0]) / omegaNum;
    size_t N = sw.getN(), Nm = sw.getNm();
    double N4 = std::pow(N, 4);
    size_t q1 = N - q[0], q2 = N - q[1];

    for (size_t t1 = 0; t1 < Nm; t1++)
    {
        for (size_t t2 = t1; t2 < Nm; t2++)
        {
            for (size_t t3 = t2; t3 < Nm; t3++)
            {
                std::vector<std::vector<double>> doss;
                doss.reserve(threadNum);
                // std::vector<std::vector<Complex>> doss;
                // doss.reserve(threadNum);
                for (size_t i = 0; i < threadNum; i++)
                {
                    doss.emplace_back(omegaNum + 1);
                    // doss.emplace_back(omegaNum + 1);
                }
                double sym;
                if (t1 == t3)
                {
                    sym = 1. / 6.;
                }
                else
                {
                    if (t1 == t2 || t2 == t3)
                    {
                        sym = 0.5;
                    }
                    else
                    {
                        sym = 1;
                    }
                }
#pragma omp parallel num_threads(threadNum)
                {
                    size_t myRank = omp_get_thread_num();
                    size_t kq1, kq2;
                    double omega;
                    for (size_t k1 = myRank; k1 < N; k1 += threadNum)
                    {
                        for (size_t k2 = 0; k2 < N; k2++)
                        {
                            for (size_t p1 = 0; p1 < N; p1++)
                            {
                                for (size_t p2 = 0; p2 < N; p2++)
                                {
                                    kq1 = (-q1 - k1 - p1 + 3 * N) % N;
                                    kq2 = (-q2 - k2 - p2 + 3 * N) % N;
                                    Complex cf = (sw.getEk(t1, k1, k2) + sw.getEk(t2, p1, p2) + sw.getEk(t3, kq1, kq2) - eps);
                                    for (size_t o = 0; o < omegaNum + 1; o++)
                                    {
                                        omega = o * delta + omegaRange[0];
                                        doss[myRank][o] += -(sym / (omega - cf)).imag() / M_PI;
                                    }
                                }
                            }
                        }
                    }
                }
                // parallel finish
                for (size_t i = 0; i < threadNum; i++)
                {
                    for (size_t o = 0; o < omegaNum + 1; o++)
                    {
                        dos[o] += doss[i][o] / N4;
                    }
                }
            }
        }
    }
    return dos;
}
GreensProp oneLoopGreenFunction(Vec2<size_t> q, size_t type, const Bosons &BB, const SpinWave &sw,
                                Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum)
{
    Vec2<size_t> qi{sw.getN() - q[0], sw.getN() - q[1]};
    Bosons B2 = BB.grep({{type, 1}, {type, 0}});
    Complex v2 = vertex(B2, {qi, q});
    std::vector<Complex> green;
    green.reserve(omegaNum + 1);
    double delta = (omegaRange[1] - omegaRange[0]) / omegaNum;
    auto sigma = oneLoop(qi, type, BB.getOrder(3), sw, omegaRange, omegaNum, eps, threadNum);
    double omega = omegaRange[0];
    Complex temp = sw.getEk(type, q[0], q[1]) + v2 - eps;
    for (size_t i = 0; i < omegaNum + 1; i++)
    {
        green.push_back(1. / (omega - temp - sigma[i]));
        omega += delta;
    }
    size_t max_index = std::distance(green.begin(), std::min_element(green.begin(), green.end(), [](Complex a, Complex b)
                                                                     { return a.imag() < b.imag(); }));
    double dispersion = sw.getEk(type, q[0], q[1]) + v2.real() + sigma[max_index].real();
    // For test
    // std::cout << "dispersion - max_index*delta = " << dispersion - (omegaRange[0] + max_index * delta) << std::endl;
    return GreensProp{green, sigma, dispersion};
}
GreensProp twoLoopGreenFunction(Vec2<size_t> q, size_t type, const Bosons &B4, const Bosons &B62, const SpinWave &sw,
                                Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum)
{
    Vec2<size_t> qi{sw.getN() - q[0], sw.getN() - q[1]};
    std::vector<Complex> green;
    green.reserve(omegaNum + 1);
    double delta = (omegaRange[1] - omegaRange[0]) / omegaNum;
    auto sigma = twoLoop(qi, type, B4, sw, omegaRange, omegaNum, eps, threadNum);
    auto B2 = (B4.getOrder(2) + B62.getOrder(2)).grep({{type, 1}, {type, 0}});
    Complex v2 = vertex(B2, {qi, q});
    Complex temp = sw.getEk(type, q[0], q[1]) + v2 - eps;
    double omega = omegaRange[0];
    for (size_t i = 0; i < omegaNum + 1; i++)
    {
        green.push_back(1. / (omega - temp - sigma[i]));
        omega += delta;
    }
    size_t max_index = std::distance(green.begin(), std::min_element(green.begin(), green.end(), [](Complex a, Complex b)
                                                                     { return a.imag() < b.imag(); }));
    double dispersion = sw.getEk(type, q[0], q[1]) + v2.real() + sigma[max_index].real();
    // For test
    // std::cout << "dispersion - max_index*delta = " << dispersion - (omegaRange[0] + max_index * delta) << std::endl;
    return GreensProp{green, sigma, dispersion};
}
void shiftMomentum(double &k)
{
    while (k < -M_PI)
    {
        k += 2 * M_PI;
    }
    while (k >= M_PI)
    {
        k -= 2 * M_PI;
    }
}
GreensProp bilinearInterp(const std::array<GreensProp *, 4> &gps, Vec2<double> p)
{
    GreensProp result;
    size_t size = gps[0][0].green.size();
    result.reserve(size);
    std::vector<Complex>::const_iterator si[4], gi[4];
    for (size_t i = 0; i < 4; i++)
    {
        gi[i] = gps[i][0].green.cbegin();
        si[i] = gps[i][0].sigma.cbegin();
        // di[i] = gps[i][0].dos.cbegin();
    }
    for (size_t i = 0; i < size; i++)
    {
        result.green.push_back(bilinearInterp(std::array<Complex, 4>{*gi[0], *gi[1], *gi[2], *gi[3]}, p));
        result.sigma.push_back(bilinearInterp(std::array<Complex, 4>{*si[0], *si[1], *si[2], *si[3]}, p));
        // result.dos.push_back(bilinearInterp(std::array<Complex, 4>{*di[0], *di[1], *di[2], *di[3]}, p));
        for (size_t j = 0; j < 4; j++)
        {
            gi[j]++;
            si[j]++;
            // di[j]++;
        }
    }
    result.dispersion = bilinearInterp(std::array<double, 4>{gps[0][0].dispersion, gps[1][0].dispersion, gps[2][0].dispersion, gps[3][0].dispersion}, p);
    return result;
}
/*****************************************************************/
// SelfEnergy private methods
double SelfEnergy::length(Vec2<double> a1, Vec2<double> a2)
{
    Vec2<double> a = {a1[0] - a2[0], (a1[1] - a2[1]) * rate};
    double l2 = std::pow(a[0], 2) + std::pow(a[1], 2) + 2 * a[0] * a[1] * cos(theta);
    return std::sqrt(l2);
}

std::array<Vec2<size_t>, 4> SelfEnergy::nearbyPoints(Vec2<double> k, Vec2<double> &p)
{
    shiftMomentum(k[0]);
    shiftMomentum(k[1]);
    k[0] += M_PI;
    k[1] += M_PI;
    double delta = (2 * M_PI) / sw.getN();
    size_t n1 = k[0] / delta;
    size_t n2 = k[1] / delta;
    std::array<Vec2<size_t>, 4> points;
    points[0] = {n1, n2};
    points[1] = {n1 + 1, n2};
    points[2] = {n1, n2 + 1};
    points[3] = {n1 + 1, n2 + 1};
    p = {(k[0] - n1 * delta) / delta, (k[1] - n2 * delta) / delta};
    return points;
}
void SelfEnergy::generateOneLoopInteraction()
{
    if (status < 2)
        checkSpinWave();
    Bosons_ BI = sw.getHPTerm(2) + sw.getHPTerm(3) + sw.getHPTerm(4) + (-1.) * sw.getHmf();
    BI.simplify();
    Bosons Bk(BI);
    Bk.merge();
    BB1loop = sw.BogoliubovTransformation(Bk);
    BB1loop.order();
}
void SelfEnergy::generateTwoLoopInteraction()
{
    if (status < 2)
        checkSpinWave();
    if (!BB1loop.size())
    {
        generateOneLoopInteraction();
    }
    Bosons Bk6 = sw.getHPTerm(6);
    BB6 = sw.BogoliubovTransformation(Bk6);
    BB6.order();
}
void SelfEnergy::generateOneLoopDatas(size_t type)
{
    if (status < 3)
    {
        generatePathPoints(sw.getN());
    }
    if (!BB1loop.size())
    {
        generateOneLoopInteraction();
    }
    if (!gp1loops[type].size())
    {
        gp1loops[type].reserve(nearbys.size());
    }
    std::vector<Vec2<size_t>> nearbysVector;
    nearbysVector.insert(nearbysVector.end(), nearbys.begin(), nearbys.end());
    std::mutex map_mutex;
    size_t n0{}, t0{};
    size_t num = nearbysVector.size();
#pragma omp parallel num_threads(threadNum)
    {
        size_t myRank = omp_get_thread_num();
        for (size_t i = myRank; i < num; i += threadNum)
        {
            if ((myRank == 0) && (n0 == 0) && print)
            {
                std::cout << "calculating selfenergy********" << 0. / nearbysVector.size() << '%' << " with type " << type << " from rank " << myRank << std::endl;
                t0++;
            }
            auto gf = firstOrderGreenFunction(nearbysVector[i], type, 1);
            std::lock_guard<std::mutex> lock(map_mutex);
            gp1loops[type][nearbysVector[i]] = gf;
            n0++;
            if ((10 * double(n0) / num >= t0) && n0 && print)
            {
                // std::lock_guard<std::mutex> lock(map_mutex);
                std::cout << "calculating selfenergy********" << 100 * double(n0) / num << '%' << " with type " << type << " from rank " << myRank << std::endl;
                t0++;
            }
        }
    }
}
void SelfEnergy::generateTwoLoopDatas(size_t type)
{
    if (status < 3)
    {
        generatePathPoints(sw.getN());
    }
    if (!BB6.size())
    {
        generateTwoLoopInteraction();
    }
    if (!gp2loops[type].size())
    {
        gp2loops[type].reserve(nearbys.size());
    }
    std::vector<Vec2<size_t>> nearbysVector;
    nearbysVector.insert(nearbysVector.end(), nearbys.begin(), nearbys.end());
    std::mutex map_mutex;
    size_t n0{}, t0{};
    size_t num = nearbysVector.size();
#pragma omp parallel num_threads(threadNum)
    {
        size_t myRank = omp_get_thread_num();
        for (size_t i = myRank; i < num; i += threadNum)
        {
            if ((myRank == 0) && (n0 == 0) && print)
            {
                std::cout << "calculating selfenergy********" << 0. / nearbysVector.size() << '%' << " with type " << type << " from rank " << myRank << std::endl;
                t0++;
            }
            auto gf = secondOrderGreenFunction(nearbysVector[i], type, 1);
            std::lock_guard<std::mutex> lock(map_mutex);
            gp2loops[type][nearbysVector[i]] = gf;
            n0++;
            if ((100 * double(n0) / num >= t0) && n0 && print)
            {
                // std::lock_guard<std::mutex> lock(map_mutex);
                std::cout << "calculating selfenergy********" << 100 * double(n0) / num << '%' << " with type " << type << " from rank " << myRank << std::endl;
                t0++;
            }
        }
    }
}
void SelfEnergy::generateTwoMagnonDoS()
{
    if (status < 3)
    {
        generatePathPoints(sw.getN());
    }
    std::vector<Vec2<size_t>> nearbysVector;
    nearbysVector.insert(nearbysVector.end(), nearbys.begin(), nearbys.end());
    std::mutex map_mutex;
    size_t n0{}, t0{};
    size_t num = nearbysVector.size();
#pragma omp parallel num_threads(threadNum)
    {
        size_t myRank = omp_get_thread_num();
        for (size_t i = myRank; i < num; i += threadNum)
        {
            if ((myRank == 0) && (n0 == 0) && print)
            {
                std::cout << "calculating two-magnon DoS********" << 0. / nearbysVector.size() << '%' << " from rank " << myRank << std::endl;
                t0++;
            }
            std::vector<double> dos = twoMagnonDos(nearbysVector[i], sw, omegaRange, omegaNum, eps, 1);

            std::lock_guard<std::mutex> lock(map_mutex);
            twoMagnonContinuumData[nearbysVector[i]] = dos;
            n0++;
            if ((10 * double(n0) / num >= t0) && n0 && print)
            {
                // std::lock_guard<std::mutex> lock(map_mutex);
                std::cout << "calculating two-magnon DoS********" << 100 * double(n0) / num << '%' << " from rank " << myRank << std::endl;
                t0++;
            }
        }
    }
}
void SelfEnergy::generateThreeMagnonDoS()
{
    if (status < 3)
    {
        generatePathPoints(sw.getN());
    }
    std::vector<Vec2<size_t>> nearbysVector;
    nearbysVector.insert(nearbysVector.end(), nearbys.begin(), nearbys.end());
    std::mutex map_mutex;
    size_t n0{}, t0{};
    size_t num = nearbysVector.size();
#pragma omp parallel num_threads(threadNum)
    {
        size_t myRank = omp_get_thread_num();
        for (size_t i = myRank; i < num; i += threadNum)
        {
            if ((myRank == 0) && (n0 == 0) && print)
            {
                std::cout << "calculating three-magnon DoS********" << 0. / nearbysVector.size() << '%' << " from rank " << myRank << std::endl;
                t0++;
            }
            std::vector<double> dos = threeMagnonDos(nearbysVector[i], sw, omegaRange, omegaNum, eps, 1);
            std::lock_guard<std::mutex> lock(map_mutex);
            threeMagnonContinuumData[nearbysVector[i]] = dos;
            n0++;
            if ((100 * double(n0) / num >= t0) && n0 && print)
            {
                // std::lock_guard<std::mutex> lock(map_mutex);
                std::cout << "calculating three-magnon DoS********" << 100 * double(n0) / num << '%' << " from rank " << myRank << std::endl;
                t0++;
            }
        }
    }
}
std::vector<GreensProp> SelfEnergy::getInterpValue(Data<GreensProp> &data) const
{
    std::vector<GreensProp> interpValue;
    interpValue.reserve(pathPoints.size());
    for (auto p : pathPoints)
    {
        const auto &nearby = pointsNearby.at(p);
        std::array<GreensProp *, 4> temp = {&data[nearby[0]], &data[nearby[1]], &data[nearby[2]], &data[nearby[3]]};
        interpValue.push_back(bilinearInterp(temp, pointsPosition.at(p)));
    }
    return interpValue;
}
std::vector<std::vector<Complex>> SelfEnergy::dynamicStructureFactor(std::vector<Data<std::vector<Complex>>> &green,
                                                                     const std::vector<Vec2<double>> &sublatticePosition, Vec3<double> weight, int higgs) const
{
    size_t Nm = sw.getNm();
    size_t N = sw.getN();
    std::vector<std::vector<Complex>> structure;
    std::unordered_map<Vec2<double>, std::vector<Complex>> sData;
    structure.reserve(pathPoints.size());
    sData.reserve(nearbysTrueCoor.size());
    for (auto point : nearbysTrueCoor)
    {
        sData[point.first] = std::vector<Complex>(omegaNum + 1);
    }
    if (higgs < 2)
    {
        for (int s = 0; s < 3; s++)
        {
            Bosons_ B;
            for (size_t t = 0; t < Nm; t++)
            {
                B = B + std::sqrt(weight[s]) * sw.SToB(spin{t, s}, rotation(sw.getMagneticOrder()[t]));
            }
            Bosons_ B1 = B.getOrder(1);
            B1.simplify();
            Bosons Bk(B1, sublatticePosition);
            Bosons BB = sw.BogoliubovTransformation(Bk);
            for (size_t t = 0; t < Nm; t++)
            {
                Bosons Bv = BB.grep({boson{t, false}});
                Bosons Bdv = BB.grep({boson{t, true}});
                for (auto point : nearbysTrueCoor)
                {
                    Complex V1 = vertex(Bv, {point.second}, point.first);
                    Complex V2 = vertex(Bdv, {{N - point.second[0], N - point.second[1]}},
                                        {-point.first[0], -point.first[1]});
                    if (std::abs(V1 - std::conj(V2)) > 1e-6)
                    {
                        std::cout << V1 << ' ' << V2 << std::endl;
                    }
                    for (size_t o = 0; o < omegaNum + 1; o++)
                    {
                        sData[point.first][o] += V1 * V2 * green[t][point.second][o];
                    }
                }
            }
        }
    }
    if (higgs)
    {
        double delta = (omegaRange[1] - omegaRange[0]) / omegaNum;
        size_t count = 0;
        size_t N = sw.getN();
        double N2 = std::pow(N, 2);
        for (int s = 0; s < 3; s++)
        {
            Bosons_ B;
            for (size_t t = 0; t < Nm; t++)
            {
                B = B + std::sqrt(weight[s]) * sw.SToB(spin{t, s}, rotation(sw.getMagneticOrder()[t]));
            }
            Bosons_ B2 = B.getOrder(2);
            B2.simplify();
            Bosons Bk(B2, sublatticePosition);
            Bosons BB = sw.BogoliubovTransformation(Bk);
            BB.order();
            size_t total = 3 * sw.getNm() * (sw.getNm() + 1) / 2 * nearbysTrueCoor.size();
            for (size_t t1 = 0; t1 < sw.getNm(); t1++)
            {
                for (size_t t2 = t1; t2 < sw.getNm(); t2++)
                {
                    Bosons Bv = BB.grep({boson{t1, false}, boson{t2, false}});
                    Bosons Bdv = BB.grep({boson{t1, true}, boson{t2, true}});
                    for (auto point : nearbysTrueCoor)
                    {
                        if (count % (total / 10) == 0)
                        {
                            std::cout << "Calculating 4 points correlation (Higgs) " << std::ceil(100. * count / total) << '%' << std::endl;
                        }
                        count++;
                        Vec2<size_t> q = point.second;
                        std::vector<std::vector<Complex>> as;
                        as.reserve(threadNum);
                        for (size_t rank = 0; rank < threadNum; rank++)
                        {
                            as.emplace_back(omegaNum + 1);
                        }
#pragma omp parallel num_threads(threadNum)
                        {
                            size_t myRank = omp_get_thread_num();
                            size_t kq1{}, kq2{};
                            Complex V1{}, V2{};
                            for (size_t k1 = myRank; k1 < N; k1 += threadNum)
                            {
                                for (size_t k2 = 0; k2 < N; k2++)
                                {
                                    kq1 = (2 * N + N / 2 + q[0] - k1) % N;
                                    kq2 = (2 * N + N / 2 + q[1] - k2) % N;
                                    V1 = vertex(Bv, {{k1, k2}, {kq1, kq2}}, point.first);
                                    V2 = vertex(Bdv, {{N - k1, N - k2}, {N - kq1, N - kq2}}, {-point.first[0], -point.first[1]});
                                    Complex V = V1 * V2;
                                    if (t1 == t2)
                                    {
                                        V1 += vertex(Bv, {{kq1, kq2}, {k1, k2}}, point.first);
                                        V2 += vertex(Bdv, {{N - kq1, N - kq2}, {N - k1, N - k2}}, {-point.first[0], -point.first[1]});
                                        V = V1 * V2 / 2.;
                                    }
                                    double omega = omegaRange[0];
                                    Complex temp = sw.getEk(t1, k1, k2) + sw.getEk(t2, kq1, kq2) - eps;
                                    for (size_t o = 0; o < omegaNum + 1; o++)
                                    {
                                        as[myRank][o] += V / (omega - temp);
                                        omega += delta;
                                    }
                                }
                            }
                        }
                        // parallel finished
                        for (size_t rank = 0; rank < threadNum; rank++)
                        {
                            for (size_t o = 0; o < omegaNum + 1; o++)
                            {
                                sData[point.first][o] += as[rank][o] / N2;
                                // std::cout << as[rank][o] / N2 << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    for (auto p : pathPoints)
    {
        auto trueCoor = pointsNearbyTrueCoor.at(p);
        structure.push_back(bilinearInterp(std::array<std::vector<Complex> *, 4>{
                                               &sData[trueCoor[0]],
                                               &sData[trueCoor[1]],
                                               &sData[trueCoor[2]],
                                               &sData[trueCoor[3]]},
                                           pointsPosition.at(p)));
    }
    return structure;
}
/********************************************************************************/
// SelfEnergy public methods
SelfEnergy::SelfEnergy(const SpinWave &sw) : sw(sw)
{
    omegaNum = sw.getN();
    gp1loops.resize(sw.getNm());
    gp2loops.resize(sw.getNm());
}

SelfEnergy::SelfEnergy(const SpinWave &sw, const std::vector<Vec2<double>> &path) : sw(sw), path(path)
{
    omegaNum = sw.getN();
    gp1loops.resize(sw.getNm());
    gp2loops.resize(sw.getNm());
}

SelfEnergy::SelfEnergy(const SpinWave &sw, const std::vector<Vec2<double>> &path, double theta, double rate) : SelfEnergy(sw, path)
{
    this->theta = theta;
    this->rate = rate;
}

void SelfEnergy::generatePathPoints(size_t n)
{
    if (status < 2)
    {
        checkSpinWave();
    }
    double delta = 2 * M_PI / sw.getN();
    pointsIndex.clear();
    pathPoints.clear();
    pointsNearby.clear();
    pointsNearbyTrueCoor.clear();
    nearbys.clear();
    nearbysTrueCoor.clear();
    pointsIndex.reserve(path.size());
    std::vector<double> ls;
    ls.reserve(path.size() - 1);
    size_t nt = 0;
    pointsIndex.push_back(nt);
    Vec2<double> pp;
    for (auto p = path.begin(); p < path.end() - 1; p++)
    {
        double l = length(*p, *(p + 1));
        ls.push_back(l);
        size_t n0 = l / M_PI * n;
        nt += n0;
        pointsIndex.push_back(nt);
    }
    pathPoints.reserve(nt + 1);
    pointsNearby.reserve(nt + 1);
    for (size_t i = 0; i < path.size() - 1; i++)
    {
        size_t nl = pointsIndex[i + 1] - pointsIndex[i];
        double deltaX = (path[i + 1].x - path[i].x) / nl;
        double deltaY = (path[i + 1].y - path[i].y) / nl;
        Vec2<double> p = path[i];
        for (size_t j = 0; j < nl; j++)
        {
            pathPoints.push_back(p);
            Vec2<long long> index = {static_cast<long long>(std::floor((p[0] + M_PI) / delta)),
                                     static_cast<long long>(std::floor((p[1] + M_PI) / delta))};
            auto nearby = nearbyPoints(p, pp);
            std::array<Vec2<double>, 4> nearbyTrueCoor = {
                Vec2<double>{-M_PI + delta * index[0], -M_PI + delta * index[1]},
                Vec2<double>{-M_PI + delta * (index[0] + 1), -M_PI + delta * index[1]},
                Vec2<double>{-M_PI + delta * index[0], -M_PI + delta * (index[1] + 1)},
                Vec2<double>{-M_PI + delta * (index[0] + 1), -M_PI + delta * (index[1] + 1)}};
            pointsNearby[p] = nearby;
            pointsNearbyTrueCoor[p] = nearbyTrueCoor;
            for (size_t k = 0; k < 4; k++)
            {
                nearbys.insert(nearby[k]);
                nearbysTrueCoor[nearbyTrueCoor[k]] = nearby[k];
            }
            pointsPosition[p] = pp;
            p[0] += deltaX;
            p[1] += deltaY;
        }
    }
    auto p = path.back();
    pathPoints.push_back(path.back());
    auto nearby = nearbyPoints(path.back(), pp);
    pointsNearby[path.back()] = nearby;
    pointsPosition[path.back()] = pp;
    Vec2<long long> index = {static_cast<long long>(std::floor((p[0] + M_PI) / delta)),
                             static_cast<long long>(std::floor((p[1] + M_PI) / delta))};
    std::array<Vec2<double>, 4> nearbyTrueCoor = {
        Vec2<double>{-M_PI + delta * index[0], -M_PI + delta * index[1]},
        Vec2<double>{-M_PI + delta * (index[0] + 1), -M_PI + delta * index[1]},
        Vec2<double>{-M_PI + delta * index[0], -M_PI + delta * (index[1] + 1)},
        Vec2<double>{-M_PI + delta * (index[0] + 1), -M_PI + delta * (index[1] + 1)}};
    pointsNearbyTrueCoor[p] = nearbyTrueCoor;
    for (size_t i = 0; i < 4; i++)
    {
        nearbys.insert(nearby[i]);
        nearbysTrueCoor[nearbyTrueCoor[i]] = nearby[i];
    }
    status = 3;
}

std::vector<double> SelfEnergy::getHmfDispersion(size_t type)
{
    if (sw.getStatus() < 3)
    {
        throw std::runtime_error("BogoliubovMatrix() of SpinWave obj must be called before getting disperion.");
    }
    if (status < 3)
    {
        generatePathPoints(sw.getN());
    }
    Data<double> ek;
    ek.reserve(nearbys.size());
    for (auto nearby : nearbys)
    {
        ek[nearby] = sw.getEk(type, nearby[0], nearby[1]);
    }
    return getInterpValue(ek);
}

void SelfEnergy::checkSpinWave()
{
    if (sw.getStatus() < 2)
    {
        throw std::runtime_error("BogoliubovMatrix() of SpinWave obj must be called before selfenergy calculation.");
    }
    double maxEs{};
    if (eps == Complex{} || omegaRange == Vec2<double>{})
    {
        maxEs = sw.maxEs();
    }
    if (eps == Complex{})
    {
        eps.imag(3 * maxEs / sw.getN());
    }
    if (omegaRange == Vec2<double>{})
    {
        omegaRange = {0, 1.5 * maxEs};
    }
    status = 2;
}

const Bosons &SelfEnergy::getOneLoopInteraction()
{
    if (!BB1loop.size())
    {
        generateOneLoopInteraction();
    }
    return BB1loop;
}
const Bosons &SelfEnergy::getTwoLoopInteraction()
{
    if (!BB6.size())
    {
        generateTwoLoopInteraction();
    }
    return BB6;
}
GreensProp SelfEnergy::firstOrderGreenFunction(Vec2<size_t> q, size_t type, size_t threadNum)
{
    return oneLoopGreenFunction(q, type, getOneLoopInteraction(), sw, omegaRange, omegaNum, eps, threadNum);
}
GreensProp SelfEnergy::secondOrderGreenFunction(Vec2<size_t> q, size_t type, size_t threadNum)
{
    return twoLoopGreenFunction(q, type, getOneLoopInteraction(), getTwoLoopInteraction(), sw, omegaRange, omegaNum, eps, threadNum);
}
std::vector<GreensProp> SelfEnergy::getFirstOrderGreenFunction(size_t type)
{
    if (!gp1loops[type].size())
    {
        generateOneLoopDatas(type);
    }
    return getInterpValue(gp1loops[type]);
}
std::vector<GreensProp> SelfEnergy::getSecondOrderGreenFunction(size_t type)
{
    if (!gp2loops[type].size())
    {
        generateTwoLoopDatas(type);
    }
    return getInterpValue(gp2loops[type]);
}
std::vector<std::vector<double>> SelfEnergy::getTwoMagnonDoS()
{
    if (!twoMagnonContinuumData.size())
    {
        generateTwoMagnonDoS();
    }
    return getInterpValue(twoMagnonContinuumData);
}
std::vector<std::vector<double>> SelfEnergy::getThreeMagnonDoS()
{
    if (!threeMagnonContinuumData.size())
    {
        generateThreeMagnonDoS();
    }
    return getInterpValue(threeMagnonContinuumData);
}
std::vector<std::vector<Complex>> SelfEnergy::getLinearStructureFactor(const std::vector<Vec2<double>> &sublatticePosition,
                                                                       Vec3<double> weight, int higgs)
{
    if (status < 3)
    {
        generatePathPoints(sw.getN());
    }
    std::vector<Data<std::vector<Complex>>> green0(sw.getNm());
    double delta = (omegaRange[1] - omegaRange[0]) / sw.getN();
    for (size_t t = 0; t < sw.getNm(); t++)
    {
        green0[t].reserve(nearbys.size());
        for (auto nearby : nearbys)
        {
            green0[t][nearby] = std::vector<Complex>(omegaNum + 1);
            double omega = omegaRange[0];
            for (size_t o = 0; o < omegaNum + 1; o++)
            {
                green0[t][nearby][o] += 1. / (omega - sw.getEk(t, nearby[0], nearby[1]) + eps);
                omega += delta;
            }
        }
    }
    return dynamicStructureFactor(green0, sublatticePosition, weight, higgs);
}

std::vector<std::vector<Complex>> SelfEnergy::getOneLoopStructureFactor(const std::vector<Vec2<double>> &sublatticePosition,
                                                                        Vec3<double> weight, int higgs)
{
    std::vector<Data<std::vector<Complex>>> green1(sw.getNm());
    for (size_t t = 0; t < sw.getNm(); t++)
    {
        std::cout << "preparing one loop green's function of type " << t << std::endl;
        if (!gp1loops[t].size())
        {
            generateOneLoopDatas(t);
        }
        for (auto nearby : nearbys)
        {
            green1[t][nearby] = gp1loops[t][nearby].green;
        }
    }
    return dynamicStructureFactor(green1, sublatticePosition, weight, higgs);
}