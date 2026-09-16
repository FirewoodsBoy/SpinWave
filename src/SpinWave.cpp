/********************************************************************
 * Implementation file for SpinWave class
 * Author: Ke Liu
 * Date: 2026.4.13
 *********************************************************************/
#include "SpinWave.h"
#include <stdexcept>

std::array<double, 9> rotation(const Vec3<double> &vec)
{
    std::array<double, 9> result;
    double l = std::sqrt(vec.dot(vec));
    if (l < 1e-8)
    {
        throw std::runtime_error("The magnetic order size is smaller than 1e-8");
    }
    for (int i = 0; i < 3; i++)
    {
        result[6 + i] = vec[i] / l;
    }
    Vec3<double> v1 = vec.cross({1., 0., 0.});
    l = std::sqrt(v1.dot(v1));
    if (l < 1e-6)
    {
        v1 = vec.cross({0., 1., 0.});
        l = std::sqrt(v1.dot(v1));
    }
    for (int i = 0; i < 3; i++)
    {
        result[3 + i] = v1[i] / l;
    }
    Vec3<double> v2 = v1.cross(vec);
    l = std::sqrt(v2.dot(v2));
    for (int i = 0; i < 3; i++)
    {
        result[i] = v2[i] / l;
    }
    return result;
}

std::unique_ptr<Complex[]> Bogoliubov(Complex *H0, size_t Nm, double *e, double eps)
{
    std::unique_ptr<Complex[]> H = copy(H0, 2 * Nm * 2 * Nm);
    std::unique_ptr<double[]> eig = std::unique_ptr<double[]>(new double[2 * Nm]);
    eigenForHermitian(H.get(), Nm * 2, eig.get(), 0);
    H = copy(H0, 2 * Nm * 2 * Nm);
    int info = cholesky(H.get(), 2 * Nm);
    if (eig[0] < 1e-7)
    {
        std::cerr << "Warning: The Hamiltonian matrix is not positive definite, adding a small gap " << eps << std::endl;
        H = copy(H0, 2 * Nm * 2 * Nm);
        for (size_t i = 0; i < 2 * Nm; i++)
        {
            H[i * 2 * Nm + i] += eps;
        }
        info = cholesky(H.get(), 2 * Nm);
        if (info)
        {
            throw std::runtime_error("Hamiltonian matrix is not positive definite after adding a small gap. Exit.");
        }
    }
    std::unique_ptr<Complex[]> J = std::make_unique<Complex[]>(2 * Nm * 2 * Nm);
    for (size_t i = 0; i < Nm; i++)
    {
        J[i * 2 * Nm + i] = 1;
        J[(i + Nm) * 2 * Nm + i + Nm] = -1;
    }
    auto J1 = matMult(H.get(), J.get(), 2 * Nm);
    auto KJK = matMult(J1.get(), H.get(), 2 * Nm, Trans::N, Trans::C);
    eigenForHermitian(KJK.get(), Nm * 2, eig.get(), true);
    std::unique_ptr<Complex[]> lambda = std::make_unique<Complex[]>(2 * Nm * 2 * Nm);
    for (size_t i = 0; i < 2 * Nm; i++)
    {
        if (eig[i] < 0)
            lambda[i * 2 * Nm + i].imag(std::sqrt(-eig[i]));
        else
            lambda[i * 2 * Nm + i] = std::sqrt(eig[i]);
    }
    triInv(H.get(), 2 * Nm);
    std::unique_ptr<Complex[]> M = std::unique_ptr<Complex[]>(new Complex[2 * Nm * 2 * Nm]);
    auto temp = matMult(H.get(), KJK.get(), 2 * Nm);
    auto T = matMult(temp.get(), lambda.get(), 2 * Nm);
    for (size_t i = 0; i < Nm; i++)
    {
        std::copy(T.get() + (2 * Nm - i - 1) * 2 * Nm, T.get() + (2 * Nm - i) * 2 * Nm, M.get() + i * 2 * Nm);
        e[i] = eig[2 * Nm - i - 1];
    }
    for (size_t i = 0; i < Nm; i++)
    {
        std::copy(T.get() + i * 2 * Nm, T.get() + (i + 1) * 2 * Nm, M.get() + (Nm + i) * 2 * Nm);
        e[Nm + i] = eig[i];
    }
    return M;
}

Complex vertex(const Bosons &Bg, const std::vector<Vec2<size_t>> &ks)
{
    Complex v{};
    for (const auto &term : Bg.terms)
    {
        Complex v0{1., 0.};
        auto k = ks.begin();
        for (const auto &bos : term.vec)
        {
            v0 = v0 * bos((*k).x, (*k).y);
            k++;
        }
        v += v0 * term.coe;
    }
    return v;
}
Complex vertex(const Bosons &Bg, const std::vector<Vec2<size_t>> &ks, Vec2<double> p)
{
    Complex v{};
    for (const auto &term : Bg.terms)
    {
        Complex v0{1., 0.};
        auto k = ks.begin();
        for (const auto &bos : term.vec)
        {
            v0 = v0 * bos((*k).x, (*k).y);
            k++;
        }
        auto pos = term.vec[0].getSubPos();
        v += v0 * term.coe * std::exp(Complex{0, 1} * (pos[0] * p[0] + pos[1] * p[1]));
    }
    return v;
}
double boseDistribution(double e, double beta)
{
    if (beta <= 0)
    {
        return 0;
    }
    return 1. / (std::exp(beta * e) - 1.);
}

// SpinWave private methods
Bosons_ SpinWave::SToB(const spin &s, const std::array<double, 9> &R) const
{
    boson_ b{s.t, 0, s.pos};
    boson_ bd(s.t, 1, s.pos);
    Bosons_ Bx, By, Bz, Bs;
    Bz.add({}, S);
    Bz.add({bd, b}, -1.);
    double sqrtS = sqrt(2. * S) / 2.;
    Complex I = Complex{0., 1.};
    Bx.add({b}, sqrtS);
    By.add({b}, -sqrtS * I);
    Bx.add({bd}, sqrtS);
    By.add({bd}, sqrtS * I);
    Bx.add({bd, b, b}, -sqrtS / 4. / S);
    By.add({bd, b, b}, sqrtS / 4. / S * I);
    Bx.add({bd, bd, b}, -sqrtS / 4. / S);
    By.add({bd, bd, b}, -sqrtS / 4. / S * I);
    double coe = sqrtS / 32. / std::pow(S, 2);
    Bx.add({bd, b, bd, b, b}, -coe);
    By.add({bd, b, bd, b, b}, coe * I);
    Bx.add({bd, bd, b, bd, b}, -coe);
    By.add({bd, bd, b, bd, b}, -coe * I);
    Bs = R[s.direction] * Bx + R[s.direction + 3] * By + R[s.direction + 6] * Bz;
    return Bs;
}

void SpinWave::setHamiltonian(const Spins &Ss)
{
    if (!order.size())
    {
        throw std::runtime_error("Magnetic order must be set first.");
    }
    Hamiltionian.clear();
    for (const auto &st : Ss.terms)
    {
        Bosons_ Bi;
        Bi.add({}, st.coe);
        for (const auto &s : st.vec)
        {
            Bi = Bi * SToB(s, rotation(order[s.t]));
        }
        Hamiltionian = Hamiltionian + Bi;
    }
    return;
}

std::unique_ptr<Complex[]> SpinWave::H0k(size_t k1, size_t k2)
{
    std::unique_ptr<Complex[]> H0 = std::make_unique<Complex[]>(2 * Nm * 2 * Nm);
    for (const auto &term : Hmf.terms)
    {
        const boson &b1 = term.vec[0];
        const boson &b2 = term.vec[1];
        H0[b1.t + Nm * (!b1.c) + (2 * Nm) * (b2.t + Nm * b2.c)] += term.coe * b1(N - k1, N - k2) * b2(k1, k2);
        H0[b2.t + Nm * (!b2.c) + (2 * Nm) * (b1.t + Nm * b1.c)] += term.coe * b2(N - k1, N - k2) * b1(k1, k2);
    }
    return H0;
}

// SpinWave public methods
SpinWave::SpinWave(const Spins &Ss, const std::vector<Vec3<double>> &order) : S(Ss.S), order(order)
{
    Nm = order.size();
    setHamiltonian(Ss);
    status = 1;
}

SpinWave::SpinWave(const Spins &Ss, const std::vector<Vec3<double>> &order, size_t N)
{
    this->N = N;
    S = Ss.S;
    this->order = order;
    Nm = order.size();
    setHamiltonian(Ss);
    status = 1;
}

// size_t SpinWave::getNm() const
// {
//     return Nm;
// }

// size_t SpinWave::getN() const
// {
//     return N;
// }

Bosons_ SpinWave::getHPTerm(size_t n) const
{
    Bosons_ Bn = Hamiltionian.getOrder(n);
    Bn.simplify();
    Bn.chop();
    return Bn;
}

void SpinWave::setN(size_t N)
{
    status = 1;
    BMs.clear();
    H0s.clear();
    Es.clear();
    Hmf.clear();
    this->N = N;
}

void SpinWave::checkGroundState()
{

    if (getHPTerm(1).size())
    {
        std::cerr << "The linear term is finite. The ground state is wrong." << std::endl;
    }
    if (status < 2)
    {
        Hmf_ = getHPTerm(2);
        Hmf = Bosons(Hmf_);
        status = 2;
    }
}

void SpinWave::setHmf(const Bosons_ &Bs)
{
    Hmf_ = Bs;
    Hmf = Bosons(Hmf_);
    status = 2;
}

void SpinWave::BogoliubovMatrix()
{
    if (status < 2)
    {
        checkGroundState();
    }
    if (BMs.empty())
    {
        Es.reserve((N + 1) * (N + 1));
        BMs.reserve(2 * Nm * 2 * Nm);
        H0s.reserve((N + 1) * (N + 1));
        for (size_t i = 0; i < (N + 1) * (N + 1); i++)
        {
            Es.emplace_back(new double[2 * Nm]);
            H0s.emplace_back(new Complex[(2 * Nm) * (2 * Nm)]);
        }
        for (size_t i = 0; i < 2 * Nm * 2 * Nm; i++)
        {
            BMs.emplace_back(new Complex[(N + 1) * (N + 1)]);
        }
    }
    for (size_t i = 0; i < N + 1; i++)
    {
        for (size_t j = 0; j < N + 1; j++)
        {
            auto H0 = H0k(i, j);
            std::copy(H0.get(), H0.get() + 2 * Nm * 2 * Nm, H0s[i * (N + 1) + j].get());
            std::unique_ptr<Complex[]> BM;
            try
            {
                BM = Bogoliubov(H0.get(), Nm, Es[i * (N + 1) + j].get(), smallGap);
            }
            catch (const std::runtime_error &e)
            {
                std::cout << i << ' ' << j << std::endl;
                matPrint(H0.get(), 2 * Nm);
            }
            for (size_t k1 = 0; k1 < Nm * 2 * Nm; k1++)
            {
                BMs[k1][i * (N + 1) + j] = BM[k1];
            }
        }
    }
    for (size_t i = 0; i < N + 1; i++)
    {
        for (size_t j = 0; j < N + 1; j++)
        {
            for (size_t k1 = 0; k1 < Nm; k1++)
            {
                for (size_t k2 = 0; k2 < Nm; k2++)
                {
                    BMs[2 * Nm * Nm + 2 * k2 * Nm + k1][i * (N + 1) + j] = std::conj(BMs[2 * k2 * Nm + k1 + Nm][(N - i) * (N + 1) + (N - j)]);
                    BMs[2 * Nm * Nm + 2 * k2 * Nm + k1 + Nm][i * (N + 1) + j] = std::conj(BMs[2 * k2 * Nm + k1][(N - i) * (N + 1) + (N - j)]);
                }
            }
        }
    }
    status = 3;
}

double SpinWave::maxEs() const
{
    if (status < 3)
    {
        throw std::runtime_error("BogoliubovMatrix() must be called before get the max Es");
    }
    auto it = std::max_element(Es.begin(), Es.end(), [](const auto &a, const auto &b)
                               { return a[0] < b[0]; });
    return (*it)[0];
}

Bosons SpinWave::BogoliubovTransformation(const Bosons &Bk) const
{
    if (status < 3)
    {
        throw std::runtime_error("BogoliubovMatrix() must be called to get Bogoliubov interaction.");
    }
    Bosons BB;
    for (auto &term : Bk.terms)
    {
        size_t n = 0;
        Bosons Bb;
        for (const auto &b0 : term.vec)
        {
            Bosons bb;
            for (size_t t = 0; t < Nm; t++)
            {
                boson b{t, false, b0.sublatticePos};
                b.setPhase(b0, BMs[t * 2 * Nm + b0.c * Nm + b0.t].get());
                if (n == 0)
                {
                    bb.add({b}, term.coe);
                }
                else
                {
                    bb.add({b});
                }
            }
            for (size_t t = 0; t < Nm; t++)
            {
                boson b{t, 1, b0.sublatticePos};
                b.setPhase(b0, BMs[(t + Nm) * 2 * Nm + b0.c * Nm + b0.t].get());
                if (n == 0)
                {
                    bb.add({b}, term.coe);
                }
                else
                {
                    bb.add({b});
                }
            }
            if (n == 0)
            {
                Bb = bb;
            }
            else
            {
                Bb = Bb * bb;
            }
            n++;
        }
        BB += Bb;
    }
    return BB;
}

std::vector<std::unique_ptr<double[]>> SpinWave::boseDist(double beta)
{
    if (status < 3)
    {
        BogoliubovMatrix();
    }
    std::vector<std::unique_ptr<double[]>> bose;
    bose.reserve(Nm);
    for (size_t t = 0; t < Nm; t++)
    {
        bose.emplace_back(new double[(N + 1) * (N + 1)]);
        for (size_t i = 0; i < (N + 1) * (N + 1); i++)
        {
            bose[t][i] = boseDistribution(Es[i][t], beta);
        }
    }
    return bose;
}

Complex SpinWave::getBMs(int t, int i) const
{
    return BMs[t][i];
}

double SpinWave::getEk(size_t t, size_t i, size_t j) const
{
    if (status < 3)
    {
        throw std::runtime_error("Function BogoliubovMatrix() must be called before get Ek.");
    }
    return Es[i * (N + 1) + j][t];
}