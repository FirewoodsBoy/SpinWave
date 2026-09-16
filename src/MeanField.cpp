/********************************************************************
 * Implementation file for MeanField class
 * Author: Ke Liu
 * Date: 2026.4.18
 *********************************************************************/
#include "MeanField.h"
#include <iostream>

// MeanField private methods
void MeanField::boseDistribution()
{
    bose = sw.boseDist(beta);
    pairValues.clear();
    status = 2;
}

void MeanField::getPossibleMeanFields(const Bosons_ &Bs)
{
    std::unordered_set<Vec2<boson_>> possiblePairs;
    for (const auto &term : Bs.terms)
    {
        size_t l = term.vec.size();
        if (l > 0 && l % 2 == 0)
        {
            for (size_t i = 0; i < l - 1; i++)
            {
                for (size_t j = i + 1; j < l; j++)
                {
                    possiblePairs.insert(make_pair(term.vec[i], term.vec[j]));
                }
            }
        }
    }
    for (const auto &pair : possiblePairs)
    {
        meanFields[pair] = Pairs_{Bs, pair};
    }
}

void MeanField::generateBImf()
{
    if (status < 2)
    {
        boseDistribution();
    }
    if (meanFields.empty())
    {
        getPossibleMeanFields(BInt);
    }
    BImf.clear();
    for (const auto &pv : meanFields)
    {
        BImf.addPair(pv.first, operaorValue(pv.second));
    }
    BImf.simplify();
    status = 3;
}

// MeanField public methods
MeanField::MeanField(SpinWave &sw) : sw(sw)
{
    BInt = sw.getHPTerm(4);
    H0 = sw.getHPTerm(2);
    status = 1;
}

void MeanField::setBeta(double beta)
{
    this->beta = beta;
    boseDistribution();
    pairValues.clear();
}

Complex MeanField::pairValue(const Vec2<boson_> &bp)
{
    if (status < 2)
    {
        boseDistribution();
    }
    size_t N = sw.getN();
    int Nm = sw.getNm();
    const boson_ &b1 = bp[0];
    const boson_ &b2 = bp[1];
    static auto ks = boson::linspace(-M_PI, M_PI);
    Complex v{};
    for (int t = 0; t < Nm; t++)
    {
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                v += sw.getBMs(b1.t + Nm * b1.c + t * (2 * Nm), i * (N + 1) + j) * sw.getBMs(b2.t + Nm * b2.c + (t + Nm) * (2 * Nm), (N - i) * (N + 1) + (N - j)) * std::exp(Complex{0., 1.} * ((b1.p[0] - b2.p[0]) * ks[i] + (b1.p[1] - b2.p[1]) * ks[j])) * (1 + bose[t][i * (N + 1) + j]);
                v += sw.getBMs(b1.t + Nm * b1.c + (t + Nm) * (2 * Nm), i * (N + 1) + j) * sw.getBMs(b2.t + Nm * b2.c + t * (2 * Nm), (N - i) * (N + 1) + (N - j)) * std::exp(Complex{0., 1.} * ((b1.p[0] - b2.p[0]) * ks[i] + (b1.p[1] - b2.p[1]) * ks[j])) * bose[t][(N - i) * (N + 1) + N - j];
            }
        }
    }
    return v / std::pow(N, 2);
}

Complex MeanField::operaorValue(const Pairs_ &ps)
{
    Complex v{};
    for (const auto &term : ps.pairListsTerms)
    {
        Complex tvs;
        for (const auto &pairs : term.vec)
        {
            Complex tv{1., 0};
            for (const auto &pair : pairs)
            {
                if (!pairValues.count(pair))
                {
                    pairValues[pair] = pairValue(pair);
                }
                tv *= pairValues[pair];
            }
            tvs += tv;
        }
        v += term.coe * tvs;
    }
    return v;
}

void MeanField::setInteraction(const Bosons_ &Bs)
{
    BInt = Bs;
    meanFields.clear();
    BImf.clear();
    status = 2;
}

const Bosons_ &MeanField::getBImf() const
{
    return BImf;
}

void MeanField::setMeanFiledHamiltonian()
{
    if (status < 3)
    {
        generateBImf();
    }
    Bosons_ Hmf = BImf + H0;
    Hmf.simplify();
    Hmf.symmetrize();
    sw.setHmf(Hmf);
    status = 1;
}

void MeanField::selfConsistent(bool print, double eps, size_t n_max)
{
    setMeanFiledHamiltonian();
    // std::cout << BImf << std::endl;
    auto BImf0 = BImf;
    sw.BogoliubovMatrix();
    for (size_t i = 0; i < n_max; i++)
    {
        setMeanFiledHamiltonian();
        // std::cout << BImf << std::endl;
        // double diffMax = 0;
        // Vec2<boson_> bp;
        // for (auto i : pv0)
        // {
        //     Complex diff = i.second - pairValues[i.first];
        //     if (diffMax < std::abs(diff.real()) + std::abs(diff.imag()))
        //     {
        //         diffMax = std::abs(diff.real()) + std::abs(diff.imag());
        //         bp = i.first;
        //     }
        // }
        // std::cout << bp << ' ' << pv0[bp] << ' ' << pairValues[bp] << std::endl;
        BImf0 += (-1.) * BImf;
        BImf0.simplify(eps);
        if (print)
        {
            std::cout << "Mean field difference: " << BImf0 << std::endl;
        }

        // std::cout << "Max mean field difference: " << diffMax << std::endl;
        // pv0 = pairValues;

        sw.BogoliubovMatrix();
        if (!BImf0.size())
        {
            return;
        }
        BImf0 = BImf;
        // if (diffMax < eps)
        // {
        //     return;
        // }
    }
}