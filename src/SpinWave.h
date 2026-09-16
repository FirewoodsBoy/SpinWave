/********************************************************************
 * This head file defines a class named SpinWave, which could do the
 * Holstein-Primekoff transformation, Bogoliubov transformation and
 * calculate the vertex of class Bosons. When the Hmf is not positive
 * definite, it will add a small gap and try to solve it. If the Hmf
 * is also not positive definite after adding the small gap, a run-
 * time error will throw.
 * Author: Ke Liu
 * Date: 2026.4.13
 *********************************************************************/
#if !defined(SPINWAVE_H)
#define SPINWAVE_H

#include "geometry.h"
#include <math.h>
#include <array>
#include <iostream>
#include <mutex>
#include "spin.h"
#include "Boson.h"
#include "lapack.h"
class SpinWave;

std::array<double, 9> rotation(const Vec3<double> &vec);
std::unique_ptr<Complex[]> Bogoliubov(Complex *H0, size_t Nm, double *e, double eps);
Complex vertex(const Bosons &Bg, const std::vector<Vec2<size_t>> &ks);
Complex vertex(const Bosons &Bg, const std::vector<Vec2<size_t>> &ks, Vec2<double> p);
double boseDistribution(double e, double beta);
// std::shared_ptr<Complex[]> sharedMult(const Complex *a, const Complex *b, size_t l);
/***********************************************************/
class SpinWave
{
    friend class MeanField;

private:
    double S = 0.5;
    Bosons_ Hamiltionian, Hmf_;
    Bosons Hmf;
    std::vector<Vec3<double>> order;
    std::vector<std::unique_ptr<Complex[]>> BMs, H0s;
    std::vector<std::unique_ptr<double[]>> Es;
    size_t Nm = 0;
    size_t &N = boson::N;

    bool ifInitial = false;
    size_t status = 0;

    void setHamiltonian(const Spins &Ss);

public:
    std::unique_ptr<Complex[]> H0k(size_t k1, size_t k2);
    // void setOrder(const std::vector<Vec3<double>> &order)
    // {
    //     this->order = order;
    //     Nm = order.size();
    // };
    double smallGap = 1e-5;
    SpinWave(const Spins &Ss, const std::vector<Vec3<double>> &order);
    SpinWave(const Spins &Ss, const std::vector<Vec3<double>> &order, size_t N);
    Bosons_ SToB(const spin &s, const std::array<double, 9> &R) const;
    size_t getNm() const
    {
        return Nm;
    }
    size_t getN() const
    {
        return N;
    }
    size_t getStatus() const
    {
        return status;
    }
    Bosons_ getHPTerm(size_t n) const;
    void setN(size_t N);
    void checkGroundState();
    void setHmf(const Bosons_ &Bs);
    void BogoliubovMatrix();
    double maxEs() const;
    Bosons BogoliubovTransformation(const Bosons &Bk) const;
    std::vector<std::unique_ptr<double[]>> boseDist(double beta);
    Complex getBMs(int t, int i) const;
    double getEk(size_t t, size_t i, size_t j) const;
    Bosons_ getHmf() const
    {
        return Hmf_;
    }
    const std::vector<Vec3<double>> &getMagneticOrder() const
    {
        return order;
    }
};

#endif // SPINWAVE_H