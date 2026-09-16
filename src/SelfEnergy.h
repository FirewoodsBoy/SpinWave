#if !defined(SELFENERGY_H)
#define SELFENERGY_H

#include "SpinWave.h"
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <algorithm>
#include <tuple>

struct GreensProp;
class SelfEnergy;

void shiftMomentum(double &k); // shift momentum to (-PI,PI]
template <typename T>
T bilinearInterp(const std::array<T, 4> &pvs, Vec2<double> p);
template <typename T>
std::vector<T> bilinearInterp(const std::array<std::vector<T> *, 4> &pvs, Vec2<double> p);
GreensProp bilinearInterp(const std::array<GreensProp *, 4> &gps, Vec2<double> p);
// Calculate the one loop self-energy contributed only by Cubic term
// return tuple(vector<Complex> sigma, vector<Complex> dos)
// B3 is the cubic Bogoliubov interaction terms (quadratic and quartic term can be included but have no effect.)
std::vector<Complex> oneLoop(Vec2<size_t> q, size_t type, const Bosons &B3, const SpinWave &sw,
                             Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum);
std::vector<Complex> twoLoop(Vec2<size_t> q, size_t type, const Bosons &B4, const SpinWave &sw,
                             Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum);
std::vector<double> twoMagnonDos(Vec2<size_t> q, const SpinWave &sw,
                                 Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum);
std::vector<double> threeMagnonDos(Vec2<size_t> q, const SpinWave &sw,
                                   Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum);
// Calculate the one loop green function to 1st order correction.
// BB is the Bogoliubov interaction terms up to quartic term.
GreensProp oneLoopGreenFunction(Vec2<size_t> q, size_t type, const Bosons &BB, const SpinWave &sw,
                                Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum);
GreensProp twoLoopGreenFunction(Vec2<size_t> q, size_t type, const Bosons &B4, const Bosons &B62, const SpinWave &sw,
                                Vec2<double> omegaRange, size_t omegaNum, Complex eps, size_t threadNum);
struct GreensProp
{
    std::vector<Complex> green;
    std::vector<Complex> sigma;
    // std::vector<Complex> dos;
    double dispersion{};
    void reserve(size_t n)
    {
        green.reserve(n);
        sigma.reserve(n);
        // dos.reserve(n);
    }
};
// class SelfEnergy is to calculate self-energy, spectral funcion, density of state and dispersion.
// status: 1. The initial value;
//         2. Check sw have been diagonlized, and eps is set (default 1.5*EMax/N);
//         3. pathPoints are genertated
class SelfEnergy
{

private:
    const SpinWave &sw;
    // The near by points with integer coordinates of true points.
    std::unordered_map<Vec2<double>, std::array<Vec2<size_t>, 4>> pointsNearby;
    std::unordered_map<Vec2<double>, std::array<Vec2<double>, 4>> pointsNearbyTrueCoor;
    // The true coordinate points position in the parallelogram (for interp)
    std::unordered_map<Vec2<double>, Vec2<double>> pointsPosition;
    // The true coordinate points in the path.
    std::vector<Vec2<double>> pathPoints;
    // The index of the endpoints.
    std::vector<size_t> pointsIndex;
    // All the points of integer coordinates need to calculate.
    std::unordered_set<Vec2<size_t>> nearbys;
    std::unordered_map<Vec2<double>, Vec2<size_t>> nearbysTrueCoor;
    template <typename T>
    using Data = std::unordered_map<Vec2<size_t>, T>;
    std::vector<Vec2<double>> path = {{0, 0}, {M_PI, M_PI}};
    double theta = M_PI_2;
    double rate = 1.;
    size_t status = 1;
    size_t threadNum = 1;
    size_t omegaNum{};
    Complex eps{};
    Vec2<double> omegaRange{};
    Bosons BB1loop, BB6;
    // std::vector<Vec2<double>> sublatticePosition;

    std::vector<Data<GreensProp>> gp1loops;
    std::vector<Data<GreensProp>> gp2loops;
    Data<std::vector<double>> twoMagnonContinuumData;
    Data<std::vector<double>> threeMagnonContinuumData;

    double length(Vec2<double> a1, Vec2<double> a2);
    std::array<Vec2<size_t>, 4> nearbyPoints(Vec2<double> k, Vec2<double> &p);
    void generateOneLoopInteraction();
    void generateTwoLoopInteraction();
    void generateOneLoopDatas(size_t type);
    void generateTwoLoopDatas(size_t type);
    void generateTwoMagnonDoS();
    void generateThreeMagnonDoS();
    std::vector<std::vector<Complex>> dynamicStructureFactor(std::vector<Data<std::vector<Complex>>> &green,
                                                             const std::vector<Vec2<double>> &sublatticePosition, Vec3<double> weight, int higgs) const;
    template <typename T>
    std::vector<std::vector<T>> getInterpValue(Data<std::vector<T>> &data) const
    {
        std::vector<std::vector<T>> interpValue;
        interpValue.reserve(pathPoints.size());
        for (auto p : pathPoints)
        {
            const auto &nearby = pointsNearby.at(p);
            std::array<std::vector<T> *, 4> temp = std::array<std::vector<T> *, 4>{&data[nearby[0]], &data[nearby[1]], &data[nearby[2]], &data[nearby[3]]};
            interpValue.push_back(bilinearInterp(temp, pointsPosition.at(p)));
        }
        return interpValue;
    }
    template <typename T>
    std::vector<T> getInterpValue(Data<T> &data) const
    {
        std::vector<T> interpValue;
        interpValue.reserve(pathPoints.size());
        for (auto p : pathPoints)
        {
            const auto &nearby = pointsNearby.at(p);
            std::array<T, 4> temp = {data[nearby[0]], data[nearby[1]], data[nearby[2]], data[nearby[3]]};
            interpValue.push_back(bilinearInterp(temp, pointsPosition.at(p)));
        }
        return interpValue;
    }
    std::vector<GreensProp> getInterpValue(Data<GreensProp> &data) const;

public:
    bool print = true;
    SelfEnergy(const SpinWave &sw);
    SelfEnergy(const SpinWave &sw, const std::vector<Vec2<double>> &path);
    SelfEnergy(const SpinWave &sw, const std::vector<Vec2<double>> &path, double theta, double rate);
    void generatePathPoints(size_t n);
    std::vector<double> getHmfDispersion(size_t type);
    void checkSpinWave();
    GreensProp firstOrderGreenFunction(Vec2<size_t> q, size_t type, size_t threadNum);
    GreensProp secondOrderGreenFunction(Vec2<size_t> q, size_t type, size_t threadNum);
    std::vector<GreensProp> getFirstOrderGreenFunction(size_t type);
    std::vector<GreensProp> getSecondOrderGreenFunction(size_t type);
    std::vector<std::vector<double>> getTwoMagnonDoS();
    std::vector<std::vector<double>> getThreeMagnonDoS();
    std::vector<std::vector<Complex>> getLinearStructureFactor(const std::vector<Vec2<double>> &sublatticePosition,
                                                               Vec3<double> weight = {1, 1, 1}, int higgs = 0);
    std::vector<std::vector<Complex>> getOneLoopStructureFactor(const std::vector<Vec2<double>> &sublatticePosition,
                                                               Vec3<double> weight = {1, 1, 1}, int higgs = 1);
    const Bosons &getOneLoopInteraction();
    const Bosons &getTwoLoopInteraction();
    // Inline method
    std::vector<size_t> getPointIndex() const
    {
        return pointsIndex;
    }
    void setPath(const std::vector<Vec2<double>> &path)
    {
        this->path = path;
        if (status == 3)
        {
            status = 2;
        }
    }
    void setThetaAndRate(double theta, double rate)
    {
        this->theta = theta;
        this->rate = rate;
        if (status == 3)
        {
            status = 2;
        }
    }
    void setEps(double eps)
    {
        this->eps.imag(eps);
    }
    Complex getEps() const
    {
        return eps;
    }
    void setOmegaRange(double omegaMax)
    {
        omegaRange = {0, omegaMax};
    }
    void setOmegaRange(Vec2<double> omegaRange)
    {
        this->omegaRange = omegaRange;
    }
    Vec2<double> getOmegaRange() const
    {
        return omegaRange;
    }
    void setOmegaNum(size_t omegaNum)
    {
        this->omegaNum = omegaNum;
    }
    size_t getOmegaNum() const
    {
        return omegaNum;
    }
    size_t getStatus() const
    {
        return status;
    }
    void setThreadNum(size_t threadNum)
    {
        this->threadNum = threadNum;
    }
    size_t getThreadNum()
    {
        return threadNum;
    }
    double getDeltaOmega()
    {
        return (omegaRange[1] - omegaRange[0]) / (omegaNum);
    }
    template <typename T>
    T onShellValue(std::vector<T> values, double omega)
    {
        if (omega < omegaRange[0] || omega > omegaRange[1])
        {
            throw std::runtime_error("The omega range is small as dipersion is out of range.");
        }
        double delta = getDeltaOmega();
        double pos = (omega - omegaRange[0]) / delta;
        double alpha = pos - std::floor(pos);
        // std::cout<<alpha <<std::endl;
        size_t floor = static_cast<size_t>(std::floor(pos));
        size_t ceil = floor + 1;
        if (ceil > omegaNum + 1)
        {
            ceil = floor;
        }
        return values[floor] + alpha * (values[ceil] - values[floor]);
    }
};

/******************************************************************************/
template <typename T>
T bilinearInterp(const std::array<T, 4> &pvs, Vec2<double> p)
{
    Vec2<double> dp{1 - p[0], 1 - p[1]};
    return dp[0] * dp[1] * pvs[0] + p[0] * dp[1] * pvs[1] + dp[0] * p[1] * pvs[2] + p[0] * p[1] * pvs[3];
}
template <typename T>
std::vector<T> bilinearInterp(const std::array<std::vector<T> *, 4> &pvs, Vec2<double> p)
{
    size_t size = pvs[0][0].size();
    std::vector<T> result;
    result.reserve(size);
    auto ips0 = pvs[0][0].cbegin();
    auto ips1 = pvs[1][0].cbegin();
    auto ips2 = pvs[2][0].cbegin();
    auto ips3 = pvs[3][0].cbegin();
    // for (size_t i = 0; i < 4; i++)
    // {
    //     ips[i] = pvs[i][0].cbegin();
    // }
    for (size_t i = 0; i < size; i++)
    {
        result.push_back(bilinearInterp(std::array<T, 4>{*ips0, *ips1, *ips2, *ips3}, p));
        ips0++;
        ips1++;
        ips2++;
        ips3++;
    }
    return result;
}

#endif // SELFENERGY_H