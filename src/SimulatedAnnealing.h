#include "lapack.h"
#include "geometry.h"
#include "spin.h"
#include <random>
#include <unordered_map>

class SimulatedAnnealing;
Vec3<double> randomVec3(std::mt19937_64 &re);

class SimulatedAnnealing
{
    using Config = std::unordered_map<Vec3<size_t>, Vec3<double>>;

public:
    Spins Ss;
    size_t Nu = 1;
    size_t seed{};
    Config fixedConfig;
    std::unordered_map<Vec2<size_t>, Spins> Sfs;
    std::mt19937_64 re;

    Vec3<double> effectiveField(Vec2<size_t> latticeSize, const Config &conf, Vec2<size_t> index, size_t t);
    double totalEnergy(Vec2<size_t> latticeSize, const Config &conf) const;
    Vec3<double> heatBathSampling(const Vec3<double> &ef, std::mt19937_64 &re, double beta) const;
    void getNewConfig(Vec2<size_t> latticeSize, Config &conf, std::mt19937_64 &re, double beta);

public:
    size_t maxL = 3;
    size_t Nmc = 50;
    double betaMin = 1;
    double betaMax = 1e16;
    double alpha = 1.2;
    bool ifPrint = true;
    SimulatedAnnealing(Spins &Ss);
    SimulatedAnnealing(Spins &Ss, size_t seed);
    Config groundState(Vec2<size_t> latticeSize);
};
