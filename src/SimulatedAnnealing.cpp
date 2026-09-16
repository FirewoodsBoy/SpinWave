#include "SimulatedAnnealing.h"
#include "SpinWave.h"
#include <time.h>
std::uniform_real_distribution<> PiDis(0, M_PI);
std::uniform_real_distribution<> OneDis(0, 1);
SimulatedAnnealing::SimulatedAnnealing(Spins &Ss) : Ss(Ss)
{
    Nu = Ss.getNu();
    seed = static_cast<size_t>(time(0));
    re.seed(seed);
}
SimulatedAnnealing::SimulatedAnnealing(Spins &Ss, size_t seed) : Ss(Ss), seed(seed)
{
    Nu = Ss.getNu();
    re.seed(seed);
}
Vec3<double> SimulatedAnnealing::effectiveField(Vec2<size_t> latticeSize, const Config &conf, Vec2<size_t> index, size_t t)
{
    Vec3<double> ef;
    for (size_t d = 0; d < 3; d++)
    {
        if (!Sfs.count({t, d}))
        {
            Sfs[{t, d}] = Ss.grep({t, static_cast<int>(d)});
        }
        Spins Sf = Sfs[{t, d}];
        Complex fd{};
        for (size_t i = 0; i < Sf.size(); i++)
        {
            Complex f = Sf[i].coe;
            for (auto s : Sf[i].vec)
            {
                size_t x = (index[0] + 100 * latticeSize[0] + s.pos[0]) % latticeSize[0];
                size_t y = (index[1] + 100 * latticeSize[1] + s.pos[1]) % latticeSize[1];
                f *= conf.at(Vec3{s.t, x, y})[s.direction] * Ss.S;
            }
            fd += f;
        }
        ef[d] = fd.real();
    }
    return ef;
}
double SimulatedAnnealing::totalEnergy(Vec2<size_t> latticeSize, const Config &conf) const
{
    Complex energy{};
    for (size_t x0 = 0; x0 < latticeSize[0]; x0++)
    {
        for (size_t y0 = 0; y0 < latticeSize[1]; y0++)
        {
            for (size_t i = 0; i < Ss.size(); i++)
            {
                Complex e = Ss[i].coe;
                for (auto s : Ss[i].vec)
                {
                    size_t x = (x0 + 100 * latticeSize[0] + s.pos[0]) % latticeSize[0];
                    size_t y = (y0 + 100 * latticeSize[1] + s.pos[1]) % latticeSize[1];
                    e *= conf.at({s.t, x, y})[s.direction] * Ss.S;
                }
                energy += e;
            }
        }
    }
    return energy.real();
}
Vec3<double> SimulatedAnnealing::heatBathSampling(const Vec3<double> &ef, std::mt19937_64 &re, double beta) const
{
    double H = sqrt(pow(ef[0], 2) + pow(ef[1], 2) + pow(ef[2], 2));
    double r = OneDis(re);
    double S = Ss.S;
    double cosTheta = log(1 - r + r * exp(-2 * beta * H * S)) / beta / H / S + 1;
    double theta = acos(cosTheta);
    double sinTheta = sin(theta);
    double phi = 2 * PiDis(re);
    std::array<double, 3> ra0 = {sinTheta * cos(phi), sinTheta * sin(phi), cosTheta};
    auto R = rotation({ef[0], ef[1], ef[2]});
    std::array<double, 3> ra;
    realMv(R.data(), ra0.data(), ra.data(), 3);
    Vec3<double> spinAxis(-ra[0], -ra[1], -ra[2]);
    return spinAxis.getNormalizedVec3();
}
void SimulatedAnnealing::getNewConfig(Vec2<size_t> latticeSize,
                                      Config &conf, std::mt19937_64 &re, double beta)
{
    for (size_t count = 0; count < Nmc; count++)
    {
        for (size_t i = 0; i < latticeSize[0]; i++)
        {
            for (size_t j = 0; j < latticeSize[1]; j++)
            {
                for (size_t t = 0; t < Nu; t++)
                {
                    if (!fixedConfig.count({t, i, j}))
                    {
                        Vec3<double> ef = effectiveField(latticeSize, conf, {i, j}, t);
                        conf[{t, i, j}] = heatBathSampling(ef, re, beta);
                    }
                }
            }
        }
    }
}

SimulatedAnnealing::Config SimulatedAnnealing::groundState(Vec2<size_t> latticeSize)
{
    double beta = betaMin / alpha;
    Config config;
    config.reserve(latticeSize[0] * latticeSize[1]);
    for (size_t i = 0; i < latticeSize[0]; i++)
    {
        for (size_t j = 0; j < latticeSize[1]; j++)
        {
            for (size_t t = 0; t < Nu; t++)
            {
                if (fixedConfig.count({t, i, j}))
                    config[{t, i, j}] = fixedConfig[{t, i, j}].getNormalizedVec3();
                else
                    config[{t, i, j}] = randomVec3(re);
            }
        }
    }
    do
    {
        beta = beta * alpha;
        getNewConfig(latticeSize, config, re, beta);
        if (ifPrint)
        {
            std::cout << "config energy in beta=" << beta << " is " << totalEnergy(latticeSize, config) / (latticeSize[0] * latticeSize[1]) << std::endl;
        }

    } while (beta < betaMax);
    return config;
}
/*****************************************************************************/
Vec3<double> randomVec3(std::mt19937_64 &re)
{
    return Vec3<double>{OneDis(re) - 0.5, OneDis(re) - 0.5, OneDis(re) - 0.5}.getNormalizedVec3();
}