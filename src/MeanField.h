/********************************************************************
 * This is a head file defines a class named MeanFileld, which can
 * get the mean-field value of boson pair Vec2<boson_>, decouping the
 * interactions, and set the Hmf to SpinWave to get the new wave-funcition.
 * The self-consistent calculation can also be processed in this class.
 * Author: Ke Liu
 * Date: 2026.4.18
 *********************************************************************/
#if !defined(SELFCONSISTENT_H)
#define SELFCONSISTENT_H

#include "SpinWave.h"
#include <unordered_map>
#include <unordered_set>
class MeanField;

class MeanField
{
private:
    SpinWave &sw;
    size_t status = 0;
    std::vector<std::unique_ptr<double[]>> bose;
    std::unordered_map<Vec2<boson_>, Complex> pairValues;
    std::unordered_map<Vec2<boson_>, Pairs_> meanFields;
    Bosons_ BImf, BInt;
    Bosons_ H0;
    double beta = 0;
    void boseDistribution();
    void getPossibleMeanFields(const Bosons_ &Bs);
    void generateBImf();

public:
    MeanField(SpinWave &sw);
    void setBeta(double beta);
    Complex pairValue(const Vec2<boson_> &bp);
    Complex operaorValue(const Pairs_ &ps);
    void setInteraction(const Bosons_ &Bs);
    const Bosons_ &getBImf() const;
    // This function push the decouping mean field to linear spin wave.
    void setMeanFiledHamiltonian();
    void selfConsistent(bool print = 1, double eps = 1e-6, size_t n_max = 100);
    std::unordered_map<Vec2<boson_>, Complex> getPairValues() const
    {
        return pairValues;
    }
};

#endif // SELFCONSISTENT_H