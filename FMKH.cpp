#include "src/SelfEnergy.h"
#include "src/MeanField.h"
#include <fstream>
#include <algorithm>
#include <chrono>
int main(int argc, char const *argv[])
{
    size_t threadNum = 12;
    double J = -1;
    double K = -2;
    auto HS = Heisenberg(0, 1, {0, 0}, J);
    HS = HS + Heisenberg(0, 1, {0, 1}, J);
    HS = HS + Heisenberg(0, 1, {-1, 1}, J);
    HS = HS + interaction(0, 2, 1, 2, {0, 0}, K);
    HS = HS + interaction(0, 0, 1, 0, {0, 1}, K);
    HS = HS + interaction(0, 1, 1, 1, {-1, 1}, K);
    HS.S = 0.5;
    std::vector<Vec3<double>> order = {{0, 0, 1}, {0, 0, 1}};
    size_t N = 36;
    SpinWave sw(HS, order, N);
    std::cout << N << std::endl;
    sw.BogoliubovMatrix();
    SelfEnergy se(sw, {{M_PI * 2 / 3, M_PI * 4 / 3}}, 2. / 3. * M_PI, 1.);
    Vec2<size_t> ks = {5 * N / 6, N / 6};
    Bosons Bk4 = sw.getHPTerm(4), Bk6 = sw.getHPTerm(6);
    Bk4.merge();
    Bosons BB4 = sw.BogoliubovTransformation(Bk4);
    BB4.order();
    Bk6.merge();
    Bosons BB6 = sw.BogoliubovTransformation(Bk6);
    BB6.order();
    auto start = std::chrono::steady_clock::now();
    auto [green, sigma, dispersion] = twoLoopGreenFunction(ks, 0, BB4, BB6, sw, se.getOmegaRange(), se.getOmegaNum(), se.getEps(), threadNum);
    auto end = std::chrono::steady_clock::now();
    std::cout << "calculation time "<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    return 0;
}