#if !defined(SPIN_H)
#define SPIN_H
#include "geometry.h"

class spin;
class Spins;

Spins operator+(const Spins &a, const Spins &b);
Spins operator*(Complex c, const Spins &s);
Spins interaction(const std::vector<size_t> &ts, const std::vector<int> &ds, const std::vector<Vec2<int>> &ps, Complex couple);
Spins interaction(size_t t1, int d1, size_t t2, int d2, Vec2<int> pos, Complex couple);
Spins Heisenberg(size_t t1, size_t t2, Vec2<int> pos, Complex couple);
Spins fieldInteraction(size_t t, Vec3<double> field);
/*****************************************************/
class spin
{

public:
    static double S;
    size_t t{};
    int direction{};
    Vec2<int> pos{};
    spin() = default;
    spin(size_t t, int d);
    spin(size_t t, int d, int x, int y);
    spin(size_t t, int d, Vec2<int> p);
    bool operator==(const spin &other) const;
    bool sameType(const spin &other) const;
    bool operator<(const spin &other) const;
};
class Spins
{
private:
public:
    double S = 0.5;
    using term = Term<spin>;
    std::vector<term> terms;
    Spins() = default;
    void add(const term &t);
    void add(const Spins &Ss, Complex coe = 1.);
    void add(spin s, Complex coe = 1.);
    size_t size() const;
    term &operator[](size_t i);
    const term &operator[](size_t i) const;
    void order();
    void shift();
    void merge();
    void chop(double eps = 1e-8);
    void simplify(double eps = 1e-8);
    size_t getNu() const;
    Spins grep(spin s) const;
};
/*****************************************************/
std::ostream &operator<<(std::ostream &os, const spin &s);
std::ostream &operator<<(std::ostream &os, const Spins &Ss);
#endif // SPIN_H