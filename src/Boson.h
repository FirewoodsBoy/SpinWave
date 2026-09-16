/********************************************************************
 * This is a head file for boson operators and their polynormal in
 * both real and momentum space. The operator polynormals can do
 * normal order and simplify such as merge and chop. These operators
 * are used for vertex calculation and "Hartree-Fock" decouping.
 * Author: Ke Liu
 * Date: 2026.4.8
 *********************************************************************/
#if !defined(BOSON_H)
#define BOSON_H

#define _USE_MATH_DEFINES

#include "geometry.h"
#include <complex>
#include <memory>
#include <vector>
#include <ostream>
#include <algorithm>

typedef std::complex<double> Complex;
namespace SpinWaveParameter
{
    const double ep = 1e-6;
} // namespace SpinWaveParameter
/***************Class defined**************/
// The "_" means that these are defined in real space.
class boson_; // The basic operator in real space.

class boson;   // The basic operator in k space.
class Bosons_; // The polynomial of class boson_.

class Bosons; // The polynomial of class boson.
class Pairs_; // Wick pairs in real space.
/********************************************** */
template <typename BOSON>
Term<BOSON> operator*(Complex c, const Term<BOSON> &term);
Bosons_ operator*(Complex c, const Bosons_ &Bs);
Bosons_ operator*(const Bosons_ &Bs1, const Bosons_ &Bs2);
Bosons operator*(const Bosons &Bk1, const Bosons &Bk2);
bool operator==(const Term<boson_> &a, const Term<boson_> &b);
// template <size_t N>
// bool operator==(const Term<boson> &a,const Term<boson> &b);
/***************ostream function*************/
std::ostream &operator<<(std::ostream &os, const boson_ &a);

std::ostream &operator<<(std::ostream &os, const boson &a);
template <typename BOSON>
std::ostream &operator<<(std::ostream &os, const std::vector<BOSON> &bos);
template <typename BOSON>
std::ostream &operator<<(std::ostream &os, const Term<BOSON> &term);
std::ostream &operator<<(std::ostream &os, const Bosons_ &Bs);

std::ostream &operator<<(std::ostream &os, const Bosons &Bs);
template <typename BOSON>
std::ostream &operator<<(std::ostream &os, const Vec2<BOSON> pair);
/********************************************************************************/

class boson_
{
    /******************
     * A class of bosonic operator in real space.
     * Contains variable of the position of the operator.
     * member variables: type t (int), creation or annihilation c (bool) and position p (Vec2<int>).
     * ***************/
public:
    size_t t{};
    bool c{};
    Vec2<int> p{};
    boson_() = default;
    boson_(size_t t, bool c);
    boson_(size_t t, bool c, int x, int y);
    boson_(size_t t, bool c, Vec2<int> p);
    void setPosition(int x, int y);
    bool operator==(const boson_ &f) const;
    bool sameType(const boson_ &f) const;
};
Vec2<boson_> make_pair(const boson_ &a, const boson_ &b);
namespace std
{
    template <>
    struct hash<boson_>
    {
        size_t operator()(const boson_ &b) const noexcept
        {
            size_t h1 = hash<size_t>{}(b.t);
            size_t h2 = hash<bool>{}(b.c);
            size_t h3 = hash<Vec2<int>>{}(b.p);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}
class Bosons_
{
public:
    using term_ = Term<boson_>;
    std::vector<term_> terms;
    Bosons_() = default;
    void clear();
    // add terms to polynomials
    void add(const std::vector<boson_> &bs, Complex co = 1.);
    void addPair(const Vec2<boson_> &pair, Complex co = 1.);
    Bosons_ operator+(const term_ &term) const;
    Bosons_ operator+(const Bosons_ &Bs) const;
    void operator+=(const Bosons_ &Bs);
    // Get the term number.
    auto size() const
    {
        return terms.size();
    }
    // Get the order (operator number) n in term i.
    auto size(size_t i) const
    {
        return terms[i].vec.size();
    }
    // Return term i.
    term_ &operator[](size_t i) { return terms[i]; }
    const term_ &operator[](size_t i) const { return terms[i]; }
    // Return terms with order n.
    Bosons_ getOrder(size_t n) const;
    // Grep all the specific terms in the object
    // e.g. t={0,1}, c={1,0} ==> grep all the terms of a^b
    Bosons_ grep(std::vector<size_t> t, std::vector<bool> c);
    // Normal order the operator
    void order();
    // Only relative coordinates matters.
    // This function set all the positions of bos[i][0] to be zero.
    void shift();
    void merge();
    void chop(double eps = SpinWaveParameter::ep);
    void simplify(double eps = SpinWaveParameter::ep);
    void symmetrize(double eps = SpinWaveParameter::ep);
    Bosons toK(Vec2<double> subPos = {0, 0});
};

class Pairs_
{
private:
    using Pair_ = Vec2<boson_>;
    using PairList = std::vector<Pair_>;
    bool WickPair(const std::vector<boson_> &list, const PairList &pairList,
                  std::vector<PairList> &pairLists);

public:
    std::vector<Term<PairList>> pairListsTerms;
    Pairs_() = default;
    Pairs_(const Bosons_ &Bs);
    Pairs_(const Bosons_ &Bs, const Pair_ &pair);
    Pairs_(Bosons_ Bs, boson_ b);
};

class boson
{

private:
    std::shared_ptr<Complex[]> phase;
    bool trivial = 0;
    // Complex &operator[](size_t i)
    // {
    //     // if (!phase.get())
    //     // {
    //     //     throw std::runtime_error("The boson phase does not initial.");
    //     // }
    //     return phase[i];
    // }
    // Complex &operator()(size_t i, size_t j)
    // {
    //     return this->operator[](i * (N + 1) + j);
    // }

    static size_t N;
    // boson(const boson& other)=default;
    // boson & operator=(const boson& other)=default;
    friend Complex vertex(const Bosons &Bg, const std::vector<Vec2<size_t>> &ks);

public:
    size_t t{};
    bool c{};
    Vec2<double> sublatticePos{};
    friend class SpinWave;
    friend class Bosons;
    static std::vector<double> linspace(double s, double e);
    boson() = default;
    boson(size_t t, bool c);
    boson(size_t t, bool c, Vec2<double> subPos);
    boson(const boson_ &br);
    boson(const boson_ &br, Vec2<double> subPos);
    void copy(boson &other) const;
    const Complex &operator[](size_t i) const
    {
        return phase[i];
    }
    const Complex &operator()(size_t i, size_t j) const
    {
        return this->operator[](i * (N + 1) + j);
    }
    // used for Bogoliubov transformation
    void setPhase(const boson &other, Complex *p);
    long use_count() const;
    bool ifTrivial() const;
    bool sameType(const boson &other) const;
    Complex contract(const boson &other) const;
    Vec2<double> getSubPos() const
    {
        return sublatticePos;
    }
};

class Bosons
{
public:
    using term = Term<boson>;
    std::vector<term> terms;
    Bosons() = default;

    void clear();
    void add(std::vector<boson> bs, Complex co = 1.);
    Bosons operator+(const term &ter) const;
    Bosons operator+(const Bosons &Bs) const;
    void operator+=(const Bosons &Bs);
    void operator*=(Complex c);
    auto size() const
    {
        return terms.size();
    }
    term &operator[](size_t i) { return terms[i]; }
    const term &operator[](size_t i) const { return terms[i]; }
    Bosons getOrder(size_t n) const;
    Bosons grep(const std::vector<boson> &bs) const;
    void order();
    void merge();
    std::vector<term> getTerms(const std::vector<boson> &bs);
    Complex constant() const;
    Bosons(const Bosons_ &Bs);
    Bosons(const Bosons_ &Bs, const std::vector<Vec2<double>> &subPositions);
    Bosons swap(std::vector<size_t> indexs) const;
    Bosons swap(size_t i, size_t j) const;
    std::unique_ptr<Complex[]> reduce(std::vector<Vec2<size_t>> ks, size_t i0, size_t j0) const;
};

/********************************************/
template <typename BOSON>
Term<BOSON> operator*(Complex c, const Term<BOSON> &term)
{
    Term<BOSON> t = term;
    t.coe *= c;
    return t;
}
/***************ostream function*************/
template <typename BOSON>
std::ostream &operator<<(std::ostream &os, const std::vector<BOSON> &bos)
{
    for (auto const &b : bos)
    {
        os << b;
    }
    return os;
}
template <typename BOSON>
std::ostream &operator<<(std::ostream &os, const Term<BOSON> &t)
{
    os << t.coe << t.vec;
    return os;
}
template <typename BOSON>
std::ostream &operator<<(std::ostream &os, const Vec2<BOSON> pair)
{
    os << '[' << pair[0] << ' ' << pair[1] << ']';
    return os;
}
#endif // BOSON_H