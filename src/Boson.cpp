/********************************************************************
 * Implementation file for boson operators
 * Author: Ke Liu
 * Date: 2026.4.8
 *********************************************************************/
#include "Boson.h"

size_t boson::N = 120;

// boson_ implementations
boson_::boson_(size_t t, bool c) : t(t), c(c) {}
boson_::boson_(size_t t, bool c, int x, int y) : t(t), c(c), p(x, y) {}
boson_::boson_(size_t t, bool c, Vec2<int> p) : t(t), c(c), p(p) {}

void boson_::setPosition(int x, int y)
{
    p[0] = x;
    p[1] = y;
}

bool boson_::operator==(const boson_ &f) const
{
    // return (this->c == f.c) && (this->t == f.t) && (p[0] == f.p[0]) && (p[1] == f.p[1]);
    return (this->c == f.c) && (this->t == f.t) && p == f.p;
}

bool boson_::sameType(const boson_ &f) const
{
    return (this->c == f.c) && (this->t == f.t);
}

Vec2<boson_> make_pair(const boson_ &a, const boson_ &b)
{
    return {{a.t, a.c, {0, 0}}, {b.t, b.c, {b.p[0] - a.p[0], b.p[1] - a.p[1]}}};
}

// Bosons_ implementations
void Bosons_::clear()
{
    terms.clear();
}

void Bosons_::add(const std::vector<boson_> &bs, Complex co)
{
    terms.emplace_back(co, bs);
}

void Bosons_::addPair(const Vec2<boson_> &pair, Complex co)
{
    add({pair[0], pair[1]}, co);
}

Bosons_ Bosons_::operator+(const term_ &term) const
{
    Bosons_ B(*this);
    B.terms.push_back(term);
    // B.simplify();
    return B;
}

Bosons_ Bosons_::operator+(const Bosons_ &Bs) const
{
    Bosons_ B(*this);
    B.terms.insert(B.terms.end(), Bs.terms.begin(), Bs.terms.end());
    // B.simplify();
    return B;
}

void Bosons_::operator+=(const Bosons_ &Bs)
{
    terms.insert(terms.end(), Bs.terms.begin(), Bs.terms.end());
    // simplify();
}

Bosons_ Bosons_::getOrder(size_t n) const
{
    Bosons_ B;
    for (size_t i = 0; i < size(); i++)
    {
        if (size(i) == n)
        {
            B.terms.push_back(terms[i]);
        }
    }
    return B;
}

Bosons_ Bosons_::grep(std::vector<size_t> t, std::vector<bool> c)
{
    size_t n = t.size();
    Bosons_ B;
    for (size_t i = 0; i < size(); i++)
    {
        auto &bos = terms[i].vec;
        if (bos.size() != n)
            continue;
        bool same = true;
        for (size_t j = 0; j < n; j++)
        {
            if (!(bos[j].c == c[j] && bos[j].t == t[j]))
            {
                same = false;
            }
        }
        if (same)
        {
            B.terms.push_back(terms[i]);
        }
    }
    return B;
}

void Bosons_::order()
{
    for (size_t i = 0; i < terms.size(); i++)
    {
        auto &bos = terms[i].vec;
        size_t m = bos.size();
        std::vector<term_> newTerms;
        for (size_t j = 0; j < m; j++)
        {
            for (size_t k = 0; k < m - j - 1; k++)
            {
                if (bos[k].c == 0 && bos[k + 1].c == 1)
                {
                    std::swap(bos[k], bos[k + 1]);
                    if (bos[k].t == bos[k + 1].t && bos[k].p == bos[k + 1].p)
                    {
                        term_ newTerm;
                        newTerm.vec = bos;
                        newTerm.vec.erase(newTerm.vec.begin() + k + 1);
                        newTerm.vec.erase(newTerm.vec.begin() + k);
                        newTerm.coe = terms[i].coe;
                        newTerms.push_back(newTerm);
                    }
                }
                else if (bos[k].c == bos[k + 1].c)
                {
                    if ((bos[k + 1].t < bos[k].t) || ((bos[k + 1].t == bos[k].t) && (bos[k + 1].p < bos[k].p)))
                    {
                        std::swap(bos[k], bos[k + 1]);
                    }
                }
            }
        }
        terms.insert(terms.end(), newTerms.begin(), newTerms.end());
    }
}

void Bosons_::shift()
{
    for (auto &term : terms)
    {
        int x, y;
        if (term.vec.size())
        {
            x = term.vec[0].p[0];
            y = term.vec[0].p[1];
            for (auto &bos : term.vec)
            {
                bos.p[0] -= x;
                bos.p[1] -= y;
            }
        }
    }
}

void Bosons_::merge()
{
    for (size_t i = 0; i < terms.size(); i++)
    {
        for (size_t j = i + 1; j < terms.size(); j++)
        {
            if (terms[i] == terms[j])
            {
                terms[i].coe += terms[j].coe;
                terms.erase(terms.begin() + j);
                j--;
            }
        }
    }
}

void Bosons_::chop(double eps)
{
    for (size_t i = 0; i < terms.size(); i++)
    {
        Complex &co = terms[i].coe;
        if (std::abs(co.real()) < eps && std::abs(co.imag()) < eps)
        {
            terms.erase(terms.begin() + i);
            i--;
            continue;
        }
        if (std::abs(co.real()) < eps)
            co.real(0);
        if (std::abs(co.imag()) < eps)
            co.imag(0);
    }
}

void Bosons_::simplify(double eps)
{
    order();
    shift();
    merge();
    chop(eps);
}

void Bosons_::symmetrize(double eps)
{
    for (auto i = terms.begin(); i != terms.end() - 1; i++)
    {
        for (auto j = i + 1; j != terms.end(); j++)
        {
            if (std::abs((*i).coe.real() - (*j).coe.real()) < eps)
                (*j).coe.real((*i).coe.real());
            if (std::abs((*i).coe.real() + (*j).coe.real()) < eps)
                (*j).coe.real(-(*i).coe.real());
            if (std::abs((*i).coe.imag() - (*j).coe.imag()) < eps)
                (*j).coe.imag((*i).coe.imag());
            if (std::abs((*i).coe.imag() + (*j).coe.imag()) < eps)
                (*j).coe.imag(-(*i).coe.imag());
        }
    }
}
// Pairs_ implementations
bool Pairs_::WickPair(const std::vector<boson_> &list, const PairList &pairList,
                      std::vector<PairList> &pairLists)
{
    if (list.size() % 2 == 1)
    {
        pairLists.clear();
        return false;
    }
    std::vector<boson_> newList;
    for (size_t i = 1; i < list.size(); i++)
    {
        PairList pairList_new = pairList;
        Pair_ a = make_pair(list[0], list[i]);
        // a.y.p[0] -= a.x.p[0];
        // a.y.p[1] -= a.x.p[1];
        // a.x.p = {0, 0};
        pairList_new.push_back(a);
        newList.clear();
        for (size_t j = 1; j < list.size(); j++)
        {
            if (j != i)
            {
                newList.push_back(list[j]);
            }
        }
        if (newList.size() == 0)
        {
            pairLists.push_back(pairList_new);
        }
        else
            WickPair(newList, pairList_new, pairLists);
    }
    return true;
}

Pairs_::Pairs_(const Bosons_ &Bs)
{
    PairList pairList;
    size_t n = Bs.size();
    pairListsTerms.reserve(n);
    for (auto &term : Bs.terms)
    {
        Term<PairList> pairLists{};
        pairLists.coe = term.coe;
        pairList.clear();
        if (term.vec.size())
        {
            WickPair(term.vec, pairList, pairLists.vec);
        }
        else
        {
            pairLists.vec.push_back(pairList);
        }
        pairListsTerms.push_back(pairLists);
    }
}

Pairs_::Pairs_(const Bosons_ &Bs, const Pair_ &pair) : Pairs_(Bs)
{
    Pair_ pairShift = make_pair(pair.x, pair.y);
    for (auto &pairListsTerm : pairListsTerms)
    {
        std::vector<PairList> pairLists;
        for (const auto &pairList : pairListsTerm.vec)
        {
            size_t num = pairList.size();
            for (size_t i = 0; i < num; i++)
            {
                if (pairList[i] == pairShift)
                {
                    PairList pairListNew = pairList;
                    pairListNew.erase(pairListNew.begin() + i);
                    pairLists.push_back(pairListNew);
                }
            }
        }
        if (pairLists.size() == 0)
        {
            pairListsTerm.coe = 0.;
        }
        pairListsTerm.vec = pairLists;
    }
}

Pairs_::Pairs_(Bosons_ Bs, boson_ b)
{
    PairList pairList;
    size_t n = Bs.size();
    pairListsTerms.reserve(n);
    for (const auto &term : Bs.terms)
    {
        Term<PairList> pairListsTerm;
        pairListsTerm.coe = term.coe;
        for (size_t i = 0; i < term.vec.size(); i++)
        {
            if (b.sameType(term.vec[i]))
            {
                pairList.clear();
                auto bos = term.vec;
                bos.erase(bos.begin() + i);
                if (bos.size())
                    WickPair(bos, pairList, pairListsTerm.vec);
                else
                    pairListsTerm.vec.push_back(pairList);
            }
        }
        pairListsTerms.push_back(pairListsTerm);
    }
}

// boson implementations
std::vector<double> boson::linspace(double s, double e)
{
    double step = (e - s) / N;
    std::vector<double> result;
    result.reserve(N + 1);
    for (size_t i = 0; i < N + 1; i++)
    {
        result.push_back(s);
        s += step;
    }
    return result;
}

boson::boson(size_t t, bool c) : t(t), c(c) {}
boson::boson(size_t t, bool c, Vec2<double> subPos) : t(t), c(c), sublatticePos(subPos) {}
boson::boson(const boson_ &br) : t(br.t), c(br.c)
{
    phase = std::shared_ptr<Complex[]>(new Complex[(N + 1) * (N + 1)]);
    if (br.p[0] == 0 && br.p[1] == 0)
    {
        std::fill(phase.get(), phase.get() + (N + 1) * (N + 1), Complex{1, 0});
        trivial = true;
        return;
    }
    auto k = linspace(-M_PI, M_PI);
    for (size_t i = 0; i < N + 1; i++)
    {
        for (size_t j = 0; j < N + 1; j++)
        {
            phase[i * (N + 1) + j] = std::exp(Complex{0., 1.} * (k[i] * br.p[0] + k[j] * br.p[1]));
        }
    }
}

boson::boson(const boson_ &br, Vec2<double> subPos) : t(br.t), c(br.c), sublatticePos(subPos)
{
    phase = std::shared_ptr<Complex[]>(new Complex[(N + 1) * (N + 1)]);
    if (br.p[0] == 0 && br.p[1] == 0)
    {
        std::fill(phase.get(), phase.get() + (N + 1) * (N + 1), Complex{1, 0});
        trivial = true;
        return;
    }
    static auto k = linspace(-M_PI, M_PI);
    for (size_t i = 0; i < N + 1; i++)
    {
        for (size_t j = 0; j < N + 1; j++)
        {
            phase[i * (N + 1) + j] = std::exp(Complex{0., 1.} * (k[i] * (br.p[0]) + k[j] * (br.p[1])));
        }
    }
}

void boson::copy(boson &other) const
{
    other.t = t;
    other.c = c;
    if ((!other.phase.get()) || other.phase.get() == phase.get())
    {
        other.phase = std::shared_ptr<Complex[]>(new Complex[(N + 1) * (N + 1)]);
    }
    std::copy(phase.get(), phase.get() + (N + 1) * (N + 1), other.phase.get());
}

// const Complex &boson::operator[](size_t i) const
// {
//     // if (!phase.get())
//     // {
//     //     throw std::runtime_error("The boson phase does not initial.");
//     // }
//     return phase[i];
// }

// const Complex &boson::operator()(size_t i, size_t j) const
// {
//     return this->operator[](i * (N + 1) + j);
// }

void boson::setPhase(const boson &other, Complex *p)
{
    phase = std::shared_ptr<Complex[]>(new Complex[(N + 1) * (N + 1)]);
    if (other.ifTrivial())
    {
        std::copy(p, p + (N + 1) * (N + 1), phase.get());
    }
    else
    {
        for (size_t i = 0; i < (N + 1) * (N + 1); i++)
        {
            phase[i] = other[i] * p[i];
        }
    }
    trivial = false;
}

long boson::use_count() const
{
    return phase.use_count();
}

bool boson::ifTrivial() const
{
    return trivial;
}

bool boson::sameType(const boson &other) const
{
    return (this->c == other.c) && (this->t == other.t);
}

Complex boson::contract(const boson &other) const
{
    Complex result{};
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            result += phase[i * (N + 1) + j] * other(N - i, N - j);
        }
    }
    return result / (double)(N * N);
}

// Bosons implementations
void Bosons::clear()
{
    terms.clear();
}

void Bosons::add(std::vector<boson> bs, Complex co)
{
    terms.emplace_back(co, bs);
}

Bosons Bosons::operator+(const term &ter) const
{
    Bosons B(*this);
    B.terms.push_back(ter);
    return B;
}

Bosons Bosons::operator+(const Bosons &Bs) const
{
    Bosons B(*this);
    B.terms.insert(B.terms.end(), Bs.terms.begin(), Bs.terms.end());
    return B;
}

void Bosons::operator+=(const Bosons &Bs)
{
    terms.insert(terms.end(), Bs.terms.begin(), Bs.terms.end());
}
void Bosons::operator*=(Complex c)
{
    for (auto &t : terms)
    {
        t.coe *= c;
    }
}

Bosons Bosons::getOrder(size_t n) const
{
    Bosons B;
    for (size_t i = 0; i < size(); i++)
    {
        if (terms[i].vec.size() == n)
        {
            B.terms.push_back(terms[i]);
        }
    }
    return B;
}

Bosons Bosons::grep(const std::vector<boson> &bs) const
{
    size_t n = bs.size();
    Bosons B;
    for (size_t i = 0; i < size(); i++)
    {
        auto &bos = terms[i].vec;
        if (bos.size() != n)
            continue;
        bool same = true;
        for (size_t j = 0; j < n; j++)
        {
            if (!(bos[j].c == bs[j].c && bos[j].t == bs[j].t))
            {
                same = false;
            }
        }
        if (same)
        {
            B.terms.push_back(terms[i]);
        }
    }
    return B;
}

void Bosons::order()
{
    for (size_t i = 0; i < terms.size(); i++)
    {
        std::vector<boson> &bos = terms[i].vec;
        size_t m = bos.size();
        std::vector<term> newTerms;
        for (size_t j = 0; j < m; j++)
        {
            for (size_t k = 0; k < m - j - 1; k++)
            {
                if (bos[k].c == 0 && bos[k + 1].c == 1)
                {
                    std::swap(bos[k], bos[k + 1]);
                    if (bos[k].t == bos[k + 1].t)
                    {
                        term newTerm;
                        newTerm.vec = bos;
                        newTerm.coe = terms[i].coe * bos[k + 1].contract(bos[k]);
                        newTerm.vec.erase(newTerm.vec.begin() + k + 1);
                        newTerm.vec.erase(newTerm.vec.begin() + k);
                        newTerms.push_back(newTerm);
                    }
                }
                else if (bos[k].c == bos[k + 1].c)
                {
                    if ((bos[k + 1].t < bos[k].t))
                    {
                        std::swap(bos[k], bos[k + 1]);
                    }
                    else if (bos[k + 1].t == bos[k].t && ((!bos[k].ifTrivial()) && bos[k + 1].ifTrivial()))
                    {
                        std::swap(bos[k], bos[k + 1]);
                    }
                }
            }
        }
        terms.insert(terms.end(), newTerms.begin(), newTerms.end());
    }
}

void Bosons::merge()
{
    for (size_t i = 0; i < terms.size(); i++)
    {
        for (size_t j = i + 1; j < terms.size(); j++)
        {
            bool same = true;
            int nontrivialNum = 0;
            size_t nontrivialIndex;
            const auto &a = terms[i].vec;
            const auto &b = terms[j].vec;
            if (a.size() != b.size())
            {
                continue;
            }
            else
            {
                for (size_t k = 0; k < a.size(); k++)
                {
                    if (!a[k].sameType(b[k]))
                    {
                        same = false;
                        break;
                    }
                    else
                    {
                        if (!(a[k].trivial && b[k].trivial))
                        {
                            nontrivialNum++;
                            nontrivialIndex = k;
                        }
                    }
                }
                if (same)
                {
                    if (nontrivialNum == 0)
                    {
                        terms[i].coe += terms[j].coe;
                        terms.erase(terms.begin() + j);
                        j--;
                    }
                    else if (nontrivialNum == 1)
                    {
                        auto &b1 = terms[i].vec[nontrivialIndex];
                        const auto &b2 = terms[j].vec[nontrivialIndex];
                        b1.trivial = false;
                        for (size_t k = 0; k < (boson::N + 1) * (boson::N + 1); k++)
                        {
                            b1.phase[k] = b1.phase[k] * terms[i].coe + b2.phase[k] * terms[j].coe;
                        }
                        terms[i].coe = 1.;
                        terms.erase(terms.begin() + j);
                        j--;
                    }
                    else
                        continue;
                }
            }
        }
    }
}

std::vector<Bosons::term> Bosons::getTerms(const std::vector<boson> &bs)
{
    std::vector<term> ts;
    size_t n = bs.size();
    for (size_t i = 0; i < size(); i++)
    {
        auto &bos = terms[i].vec;
        if (bos.size() != n)
            continue;
        bool same = true;
        for (size_t j = 0; j < n; j++)
        {
            if (!(bos[j].c == bs[j].c && bos[j].t == bs[j].t))
            {
                same = false;
            }
        }
        if (same)
        {
            ts.push_back(terms[i]);
        }
    }
    return ts;
}

Complex Bosons::constant() const
{
    Complex c{};
    for (const auto &term : terms)
    {
        if (term.vec.size() == 0)
        {
            c += term.coe;
        }
    }
    return c;
}

Bosons::Bosons(const Bosons_ &Bs)
{
    for (const auto &ter : Bs.terms)
    {
        std::vector<boson> bos;
        for (const auto &br : ter.vec)
        {
            bos.emplace_back(br);
        }
        terms.emplace_back(ter.coe, bos);
    }
    merge();
}

Bosons::Bosons(const Bosons_ &Bs, const std::vector<Vec2<double>> &subPositions)
{
    for (const auto &ter : Bs.terms)
    {
        std::vector<boson> bos;
        for (const auto &br : ter.vec)
        {
            bos.emplace_back(br, subPositions[br.t]);
        }
        terms.emplace_back(ter.coe, bos);
    }
    merge();
}
Bosons Bosons::swap(std::vector<size_t> indexs) const
{
    Bosons BB;
    std::sort(indexs.begin(), indexs.end());
    std::vector<size_t> perm = indexs;
    do
    {
        for (const auto &term : terms)
        {
            auto te = term;
            for (size_t i = 0; i < indexs.size(); i++)
            {
                te.vec[indexs[i]] = term.vec[perm[i]];
            }
            BB.terms.push_back(te);
        }
    } while (std::next_permutation(perm.begin(), perm.end()));
    return BB;
}
Bosons Bosons::swap(size_t i, size_t j) const
{
    Bosons BB;
    for (const auto &term : terms)
    {
        auto te = term;
        std::swap(te.vec[i], te.vec[j]);
        BB.terms.push_back(te);
    }
    return BB;
}
std::unique_ptr<Complex[]> Bosons::reduce(std::vector<Vec2<size_t>> ks, size_t i0, size_t j0) const
{
    boson b;
    size_t N = boson::N;
    size_t m = ks.size();
    std::unique_ptr<Complex[]> result = std::make_unique<Complex[]>(N*N);
    std::vector<size_t> reIndex;
    Vec2<size_t> sum = {0, 0};
    for (size_t in = 0; in < ks.size(); in++)
    {
        if (in != i0 && in != j0)
        {
            reIndex.push_back(in);
            sum.x += ks[in].x;
            sum.y += ks[in].y;
        }
    }
    for (auto &t : terms)
    {
        Complex co = t.coe;
        for (auto rei : reIndex)
        {
            co *= t.vec[rei](ks[rei].x, ks[rei].y);
        }
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                size_t ipq = (m * N / 2 + m * N - sum.x - i)%N;
                size_t jpq = (m * N / 2 + m * N - sum.y - j)%N;
                result[i * N + j] += co * t.vec[i0](i, j) * t.vec[j0](ipq, jpq);
            }
        }
    }
    return result;
}

// Operator overloads
Bosons_ operator*(Complex c, const Bosons_ &Bs)
{
    Bosons_ B = Bs;
    for (Term<boson_> &term : B.terms)
    {
        term.coe *= c;
    }
    return B;
}

Bosons_ operator*(const Bosons_ &Bs1, const Bosons_ &Bs2)
{
    Bosons_ B;
    for (size_t i = 0; i < Bs1.size(); i++)
    {
        for (size_t j = 0; j < Bs2.size(); j++)
        {
            auto bos = Bs1.terms[i].vec;
            bos.insert(bos.end(), Bs2.terms[j].vec.begin(), Bs2.terms[j].vec.end());
            B.terms.emplace_back(Bs1.terms[i].coe * Bs2.terms[j].coe, bos);
        }
    }
    return B;
}

Bosons operator*(const Bosons &Bk1, const Bosons &Bk2)
{
    Bosons Bk;
    for (size_t i = 0; i < Bk1.size(); i++)
    {
        for (size_t j = 0; j < Bk2.size(); j++)
        {
            auto bos = Bk1.terms[i].vec;
            bos.insert(bos.end(), Bk2.terms[j].vec.begin(), Bk2.terms[j].vec.end());
            Bk.terms.emplace_back(Bk1.terms[i].coe * Bk2.terms[j].coe, bos);
        }
    }
    return Bk;
}

bool operator==(const Term<boson_> &a, const Term<boson_> &b)
{
    if (a.vec.size() == b.vec.size())
    {
        for (size_t i = 0; i < a.vec.size(); i++)
        {
            if (!(a.vec[i] == b.vec[i]))
            {
                return false;
            }
        }
        return true;
    }
    return false;
}

// template<size_t N>
// bool operator==(const Term<boson> &a, const Term<boson> &b)
// {
//     if (a.vec.size() == b.vec.size())
//     {
//         for (size_t i = 0; i < a.vec.size(); i++)
//         {
//             if (!(a.vec[i].sameType(b.vec[i])))
//             {
//                 return false;
//             }
//         }
//         return true;
//     }
//     return false;
// }

/***************ostream function*************/
std::ostream &operator<<(std::ostream &os, const boson_ &a)
{
    if (a.c)
        os << (char)(a.t + 97) << '^' << '(' << a.p[0] << ' ' << a.p[1] << ')';
    else
    {
        os << (char)(a.t + 97) << '(' << a.p[0] << ' ' << a.p[1] << ')';
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const boson &a)
{
    if (a.c)
        os << (char)(a.t + 97) << '^';
    else
    {
        os << (char)(a.t + 97);
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const Term<std::vector<Vec2<boson_>>> &t)
{
    os << t.coe << '\n';
    for (const auto &pairs : t.vec)
    {
        for (const auto &pair : pairs)
        {
            os << pair;
        }
        os << '\n';
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const Bosons_ &Bs)
{
    for (auto i = Bs.terms.begin(); i != Bs.terms.end(); i++)
    {
        os << *i;
        if (i + 1 != Bs.terms.end())
        {
            os << '+';
        }
    }
    return os;
}

std::ostream &operator<<(std::ostream &os, const Bosons &Bs)
{
    for (auto i = Bs.terms.begin(); i != Bs.terms.end(); i++)
    {
        os << *i;
        if (i + 1 != Bs.terms.end())
        {
            os << '+';
        }
    }
    return os;
}