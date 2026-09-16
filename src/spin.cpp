#include "spin.h"
#include <stdexcept>
#include <set>
#include <iostream>

double spin::S = 0.5;

// spin implementations
spin::spin(size_t t, int d) : t(t)
{
    if (d < 3 && d >= 0)
        direction = d;
    else
        throw std::out_of_range("Spin direction must be 0, 1, 2 (corresponding x, y, z)");
}

spin::spin(size_t t, int d, int x, int y) : spin(t, d)
{
    pos.x = x;
    pos.y = y;
}

spin::spin(size_t t, int d, Vec2<int> p) : spin(t, d)
{
    pos = p;
}

bool spin::operator==(const spin &other) const
{
    return (direction == other.direction && t == other.t && pos == other.pos);
}

bool spin::sameType(const spin &other) const
{
    return (direction == other.direction && t == other.t);
}

bool spin::operator<(const spin &other) const
{
    return (t < other.t) || (t == other.t && direction < other.direction) || (sameType(other) && pos < other.pos);
}

// Spins implementations
void Spins::add(const term &t)
{
    terms.push_back(t);
}

void Spins::add(const Spins &Ss, Complex coe)
{
    for (const auto &term : Ss.terms)
    {
        terms.push_back({coe * term.coe, term.vec});
    }
}

void Spins::add(spin s, Complex coe)
{
    terms.push_back({coe, {s}});
}

size_t Spins::size() const
{
    return terms.size();
}

Spins::term &Spins::operator[](size_t i)
{
    return terms[i];
}

const Spins::term &Spins::operator[](size_t i) const
{
    return terms[i];
}

void Spins::order()
{
    for (size_t i = 0; i < terms.size(); i++)
    {
        auto &vec = terms[i].vec;
        size_t m = vec.size();
        for (size_t j = 0; j < m; j++)
        {
            for (size_t k = 0; k < m - j - 1; k++)
            {
                if (vec[k + 1] < vec[k])
                {
                    std::swap(vec[k], vec[k + 1]);
                }
            }
        }
    }
}

void Spins::shift()
{
    for (auto &term : terms)
    {
        int x, y;
        if (term.vec.size())
        {
            x = term.vec[0].pos[0];
            y = term.vec[0].pos[1];
            for (auto &vec : term.vec)
            {
                vec.pos[0] -= x;
                vec.pos[1] -= y;
            }
        }
    }
}

void Spins::merge()
{
    for (size_t i = 0; i < terms.size(); i++)
    {
        for (size_t j = i + 1; j < terms.size(); j++)
        {
            bool same = true;
            const auto &a = terms[i].vec;
            const auto &b = terms[j].vec;
            if (a.size() != b.size())
                continue;
            else
            {
                for (size_t k = 0; k < a.size(); k++)
                {
                    if (!(a[k] == b[k]))
                    {
                        same = false;
                        break;
                    }
                }
                if (same)
                {
                    terms[i].coe += terms[j].coe;
                    terms.erase(terms.begin() + j);
                    j--;
                }
            }
        }
    }
}

void Spins::chop(double eps)
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

void Spins::simplify(double eps)
{
    order();
    shift();
    merge();
    chop(eps);
}
size_t Spins::getNu() const
{
    std::set<size_t> flavors;
    for (auto term : terms)
    {
        for (auto spin : term.vec)
        {
            flavors.insert(spin.t);
        }
    }
    size_t f = 0;
    for (auto i : flavors)
    {
        if (i != f)
        {
            std::cout << "The types of spin in Spins are: ";
            for (auto j : flavors)
            {
                std::cout << "j ";
            }
            std::cout << std::endl;
            std::runtime_error("The type of spin in Spins must be 0 to Nu-1");
        }
        f++;
    }
    return flavors.size();
}
Spins Spins::grep(spin s) const
{
    Spins SR;
    for (auto te : terms)
    {
        for (size_t i = 0; i < te.vec.size(); i++)
        {
            if (te.vec[i].sameType(s))
            {
                term ss = te;
                ss.vec.erase(ss.vec.begin() + i);
                for (auto &v : ss.vec)
                {
                    v.pos[0] -= te.vec[i].pos[0];
                    v.pos[1] -= te.vec[i].pos[1];
                }
                SR.terms.push_back(ss);
            }
        }
    }
    return SR;
}

/*****************************************************/
Spins operator+(const Spins &a, const Spins &b)
{
    Spins Ss(a);
    Ss.add(b);
    return Ss;
}

Spins operator*(Complex c, const Spins &s)
{
    Spins Ss(s);
    for (auto &term : Ss.terms)
    {
        term.coe *= c;
    }
    return Ss;
}

Spins interaction(const std::vector<size_t> &ts, const std::vector<int> &ds, const std::vector<Vec2<int>> &ps, Complex couple)
{
    if (ts.size() != ds.size() || ts.size() != ps.size())
    {
        throw std::runtime_error("Spin interaction input is invalid.");
    }
    Spins Ss;
    Term<spin> term;
    term.coe = couple;
    for (size_t i = 0; i < ts.size(); i++)
    {
        term.vec.emplace_back(ts[i], ds[i], ps[i]);
    }
    Ss.terms.push_back(term);
    return Ss;
}

Spins interaction(size_t t1, int d1, size_t t2, int d2, Vec2<int> pos, Complex couple)
{
    Spins Ss;
    Ss.terms.emplace_back(couple, std::vector<spin>{{t1, d1}, {t2, d2, pos}});
    return Ss;
}

Spins Heisenberg(size_t t1, size_t t2, Vec2<int> pos, Complex couple)
{
    Spins Ss;
    for (int i = 0; i < 3; i++)
    {
        Ss.terms.emplace_back(couple, std::vector<spin>{{t1, i}, {t2, i, pos}});
    }
    return Ss;
}

Spins fieldInteraction(size_t t, Vec3<double> field)
{
    Spins Ss;
    for (int d = 0; d < 3; d++)
    {
        spin sa(t, d);
        Ss.terms.emplace_back(field[d], std::vector<spin>{sa});
    }
    return Ss;
}

/*****************************************************/
std::ostream &operator<<(std::ostream &os, const spin &s)
{
    os << 'S' << char('a' + s.t) << char('x' + s.direction) << '(' << s.pos.x << ',' << s.pos.y << ')';
    return os;
}

std::ostream &operator<<(std::ostream &os, const Spins &Ss)
{
    bool sign = false;
    for (const auto &term : Ss.terms)
    {
        if (!sign)
        {
            sign = true;
        }
        else
        {
            os << '+';
        }
        os << term.coe;
        for (const auto &s : term.vec)
        {
            os << s;
        }
    }
    return os;
}