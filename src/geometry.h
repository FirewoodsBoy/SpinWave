#if !defined(GEOMETRY_H)
#define GEOMETRY_H
#define _USE_MATH_DEFINES
#include <stdexcept>
#include <vector>
#include <complex>
typedef std::complex<double> Complex;
template <typename T>
class Vec2
{
public:
    T x, y;
    Vec2() : x{}, y{} {};
    Vec2(T x, T y) : x(x), y(y) {};
    T &operator[](size_t i)
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            throw std::out_of_range("Vec2 index out of range");
    }
    const T &operator[](size_t i) const
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else
            throw std::out_of_range("Vec2 index out of range");
    }
    bool operator==(const Vec2<T> &other) const
    {
        return (x == other.x) && (y == other.y);
    }
    bool operator<(const Vec2<T> &other) const
    {
        return (x < other.x) || (x == other.x && y < other.y);
    }
};
template <typename T>
class Vec3
{
public:
    T x, y, z;
    Vec3() : x{}, y{}, z{} {};
    Vec3(T x, T y, T z) : x(x), y(y), z(z) {};
    T &operator[](size_t i)
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else if (i == 2)
            return z;
        else
            throw std::out_of_range("Vec3 index out of range");
    }
    const T &operator[](size_t i) const
    {
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        else if (i == 2)
            return z;
        else
            throw std::out_of_range("Vec3 index out of range");
    }
    bool operator==(const Vec3<T> &other) const
    {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }
    bool operator<(const Vec3<T> &other) const
    {
        return (x < other.x) || (x == other.x && (y < other.y || (y == other.y && z < other.z)));
    }
    Vec3<T> cross(const Vec3<T> &other) const
    {
        Vec3<T> result{
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x};
        return result;
    }
    T dot(const Vec3<T> &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }
    Vec3 getNormalizedVec3() const
    {
        double l = std::sqrt(x * x + y * y + z * z);
        return Vec3<T>{x / l, y / l, z / l};
    }
};
template <typename T>
struct Term
{
    Complex coe{};
    std::vector<T> vec;
    Term() = default;
    Term(Complex coe, const std::vector<T> &vec) : coe(coe), vec(vec) {};
};
template <typename T>
struct Term2
{
    Complex coe{};
    Vec2<T> vec;
    Term2() = default;
    Term2(Complex coe, Vec2<T> vec) : coe(coe), vec(vec) {};
};
namespace std
{
    template <typename T>
    struct hash<Vec2<T>>
    {
        size_t operator()(const Vec2<T> &v) const noexcept
        {
            return hash<T>{}(v.x) ^ (hash<T>{}(v.y) << 1);
        }
    };
}
namespace std
{
    template <typename T>
    struct hash<Vec3<T>>
    {
        size_t operator()(const Vec3<T> &v) const noexcept
        {
            return hash<T>{}(v.x) ^ (hash<T>{}(v.y) << 1) ^ (hash<T>{}(v.z) << 2);
        }
    };
}
template <typename T>
struct Term3
{
    Complex coe{};
    Vec3<T> vec;
    Term3() = default;
    Term3(Complex coe, Vec3<T> vec) : coe(coe), vec(vec) {};
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const Vec3<T> &s)
{
    os << '[' << s[0] << ',' << s[1] << ',' << s[2] << ']';
    return os;
}
#endif // GEOMETRY_H
