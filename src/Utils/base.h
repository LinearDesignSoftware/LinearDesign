#ifndef base_h
#define base_h

#include <type_traits>
#include <utility>
#include <memory>
#include <string>
#include <vector>
#include <sstream>

#if defined(__GNUC__) || defined(__clang__)
#define LINEAR_DESIGN_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define LINEAR_DESIGN_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: function deprecated")
#define LINEAR_DESIGN_DEPRECATED
#endif

#if defined(__GNUC__) || defined(__clang__)
#define LINEAR_DESIGN_INLINE inline __attribute__((always_inline))
#else
#define LINEAR_DESIGN_INLINE inline
#endif

#define LINEAR_DESIGN_CACHELINE 64

template <bool B, typename T=void>
using enable_if_t = typename std::enable_if<B, T>::type;

template <typename T, enable_if_t<std::is_integral<T>::value, int> = 0>
std::ostream& operator<< (std::ostream& out, const std::pair<T, T>& rhs) {
    out << "(" << rhs.first << ", " << rhs.second << ")";
    return out;
}

template <typename T, enable_if_t<std::is_integral<T>::value, int> = 0>
std::ostream& operator<< (std::ostream& out, const std::vector<std::pair<T, T>>& rhs) {
    out << "[";
    for (size_t i = 0; i < rhs.size(); ++i) {
        out << rhs[i];
        if (i < rhs.size() - 1) out << ",";
    }
    out << "]";
    return out;
}

namespace LinearDesign {

namespace util {
    std::vector<std::string> split(const std::string &s, char delim) {
        std::vector<std::string> result;
        std::stringstream ss(s);
        std::string item;
        while (getline(ss, item, delim)) 
            result.push_back(item);
        return result;
    }

    template <typename T>
    constexpr T value_min() {
        static_assert(std::is_integral<T>::value ||
            std::is_floating_point<T>::value, "Int or float required.");
        return std::numeric_limits<T>::lowest();
    }

    template <typename T>
    constexpr T value_max() {
        static_assert(std::is_integral<T>::value ||
            std::is_floating_point<T>::value, "Int or float required.");
        return std::numeric_limits<T>::max();
    }
} /* util */


// template <bool...> struct is_any;
// template <> struct is_any<> : std::false_type {};
// template <bool First, bool... Rest> struct is_any<First, Rest...> {
//     constexpr static bool value = First || is_any<Rest...>::value;
// };

struct hash_pair_pair {
    template <class T1, class T2, class T3>
    size_t operator()(const std::pair<std::pair<T1, T2>, T3>& p) const {
        auto hash1 = std::hash<T1>{}(p.first.first);
        auto hash2 = std::hash<T2>{}(p.first.second);
        auto hash3 = std::hash<T3>{}(p.second);
        return hash1 ^ hash2 ^ hash3;
    }
};

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

}


namespace Hash {
    template <class T>
    LINEAR_DESIGN_INLINE size_t hash_combine(size_t left_seed, const T& right) {
        return left_seed ^ (std::hash<T>{}(right) << 1);
    }

    template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
    struct TupleHashImpl {
        static size_t impl(size_t seed, const Tuple& tuple) {
            size_t h = hash_combine(seed, std::get<Index>(tuple));
            return TupleHashImpl<Tuple, Index-1>::impl(h, tuple);
        }
    };

    template <class Tuple>
    struct TupleHashImpl<Tuple, 0> {
        static size_t impl(size_t seed, const Tuple& tuple) {
            return hash_combine(seed, std::get<0>(tuple));
        }
    };
}


template <class... Ts>
struct std::hash<std::tuple<Ts...>> {
    size_t operator()(const std::tuple<Ts...>& ts) const {
        return Hash::TupleHashImpl<std::tuple<Ts...>>::impl(0, ts);
    }
};

template <class T1, class T2>
struct std::hash<std::pair<T1, T2>> {
    size_t operator()(const std::pair<T1, T2>& p) const {
        size_t h = std::hash<T1>{}(p.first);
        return Hash::hash_combine(h, p.second);
    }
};


#endif