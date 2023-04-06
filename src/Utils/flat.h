#ifndef flat_h
#define flat_h

#include <type_traits>
#include <vector>

#include "base.h"

namespace detail {
    template <class Key>
    struct DefaultIndex {
        inline size_t operator()(const Key key) const {
            return static_cast<size_t>(key);
        }
    };
}

template <class Key, class T, class IndexFn = detail::DefaultIndex<Key>>
class Flat {
public:
    using self_type       = Flat;
    using storage_type    = std::vector<T>;
    using key_type        = Key;
    using reference       = T&;
    using const_reference = const T&;
    using iterator        = typename storage_type::iterator;


    LINEAR_DESIGN_INLINE iterator begin() {
        return data_.begin();
    }

    LINEAR_DESIGN_INLINE iterator end() {
        return data_.end();
    }

    LINEAR_DESIGN_INLINE bool empty() const {
        return false;
    }
    
    LINEAR_DESIGN_INLINE void reserve(const size_t n) {
        data_.reserve(n);
    }

    LINEAR_DESIGN_INLINE void resize(const size_t n) {
        data_.resize(n);
    }

    template <enable_if_t<!std::is_integral<key_type>::value, int> = 0>
	LINEAR_DESIGN_INLINE reference operator[](size_t index) {
        return data_[index];
    }

    template <enable_if_t<!std::is_integral<key_type>::value, int> = 0>
    LINEAR_DESIGN_INLINE const_reference operator[](size_t index) const {
        return data_[index];
    }

    LINEAR_DESIGN_INLINE reference operator[](key_type key) {
        return data_[index_(key)];
    }
    
    LINEAR_DESIGN_INLINE const_reference operator[](key_type key) const {
        return data_[index_(key)];
    }

    LINEAR_DESIGN_INLINE size_t size() const {
        return data_.size();
    }
    
private:
    IndexFn index_;
    storage_type data_;
};

#endif /* flat_h */
