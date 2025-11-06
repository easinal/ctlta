#pragma once

#include <iterator>
#include <cassert>

// Range of iterators that satisfies iterability and allows random access.
// User needs to make sure that begin and end are well-defined.
template<typename IteratorT>
class IteratorRange {
public:
    using ValueType = typename std::iterator_traits<IteratorT>::value_type;
    using Iterator = IteratorT;

    IteratorRange(const IteratorT startingPoint, const IteratorT endingPoint)
            : startingPoint(startingPoint), endingPoint(endingPoint) {}

    Iterator begin() const noexcept {
        return startingPoint;
    }

    Iterator end() const noexcept {
        return endingPoint;
    }

    const ValueType &operator[](const int col) const {
        assert(col >= 0);
        assert(col < endingPoint - startingPoint);
        return startingPoint[col];
    }

    size_t size() const noexcept {
        return endingPoint - startingPoint;
    }

    bool empty() const noexcept {
        return endingPoint == startingPoint;
    }

private:
    IteratorT startingPoint;
    IteratorT endingPoint;
};


template<typename T, template<typename> typename VecT = std::vector>
class ConstantVectorRange : public IteratorRange<typename VecT<T>::const_iterator> {
    using IteratorRange<typename VecT<T>::const_iterator>::IteratorRange;
};