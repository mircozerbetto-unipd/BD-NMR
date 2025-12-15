#ifndef FLOOR_H
#define FLOOR_H

template<typename R>
struct Floor : StaticInterface<Floor<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Floor(const R& right) :
        UnaryOperation<R> (right, floor(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {}
};

template<typename R>INLINE_MODE
const Floor<R> floor(const StaticInterface<R>& right) {
    return Floor<R> (right.Self());
}

#endif  // FLOOR_H
