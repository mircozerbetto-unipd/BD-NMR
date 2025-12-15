#ifndef CEIL_H
#define CEIL_H

template<typename R>
struct Ceil : StaticInterface<Ceil<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Ceil(const R& right) :
        UnaryOperation<R> (right, ceil(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {}
};

template<typename R>INLINE_MODE
const Ceil<R> ceil(const StaticInterface<R>& right) {
    return Ceil<R> (right.Self());
}

#endif  // CEIL_H
