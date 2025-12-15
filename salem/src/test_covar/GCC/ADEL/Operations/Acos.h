#ifndef ACOS_H
#define ACOS_H

template<typename R>
struct Acos : StaticInterface<Acos<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Acos(const R& right) :
        UnaryOperation<R> (right, acos(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. Gradient<GradLevel> (gradData, index, baseCount, -gradient / sqrt(1.0 - sqr(Right)));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right. Gradient(gradData, index, baseCount, -gradient / sqrt(1.0 - sqr(Right.Value)));
    }
};

template<typename R>INLINE_MODE
const Acos<R> acos(const StaticInterface<R>& right) {
    return Acos<R> (right.Self());
}

#endif  // ACOS_H
