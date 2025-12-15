#ifndef COSH_H
#define COSH_H

template<typename R>
struct Cosh : StaticInterface<Cosh<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Cosh(const R& right) :
        UnaryOperation<R> (right, cosh(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, gradient * sinh(Right));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, gradient * sinh(Right.Value));
    }
};

template<typename R>INLINE_MODE
const Cosh<R> cosh(const StaticInterface<R>& right) {
    return Cosh<R> (right.Self());
}

#endif  // COSH_H
