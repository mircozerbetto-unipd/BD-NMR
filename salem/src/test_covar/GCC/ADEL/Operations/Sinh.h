#ifndef SINH_H
#define SINH_H

template<typename R>
struct Sinh : StaticInterface<Sinh<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Sinh(const R& right) :
        UnaryOperation<R> (right, sinh(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, gradient * cosh(Right));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, gradient * cosh(Right.Value));
    }
};

template<typename R>INLINE_MODE
const Sinh<R> sinh(const StaticInterface<R>& right) {
    return Sinh<R> (right.Self());
}

#endif  // SINH_H
