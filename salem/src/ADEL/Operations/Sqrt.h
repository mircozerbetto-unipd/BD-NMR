#ifndef SQRT_H
#define SQRT_H

template<typename R>
struct Sqrt : StaticInterface<Sqrt<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;
    using UnaryOperation<R>::Operation::Value;

    INLINE_MODE
    Sqrt(const R& right) :
        UnaryOperation<R> (right, sqrt(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, 0.5 * gradient / (*this));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, 0.5 * gradient / Value);
    }
};

template<typename R>INLINE_MODE
const Sqrt<R> sqrt(const StaticInterface<R>& right) {
    return Sqrt<R> (right.Self());
}

#endif  // SQRT_H
