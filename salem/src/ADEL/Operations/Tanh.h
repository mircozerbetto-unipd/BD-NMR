#ifndef TANH_H
#define TANH_H

template<typename R>
struct Tanh : StaticInterface<Tanh<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;
    using UnaryOperation<R>::Operation::Value;

    INLINE_MODE
    Tanh(const R& right) :
        UnaryOperation<R> (right, tanh(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, gradient * (1.0 - sqr(*this)));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, gradient * (1.0 - sqr(Value)));
    }
};

template<typename R>INLINE_MODE
const Tanh<R> tanh(const StaticInterface<R>& right) {
    return Tanh<R> (right.Self());
}

#endif  // TANH_H
