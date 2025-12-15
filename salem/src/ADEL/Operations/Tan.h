#ifndef TAN_H
#define TAN_H

template<typename R>
struct Tan : StaticInterface<Tan<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;
    using UnaryOperation<R>::Operation::Value;

    INLINE_MODE
    Tan(const R& right) :
        UnaryOperation<R> (right, tan(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, gradient * (1.0 + sqr(*this)));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, gradient * (1.0 + sqr(Value)));
    }
};

template<typename R>INLINE_MODE
const Tan<R> tan(const StaticInterface<R>& right) {
    return Tan<R> (right.Self());
}

#endif  // TAN_H
