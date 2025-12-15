#ifndef EXP_H
#define EXP_H

template<typename R>
struct Exp : StaticInterface<Exp<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;
    using UnaryOperation<R>::Operation::Value;

    INLINE_MODE
    Exp(const R& right) :
        UnaryOperation<R> (right, exp(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, gradient * (*this));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, gradient * Value);
    }
};

template<typename R>INLINE_MODE
const Exp<R> exp(const StaticInterface<R>& right) {
    return Exp<R> (right.Self());
}

#endif  // EXP_H
