#ifndef UNARYMINUS_H
#define UNARYMINUS_H

template<typename R>
struct UnaryMinus : StaticInterface<UnaryMinus<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    UnaryMinus(const R& right) :
        UnaryOperation<R> (right, -right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, -gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, -gradient);
    }
};

template<typename R>INLINE_MODE
const UnaryMinus<R> operator - (const StaticInterface<R>& right) {
    return UnaryMinus<R> (right.Self());
}

#endif  // UNARYMINUS_H
