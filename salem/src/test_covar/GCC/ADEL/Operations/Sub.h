#ifndef SUB_H
#define SUB_H

template<typename L, typename R>
struct Sub : StaticInterface<Sub<L, R> >, BinaryOperation<L, R> {
    using BinaryOperation<L, R>::Right;
    using BinaryOperation<L, R>::Left;

    INLINE_MODE
    Sub(const L& left, const R& right) :
        BinaryOperation<L, R> (left, right, left.Value - right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        BinaryOperation<L, R>::Left.template Gradient<GradLevel> (gradData, index, baseCount, gradient);
        BinaryOperation<L, R>::Right.template Gradient<GradLevel> (gradData, index, baseCount, -gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        BinaryOperation<L, R>::Left.template Gradient(gradData, index, baseCount, gradient);
        BinaryOperation<L, R>::Right.template Gradient(gradData, index, baseCount, -gradient);
    }
};

template<typename L, typename R>INLINE_MODE
const Sub<L, R> operator - (const StaticInterface<L>& left, const StaticInterface<R>& right) {
    return Sub<L, R> (left.Self(), right.Self());
}

template<typename R>
struct Sub<real, R> : StaticInterface<Sub<real, R> >, BinaryOperation<real, R> {
    using BinaryOperation<real, R>::Right;
    using BinaryOperation<real, R>::Left;

    INLINE_MODE
    Sub(const real& left, const R& right) :
        BinaryOperation<real, R> (left, right, left - right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right.template Gradient<GradLevel> (gradData, index, baseCount, -gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.template Gradient(gradData, index, baseCount, -gradient);
    }
};

template<typename R>INLINE_MODE
const Sub<real, R> operator - (const real& left, const StaticInterface<R>& right) {
    return Sub<real, R> (left, right.Self());
}

template<typename L>
struct Sub<L, real> : StaticInterface<Sub<L, real> >, BinaryOperation<L, real> {
    using BinaryOperation<L, real>::Right;
    using BinaryOperation<L, real>::Left;

    INLINE_MODE
    Sub(const L& left, const real& right) :
        BinaryOperation<L, real> (left, right, left.Value - right) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Left. template Gradient<GradLevel> (gradData, index, baseCount, gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Left.Gradient(gradData, index, baseCount, gradient);
    }
};

template<typename L>INLINE_MODE
const Sub<L, real> operator - (const StaticInterface<L>& left, const real& right) {
    return Sub<L, real> (left.Self(), right);
}

#endif  // SUB_H
