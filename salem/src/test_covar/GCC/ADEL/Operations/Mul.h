#ifndef MUL_H
#define MUL_H

template<typename L, typename R>
struct Mul : StaticInterface<Mul<L, R> >, BinaryOperation<L, R> {
    using BinaryOperation<L, R>::Right;
    using BinaryOperation<L, R>::Left;

    INLINE_MODE
    Mul(const L& left, const R& right) :
        BinaryOperation<L, R> (left, right, left.Value* right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        BinaryOperation<L, R>::Left. template Gradient<GradLevel> (gradData, index, baseCount, gradient * BinaryOperation<L, R>::Right);
        BinaryOperation<L, R>::Right. template Gradient<GradLevel> (gradData, index, baseCount, BinaryOperation<L, R>::Left * gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        BinaryOperation<L, R>::Left.Gradient(gradData, index, baseCount, gradient * BinaryOperation<L, R>::Right.Value);
        BinaryOperation<L, R>::Right.Gradient(gradData, index, baseCount, BinaryOperation<L, R>::Left.Value * gradient);
    }
};

template<typename L, typename R>INLINE_MODE
const Mul<L, R> operator * (const StaticInterface<L>& left, const StaticInterface<R>& right) {
    return Mul<L, R> (left.Self(), right.Self());
}

template<typename R>
struct Mul<real, R> : StaticInterface<Mul<real, R> >, BinaryOperation<real, R> {
    using BinaryOperation<real, R>::Right;
    using BinaryOperation<real, R>::Left;

    INLINE_MODE
    Mul(const real& left, const R& right) :
        BinaryOperation<real, R> (left, right, left* right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, Left * gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, Left * gradient);
    }
};

template<typename R>INLINE_MODE
const Mul<real, R> operator * (const real& left, const StaticInterface<R>& right) {
    return Mul<real, R> (left, right.Self());
}

template<typename L>
struct Mul<L, real> : StaticInterface<Mul<L, real> >, BinaryOperation<L, real> {
    using BinaryOperation<L, real>::Right;
    using BinaryOperation<L, real>::Left;

    INLINE_MODE
    Mul(const L& left, const real& right) :
        BinaryOperation<L, real> (left, right, left.Value* right) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Left. template Gradient<GradLevel> (gradData, index, baseCount, gradient * Right);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Left.Gradient(gradData, index, baseCount, gradient * Right);
    }
};

template<typename L>INLINE_MODE
const Mul<L, real> operator * (const StaticInterface<L>& left, const real& right) {
    return Mul<L, real> (left.Self(), right);
}

#endif  // MUL_H
