#ifndef ATAN2_H
#define ATAN2_H

template<typename L, typename R>
struct Atan2 : StaticInterface<Atan2<L, R> >, BinaryOperation<L, R> {
    using BinaryOperation<L, R>::Right;
    using BinaryOperation<L, R>::Left;

    INLINE_MODE
    Atan2(const L& left, const R& right) :
        BinaryOperation<L, R> (left, right, atan2(left.Value, right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        BinaryOperation<L, R>::Left. template Gradient<GradLevel> (gradData, index, baseCount, gradient * BinaryOperation<L, R>::Right / (sqr(BinaryOperation<L, R>::Left) + sqr(BinaryOperation<L, R>::Right)));
        BinaryOperation<L, R>::Right. template Gradient<GradLevel> (gradData, index, baseCount, -gradient * BinaryOperation<L, R>::Left / (sqr(BinaryOperation<L, R>::Left) + sqr(BinaryOperation<L, R>::Right)));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        BinaryOperation<L, R>::Left.Gradient(gradData, index, baseCount, gradient * BinaryOperation<L, R>::Right.Value / (sqr(BinaryOperation<L, R>::Left.Value) + sqr(BinaryOperation<L, R>::Right.Value)));
        BinaryOperation<L, R>::Right.Gradient(gradData, index, baseCount, -gradient * BinaryOperation<L, R>::Left.Value / (sqr(BinaryOperation<L, R>::Left.Value) + sqr(BinaryOperation<L, R>::Right.Value)));
    }
};

template<typename L, typename R>INLINE_MODE
const Atan2<L, R> atan2(const StaticInterface<L>& left, const StaticInterface<R>& right) {
    return Atan2<L, R> (left.Self(), right.Self());
}

template<typename R>
struct Atan2<real, R> : StaticInterface<Atan2<real, R> >, BinaryOperation<real, R> {
    using BinaryOperation<real, R>::Right;
    using BinaryOperation<real, R>::Left;

    INLINE_MODE
    Atan2(const real& left, const R& right) :
        BinaryOperation<real, R> (left, right, atan2(left, right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, -gradient * Left / (sqr(Left) + sqr(Right)));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, -gradient * Left / (sqr(Left) + sqr(Right.Value)));
    }
};

template<typename R>INLINE_MODE
const Atan2<real, R> atan2(const real& left, const StaticInterface<R>& right) {
    return Atan2<real, R> (left, right.Self());
}

template<typename L>
struct Atan2<L, real> : StaticInterface<Atan2<L, real> >, BinaryOperation<L, real> {
    using BinaryOperation<L, real>::Right;
    using BinaryOperation<L, real>::Left;

    INLINE_MODE
    Atan2(const L& left, const real& right) :
        BinaryOperation<L, real> (left, right, atan2(left.Value, right)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Left. template Gradient<GradLevel> (gradData, index, baseCount, gradient * Right / (sqr(Left) + sqr(Right)));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Left.Gradient(gradData, index, baseCount, gradient * Right / (sqr(Left.Value) + sqr(Right)));
    }
};

template<typename L>INLINE_MODE
const Atan2<L, real> atan2(const StaticInterface<L>& left, const real& right) {
    return Atan2<L, real> (left.Self(), right);
}

#endif  // ATAN2_H
