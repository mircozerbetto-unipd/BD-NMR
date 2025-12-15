#ifndef DIV_H
#define DIV_H

template<typename L, typename R>
struct Div : StaticInterface<Div<L, R> >, BinaryOperation<L, R> {
    using BinaryOperation<L, R>::Right;
    using BinaryOperation<L, R>::Left;

    INLINE_MODE
    Div(const L& left, const R& right) :
        BinaryOperation<L, R> (left, right, left.Value / right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Left. template Gradient<GradLevel> (gradData, index, baseCount, gradient / Right);
        Right. template Gradient<GradLevel> (gradData, index, baseCount, - (gradient * Left) / sqr(Right));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Left.Gradient(gradData, index, baseCount, gradient / Right.Value);
        Right.Gradient(gradData, index, baseCount, - (gradient * Left.Value) / sqr(Right.Value));
    }
};

template<typename L, typename R>INLINE_MODE
const Div<L, R> operator / (const StaticInterface<L>& left, const StaticInterface<R>& right) {
    return Div<L, R> (left.Self(), right.Self());
}

template<typename R>
struct Div<real, R> : StaticInterface<Div<real, R> >, BinaryOperation<real, R> {
    using BinaryOperation<real, R>::Right;
    using BinaryOperation<real, R>::Left;

    INLINE_MODE
    Div(const real& left, const R& right) :
        BinaryOperation<real, R> (left, right, left / right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, - (gradient * Left) / sqr(Right));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, - (gradient *  Left) / sqr(Right.Value));
    }
};

template<typename R>INLINE_MODE
const Div<real, R> operator / (const real& left, const StaticInterface<R>& right) {
    return Div<real, R> (left, right.Self());
}

template<typename L>
struct Div<L, real> : StaticInterface<Div<L, real> >, BinaryOperation<L, real> {
    using BinaryOperation<L, real>::Right;
    using BinaryOperation<L, real>::Left;

    INLINE_MODE
    Div(const L& left, const real& right) :
        BinaryOperation<L, real> (left, right, left.Value / right) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Left. template Gradient<GradLevel> (gradData, index, baseCount, gradient / Right);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Left.Gradient(gradData, index, baseCount, gradient / Right);
    }
};

template<typename L>INLINE_MODE
const Div<L, real> operator / (const StaticInterface<L>& left, const real& right) {
    return Div<L, real> (left.Self(), right);
}

#endif  // DIV_H
