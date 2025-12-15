#ifndef POW_H
#define POW_H

template<typename B, typename E>
struct Pow : StaticInterface<Pow<B, E> >, BinaryOperation<B, E> {
    using BinaryOperation<B, E>::Right;
    using BinaryOperation<B, E>::Left;
    using BinaryOperation<B, E>::Operation::Value;

    INLINE_MODE
    Pow(const B& base, const E& exponent) :
        BinaryOperation<B, E> (base, exponent, pow(base.Value, exponent.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        BinaryOperation<B, E>::Left. template Gradient<GradLevel> (gradData, index, baseCount, (*this) * BinaryOperation<B, E>::Right * gradient / BinaryOperation<B, E>::Left);
        BinaryOperation<B, E>::Right. template Gradient<GradLevel> (gradData, index, baseCount, (*this) * log(BinaryOperation<B, E>::Left) * gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        BinaryOperation<B, E>::Left.Gradient(gradData, index, baseCount, Value * BinaryOperation<B, E>::Right.Value * gradient / BinaryOperation<B, E>::Left.Value);
        BinaryOperation<B, E>::Right.Gradient(gradData, index, baseCount, Value * log(BinaryOperation<B, E>::Left.Value) * gradient);
    }
};

template<typename B, typename E>INLINE_MODE
const Pow<B, E> pow(const StaticInterface<B>& base, const StaticInterface<E>& exponent) {
    return Pow<B, E> (base.Self(), exponent.Self());
}

template<typename E>
struct Pow<real, E> : StaticInterface<Pow<real, E> >, BinaryOperation<real, E> {
    using BinaryOperation<real, E>::Right;
    using BinaryOperation<real, E>::Left;
    using BinaryOperation<real, E>::Operation::Value;

    INLINE_MODE
    Pow(const real& base, const E& exponent) :
        BinaryOperation<real, E> (base, exponent, pow(base, exponent.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        BinaryOperation<real, E>::Right. template Gradient<GradLevel> (gradData, index, baseCount, (*this) * log(BinaryOperation<real, E>::Left) * gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        BinaryOperation<real, E>::Right.Gradient(gradData, index, baseCount, Value * log(BinaryOperation<real, E>::Left) * gradient);
    }
};

template<typename E>INLINE_MODE
const Pow<real, E> pow(const real& base, const StaticInterface<E>& exponent) {
    return Pow<real, E> (base, exponent.Self());
}

template<typename B>
struct Pow<B, real> : StaticInterface<Pow<B, real> >, BinaryOperation<B, real> {
    using BinaryOperation<B, real>::Right;
    using BinaryOperation<B, real>::Left;
    using BinaryOperation<B, real>::Operation::Value;

    INLINE_MODE
    Pow(const B& base, const real& exponent) :
        BinaryOperation<B, real> (base, exponent, pow(base.Value, exponent)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        BinaryOperation<B, real>::Left. template Gradient<GradLevel> (gradData, index, baseCount, (*this) *  BinaryOperation<B, real>::Right * gradient /  BinaryOperation<B, real>::Left);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        BinaryOperation<B, real>::Left.Gradient(gradData, index, baseCount, Value * BinaryOperation<B, real>::Right * gradient / BinaryOperation<B, real>::Left.Value);
    }
};

template<typename B>INLINE_MODE
const Pow<B, real> pow(const StaticInterface<B>& base, const real& exponent) {
    return Pow<B, real> (base.Self(), exponent);
}

#endif  // POW_H
