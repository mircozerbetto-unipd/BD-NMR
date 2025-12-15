#ifndef SQR_H
#define SQR_H

template<typename R>
struct Sqr : StaticInterface<Sqr<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Sqr(const R& right) :
        UnaryOperation<R> (right, right.Value* right.Value) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, 2.0 * gradient * Right);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, 2.0 * gradient * Right.Value);
    }
};

template<typename R>INLINE_MODE
const Sqr<R> sqr(const StaticInterface<R>& right) {
    return Sqr<R> (right.Self());
}

INLINE_MODE
real sqr(const real& right) {
    return right * right;
}

#endif  // SQR_H
