#ifndef ABS_H
#define ABS_H

template<typename R>
struct Abs : StaticInterface<Abs<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Abs(const R& right) :
        UnaryOperation<R> (right, abs(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. Gradient<GradLevel> (gradData, index, baseCount, gradient * Right.Value / abs(Right.Value));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right. Gradient(gradData, index, baseCount, gradient * Right.Value / abs(Right.Value));
    }
};

template<typename R>INLINE_MODE
const Abs<R> abs(const StaticInterface<R>& right) {
    return Abs<R> (right.Self());
}

template<typename R>INLINE_MODE
const Abs<R> fabs(const StaticInterface<R>& right) {
    return Abs<R> (right.Self());
}

INLINE_MODE
real sign(const real& right) {
    return right / abs(right);
}

#endif  // ABS_H
