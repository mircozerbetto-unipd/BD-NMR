#ifndef COS_H
#define COS_H

template<typename R>
struct Cos : StaticInterface<Cos<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Cos(const R& right) :
        UnaryOperation<R> (right, cos(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, -gradient * sin(Right));
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, -gradient * sin(Right.Value));
    }
};

template<typename R>INLINE_MODE
const Cos<R> cos(const StaticInterface<R>& right) {
    return Cos<R> (right.Self());
}

#endif  // COS_H
