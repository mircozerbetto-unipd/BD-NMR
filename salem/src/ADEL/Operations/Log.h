#ifndef LOG_H
#define LOG_H

template<typename R>
struct Log : StaticInterface<Log<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Log(const R& right) :
        UnaryOperation<R> (right, log(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, gradient / Right);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, gradient / Right.Value);
    }
};

template<typename R>INLINE_MODE
const Log<R> log(const StaticInterface<R>& right) {
    return Log<R> (right.Self());
}

#endif  // LOG_H
