#ifndef LOG10_H
#define LOG10_H

const real BaseLog = log(10.0);
const real InverseBaseLog = 1.0 / BaseLog;

template<typename R>
struct Log10 : StaticInterface<Log10<R> >, UnaryOperation<R> {
    using UnaryOperation<R>::Right;

    INLINE_MODE
    Log10(const R& right) :
        UnaryOperation<R> (right, log10(right.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        Right. template Gradient<GradLevel> (gradData, index, baseCount, InverseBaseLog * gradient / Right);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        Right.Gradient(gradData, index, baseCount, InverseBaseLog * gradient / Right.Value);
    }
};

template<typename R>INLINE_MODE
const Log10<R> log10(const StaticInterface<R>& right) {
    return Log10<R> (right.Self());
}

#endif  // LOG10_H
