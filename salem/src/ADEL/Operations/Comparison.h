#ifndef COMPARISON_H
#define COMPARISON_H

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator== (const ADVariable<Level, Nvars>& left, const ADVariable<Level, Nvars>& right) {
    return left.Value == right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator!= (const ADVariable<Level, Nvars>& left, const ADVariable<Level, Nvars>& right) {
    return left.Value != right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator< (const ADVariable<Level, Nvars>& left, const ADVariable<Level, Nvars>& right) {
    return left.Value < right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator> (const ADVariable<Level, Nvars>& left, const ADVariable<Level, Nvars>& right) {
    return left.Value > right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator<= (const ADVariable<Level, Nvars>& left, const ADVariable<Level, Nvars>& right) {
    return left.Value <= right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator>= (const ADVariable<Level, Nvars>& left, const ADVariable<Level, Nvars>& right) {
    return left.Value >= right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator== (const ADVariable<Level, Nvars>& left, const real& right) {
    return left.Value == right;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator!= (const ADVariable<Level, Nvars>& left, const real& right) {
    return left.Value != right;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator< (const ADVariable<Level, Nvars>& left, const real& right) {
    return left.Value < right;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator> (const ADVariable<Level, Nvars>& left, const real& right) {
    return left.Value > right;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator<= (const ADVariable<Level, Nvars>& left, const real& right) {
    return left.Value <= right;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator>= (const ADVariable<Level, Nvars>& left, const real& right) {
    return left.Value >= right;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator== (const real& left, const ADVariable<Level, Nvars>& right) {
    return left == right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator!= (const real& left, const ADVariable<Level, Nvars>& right) {
    return left != right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator< (const real& left, const ADVariable<Level, Nvars>& right) {
    return left < right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator> (const real& left, const ADVariable<Level, Nvars>& right) {
    return left > right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator<= (const real& left, const ADVariable<Level, Nvars>& right) {
    return left <= right.Value;
}

template<unsigned Level, unsigned Nvars> INLINE_MODE
bool operator>= (const real& left, const ADVariable<Level, Nvars>& right) {
    return left >= right.Value;
}

#endif  // COMPARISON_H
