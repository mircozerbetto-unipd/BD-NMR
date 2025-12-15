#ifndef OPERATIONS_H
#define OPERATIONS_H

struct Operation {
    const real Value;

    INLINE_MODE
    explicit Operation(const real& value) : Value(value) {}
};

template<typename L, typename R>
struct BinaryOperation : Operation {
    const L& Left;
    const R& Right;

    INLINE_MODE
    BinaryOperation(const L& left, const R& right, const real& value) : Operation(value), Left(left), Right(right) {}

    template<unsigned Level, unsigned Nvars>INLINE_MODE
    void GetTemporaries(std::vector<const ADVariable<Level, Nvars>*>& temporaries) const {
        Left.template GetTemporaries<Level, Nvars> (temporaries);
        Right.template GetTemporaries<Level, Nvars> (temporaries);
    }
};

template<typename R>
struct BinaryOperation<real, R> : Operation {
    const real& Left;
    const R& Right;

    INLINE_MODE
    BinaryOperation(const real& left, const R& right, const real& value) : Operation(value), Left(left), Right(right) {}

    template<unsigned Level, unsigned Nvars>INLINE_MODE
    void GetTemporaries(std::vector<const ADVariable<Level, Nvars>*>& temporaries) const {
        Right.template GetTemporaries<Level, Nvars> (temporaries);
    }
};

template<typename L>
struct BinaryOperation<L, real> : Operation {
    const L& Left;
    const real& Right;

    INLINE_MODE
    BinaryOperation(const L& left, const real& right, const real& value) : Operation(value), Left(left), Right(right) {}

    template<unsigned Level, unsigned Nvars>INLINE_MODE
    void GetTemporaries(std::vector<const ADVariable<Level, Nvars>*>& temporaries) const {
        Left.template GetTemporaries<Level, Nvars> (temporaries);
    }
};

template<typename R>
struct UnaryOperation : Operation {
    const R& Right;

    INLINE_MODE
    UnaryOperation(const R& right, const real& value) : Operation(value), Right(right) {}

    template<unsigned Level, unsigned Nvars>INLINE_MODE
    void GetTemporaries(std::vector<const ADVariable<Level, Nvars>*>& temporaries) const {
        Right.template GetTemporaries<Level, Nvars> (temporaries);
    }
};

template<typename A, typename B, typename C, typename D>
struct CartesianOperation : Operation {
    A& AIN;
    const B& BIN;
    const C& CIN;
    const D& DIN;

    INLINE_MODE
    CartesianOperation(A& aIN, const B& bIN, const C& cIN, const D& dIN, const real& value) : Operation(value), AIN(aIN), BIN(bIN), CIN(cIN), DIN(dIN) {}

    template<unsigned Level, unsigned Nvars>INLINE_MODE
    void GetTemporaries(std::vector<const ADVariable<Level, Nvars>*>& temporaries) const {
        AIN.template GetTemporaries<Level, Nvars> (temporaries);
        BIN.template GetTemporaries<Level, Nvars> (temporaries);
        CIN.template GetTemporaries<Level, Nvars> (temporaries);
        DIN.template GetTemporaries<Level, Nvars> (temporaries);
    }
};

template<typename C, typename D>
struct CartesianOperation<molec*, int, C, D> : Operation {
    molec*& AIN;
    const int&    BIN;
    const C&      CIN;
    const D&      DIN;

    INLINE_MODE
    CartesianOperation(molec*& aIN, const int& bIN, const C& cIN, const D& dIN, const real& value) : Operation(value), AIN(aIN), BIN(bIN), CIN(cIN), DIN(dIN) {}

    template<unsigned Level, unsigned Nvars>INLINE_MODE
    void GetTemporaries(std::vector<const ADVariable<Level, Nvars>*>& temporaries) const {
        CIN.template GetTemporaries<Level, Nvars> (temporaries);
        DIN.template GetTemporaries<Level, Nvars> (temporaries);
    }
};

#endif  // OPERATIONS_H
