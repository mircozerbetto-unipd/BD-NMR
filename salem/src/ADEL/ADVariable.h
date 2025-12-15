#ifndef ADVARIABLE_H
#define ADVARIABLE_H

#include <string.h>
#include <stdlib.h>
#include <vector>
#include "ADHelper.h"

template<unsigned Level, typename T, unsigned GradLevel> struct ADVariableSpec;

template<unsigned Level, unsigned Nvars>
struct ADVariable : StaticInterface<ADVariable<Level, Nvars> > {
    static const unsigned CurrentSize = ADVariable < Level - 1, Nvars >::CurrentSize * Nvars;
    static const unsigned TotalSize = ADVariable < Level - 1, Nvars >::TotalSize + CurrentSize;
    static const unsigned ExtensionFactor = 2;
    static const unsigned ExtendedCurrentSize = ADVariable < Level - 1, Nvars >::ExtendedCurrentSize * Nvars * ExtensionFactor;
    static const unsigned ExtendedTotalSize = ADVariable < Level - 1, Nvars >::ExtendedTotalSize + ExtendedCurrentSize;

    real Value;
    real Data[TotalSize];
    real DataCopy[TotalSize];
    static real ExtendedData[ExtendedTotalSize];
    std::vector<const ADVariable<Level, Nvars>*> InjectionData;

    mutable unsigned ID;
    mutable bool Base;

    template<typename R>
    ADVariable(const StaticInterface<R>& expression) : Value(expression.Self().Value), ID(-1), Base(false) {
        expression.Self().GetTemporaries(InjectionData);
        unsigned tmpCount = static_cast<unsigned>(InjectionData.size());
        if(tmpCount == 0) {
            memset(Data, 0, TotalSize * sizeof(real));
            expression.Self().template Gradient<Level> (Data, 0, Nvars, 1.0);
        } else {
            unsigned tmpSize = GetDataSize(Level, Nvars + tmpCount);
            //ExtendedData = static_cast<real*>(malloc(tmpSize * sizeof(real)));
            memset(ExtendedData, 0, tmpSize * sizeof(real));
            expression.Self().template Gradient<Level> (ExtendedData, 0, Nvars + tmpCount, 1.0);
            Cropper<Level, Nvars, Level>::Crop(ExtendedData, Data, 0, 0, tmpCount);
            Injector<Level, Nvars, Level>::Inject(InjectionData, ExtendedData, Data);
            for(unsigned i = 0; i < tmpCount; i++) {
                InjectionData[i]->Base = false;
            }
        }
    }

    template<typename R>
    ADVariable& operator= (const StaticInterface<R>& expression) {
        InjectionData.clear();
        expression.Self().GetTemporaries(InjectionData);
        unsigned tmpCount = static_cast<unsigned>(InjectionData.size());
        if(tmpCount == 0) {
            memset(Data, 0, TotalSize * sizeof(real));
            expression.Self().template Gradient<Level> (Data, 0, Nvars, 1.0);
        } else {
            unsigned tmpSize = GetDataSize(Level, Nvars + tmpCount);
            //ExtendedData = static_cast<real*>(realloc(ExtendedData, tmpSize * sizeof(real)));
            memset(ExtendedData, 0, tmpSize * sizeof(real));
            (expression.Self()).template Gradient<Level> (ExtendedData, 0, Nvars + tmpCount, 1.0);
            Cropper<Level, Nvars, Level>::Crop(ExtendedData, DataCopy, 0, 0, tmpCount);
            Injector<Level, Nvars, Level>::Inject(InjectionData, ExtendedData, DataCopy);
            memcpy(Data, DataCopy, TotalSize * sizeof(real));
            for(unsigned i = 0; i < tmpCount; i++) {
                InjectionData[i]->Base = false;
            }
        }
        Value = expression.Self().Value;
        return *this;
    }

    ADVariable() : Value(0.0), ID(-1), Base(false) {
        memset(Data, 0, TotalSize * sizeof(real));
    }

    ADVariable(const real& value) : Value(value), ID(-1), Base(false) {
        memset(Data, 0, TotalSize * sizeof(real));
    }

    ADVariable& operator= (const real& value) {
        Value = value;
        memset(Data, 0, TotalSize * sizeof(real));
        return *this;
    }

    ADVariable(const ADVariable& other) : Value(other.Value), ID(-1), Base(false) {
        if(other.Base) {
            memset(Data, 0, TotalSize * sizeof(real));
            Data[other.ID] = 1.0;
        } else {
            memcpy(Data, other.Data, TotalSize * sizeof(real));
        }
    }

    ADVariable& operator= (const ADVariable& other) {
        Value = other.Value;
        ID = -1;
        Base = false;
        if(other.Base) {
            memset(Data, 0, TotalSize * sizeof(real));
            Data[other.ID] = 1.0;
        } else {
            memcpy(Data, other.Data, TotalSize * sizeof(real));
        }
        return *this;
    }

    ~ADVariable() {
    }

    template<unsigned TempLevel, unsigned TempNvars>INLINE_MODE
    void GetTemporaries(std::vector<const ADVariable<TempLevel, TempNvars>*>& temporaries) const {
        if(!Base) {
            ID = TempNvars + static_cast<unsigned>(temporaries.size());
            temporaries.push_back(this);
            Base = true;
        }
    }

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        ADVariableSpec<Level, T, GradLevel>::Compute(gradData, index, ID, baseCount, gradient);
    }

    template<unsigned GradLevel>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        gradData[GetDataSize(Level - GradLevel, baseCount) + index + ID] += gradient;
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        gradData[GetDataSize(Level - 1, baseCount) + index + ID] += gradient;
    }
};

template<unsigned Level, unsigned Nvars>
real ADVariable<Level, Nvars>::ExtendedData[ExtendedTotalSize] = {0};

template<unsigned Level, typename T, unsigned GradLevel>
struct ADVariableSpec {
    INLINE_MODE
    static void Compute(real gradData[], unsigned index, unsigned ID, const unsigned baseCount, const T& gradient) {
        gradData[GetDataSize(Level - GradLevel, baseCount) + index + ID] += gradient.Self().Value;
        gradient.Self().template Gradient < GradLevel - 1 > (gradData, (index + ID) * baseCount, baseCount, 1.0);
    }
};

template<unsigned Level, typename T>
struct ADVariableSpec<Level, T, 2> {
    INLINE_MODE
    static void Compute(real gradData[], unsigned index, unsigned ID, const unsigned baseCount, const T& gradient) {
        gradData[GetDataSize(Level - 2, baseCount) + index + ID] += gradient.Self().Value;
        gradient.Self().Gradient(gradData, (index + ID) * baseCount, baseCount, 1.0);
    }
};

template<unsigned Level, typename T>
struct ADVariableSpec<Level, T, 1> {
    INLINE_MODE
    static void Compute(real gradData[], unsigned index, unsigned ID, const unsigned baseCount, const T& gradient) {
        gradData[GetDataSize(Level - 1, baseCount) + index + ID] += gradient.Self().Value;
    }
};

template<unsigned Level, typename T>
struct ADVariableSpec<Level, T, 0> {
    INLINE_MODE
    static void Compute(real gradData[], unsigned index, unsigned ID, const unsigned baseCount, const T& gradient) {}
};

template<unsigned Nvars>
struct ADVariable<0, Nvars> {
    static const unsigned CurrentSize = 1;
    static const unsigned TotalSize = 0;
    static const unsigned ExtendedCurrentSize = 1;
    static const unsigned ExtendedTotalSize = 0;
};

#endif  // ADVARIABLE_H
