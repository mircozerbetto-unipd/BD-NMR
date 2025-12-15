#ifndef ADHELPER_H
#define ADHELPER_H

template<unsigned Level, unsigned Nvars>
struct ADVariable;

static unsigned GetDataSize(const unsigned level, const unsigned nvars) {
    unsigned totalSize = 0;
    unsigned currentSize = 1;
    for(unsigned i = 0; i < level; i++) {
        currentSize *= nvars;
        totalSize += currentSize;
    }
    return totalSize;
}

template<unsigned Level, unsigned Nvars, unsigned CropLevel>
struct Cropper {
    static void Crop(const real from[], real to[], unsigned fromIndex, unsigned toIndex, unsigned tmpCount) {
        const unsigned totalCount = Nvars + tmpCount;
        const unsigned toOffset = GetDataSize(Level - CropLevel, Nvars) + toIndex;
        const unsigned fromOffset = GetDataSize(Level - CropLevel, totalCount) + fromIndex;
        for(unsigned i = 0; i < Nvars; i++) {
            to[toOffset + i] = from[fromOffset + i];
            Cropper < Level, Nvars, CropLevel - 1 >::Crop(from, to, (fromIndex + i) * totalCount, (toIndex + i) * (Nvars), tmpCount);
        }
    }
};

template<unsigned Level, unsigned Nvars>
struct Cropper<Level, Nvars, 1> {
    static void Crop(const real from[], real to[], unsigned fromIndex, unsigned toIndex, unsigned tmpCount) {
        const unsigned totalCount = Nvars + tmpCount;
        const unsigned toOffset = GetDataSize(Level - 1, Nvars) + toIndex;
        const unsigned fromOffset = GetDataSize(Level - 1, totalCount) + fromIndex;
        for(unsigned i = 0; i < Nvars; i++) {
            to[toOffset + i] = from[fromOffset + i];
        }
    }
};

template<unsigned Level, unsigned Nvars, unsigned InjectionLevel>
struct Injector {
    static void Inject(const std::vector<const ADVariable<Level, Nvars>*>& dependencies, const real from[], real to[]) {
        Injector < Level, Nvars, InjectionLevel - 1 >::Inject(dependencies, from, to);
    }
};

template<unsigned Level, unsigned Nvars>
struct Injector<Level, Nvars, 2> {
    static void Inject(const std::vector<const ADVariable<Level, Nvars>*>& dependencies, const real from[], real to[]) {
        const unsigned tmpCount = static_cast<unsigned>(dependencies.size());
        const unsigned totalCount = Nvars + tmpCount;
        for(unsigned x = 0; x < tmpCount; x++) {
            const unsigned xID = dependencies[x]->ID;
            const unsigned d2X_ = totalCount + xID * totalCount;
            const real* xData = dependencies[x]->Data;
            for(unsigned j = 0; j < Nvars; j++) {
                const unsigned d2j_ = Nvars + j * Nvars;
                real sum = from[d2X_ + j];
                for(unsigned y = 0; y < tmpCount; y++) {
                    sum += from[d2X_ + dependencies[y]->ID] * dependencies[y]->Data[j];
                }
                for(unsigned i = j; i < Nvars; i++) {
                    to[d2j_ + i] +=  sum * xData[i] + from[totalCount + i * totalCount + xID] * xData[j] + from[xID] * xData[d2j_ + i];
                }
            }
        }
        for(unsigned j = 0; j < Nvars; j++) {
            const unsigned d2j_ = Nvars + j * Nvars;
            for(unsigned i = 0; i < j; i++) {
                to[d2j_ + i] = to[Nvars + i * Nvars + j];
            }
        }
        Injector<Level, Nvars, 1>::Inject(dependencies, from, to);
    }
};

template<unsigned Level, unsigned Nvars>
struct Injector<Level, Nvars, 1> {
    static void Inject(const std::vector<const ADVariable<Level, Nvars>*>& dependencies, const real from[], real to[]) {
        const unsigned tmpCount = static_cast<unsigned>(dependencies.size());
        for(unsigned x = 0; x < tmpCount; x++) {
            const unsigned xID = dependencies[x]->ID;
            const real* xData = dependencies[x]->Data;
            for(unsigned i = 0; i < Nvars; i++) {
                to[i] += from[xID] * xData[i];
            }
        }
    }
};



#endif  // ADHELPER_H
