#ifndef STATICINTERFACE_H
#define STATICINTERFACE_H

template<typename StaticBase>
struct StaticInterface {
    INLINE_MODE
    const StaticBase& Self() const {
        return static_cast<const StaticBase&>(*this);
    }
};

#endif  // STATICINTERFACE_H
