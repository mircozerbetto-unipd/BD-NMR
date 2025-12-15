#ifndef CARTESIANFD_H
#define CARTESIANFD_H

// coord: 0-based coordinate ID in range 0 <= coord < 3Natoms
// coordinates: x = 0, y = 1, z = 2
// example: coordinate y of atomID 8: coord = (atomID-1) * 3 * Natoms + 1
//
// qid: 0-based internal coordinate ID in range 0 <= coord <= 3Natoms
// internal coordinates: d = 0, theta = 1, phi = 2
// example: coordinate phi of atomID 5: qid = (atomID - 1) * 3 * Natoms + 2
double getCartesianFirstDerivative(molec* mol, int coord, real q1, real q2);

// coord: 0-based coordinate ID in range 0 <= coord < 3Natoms
// coordinates: x = 0, y = 1, z = 2
// example: coordinate y of atomID 8: coord = (atomID-1) * 3 * Natoms + 1
//
// qid: 0-based internal coordinate ID in range 0 <= coord <= 3Natoms
// internal coordinates: d = 0, theta = 1, phi = 2
// example: coordinate phi of atomID 5: qid = (atomID - 1) * 3 * Natoms + 2
double getCartesianSecondDerivative(molec* mol, int coord, real q1, real q2);

template<typename A, typename B, typename C, typename D>
struct CartesianFirstDerivative : StaticInterface<CartesianFirstDerivative<A, B, C, D> >, CartesianOperation<A, B, C, D> {
    using CartesianOperation<A, B, C, D>::AIN;
    using CartesianOperation<A, B, C, D>::BIN;
    using CartesianOperation<A, B, C, D>::CIN;
    using CartesianOperation<A, B, C, D>::DIN;
    using CartesianOperation<A, B, C, D>::Operation::Value;

    INLINE_MODE
    CartesianFirstDerivative(A& mol, const B& coor, const C& q1, const D& q2) :
        CartesianOperation<A, B, C, D> (mol, coor, q1, q2, getCartesianFirstDerivative(mol.Value, coor.Value, q1.Value, q2.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        CartesianOperation<A, B, C, D>::AIN. template Gradient<GradLevel> (gradData, index, baseCount, (*this));
        CartesianOperation<A, B, C, D>::BIN. template Gradient<GradLevel> (gradData, index, baseCount, (*this));
        CartesianOperation<A, B, C, D>::CIN. template Gradient<GradLevel> (gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<A, B, C, D>::AIN, CartesianOperation<A, B, C, D>::BIN, CartesianOperation<A, B, C, D>::CIN, CartesianOperation<A, B, C, D>::CIN) * gradient);
        CartesianOperation<A, B, C, D>::DIN. template Gradient<GradLevel> (gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<A, B, C, D>::AIN, CartesianOperation<A, B, C, D>::BIN, CartesianOperation<A, B, C, D>::CIN, CartesianOperation<A, B, C, D>::DIN) * gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        CartesianOperation<A, B, C, D>::AIN.Gradient(gradData, index, baseCount, 0.0);
        CartesianOperation<A, B, C, D>::BIN.Gradient(gradData, index, baseCount, 0.0);
        CartesianOperation<A, B, C, D>::CIN.Gradient(gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<A, B, C, D>::AIN.Value, CartesianOperation<A, B, C, D>::BIN.Value, CartesianOperation<A, B, C, D>::CIN, CartesianOperation<A, B, C, D>::CIN) * gradient);
        CartesianOperation<A, B, C, D>::DIN.Gradient(gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<A, B, C, D>::AIN.Value, CartesianOperation<A, B, C, D>::BIN.Value, CartesianOperation<A, B, C, D>::CIN, CartesianOperation<A, B, C, D>::DIN) * gradient);
    }
};

template<typename C, typename D>
struct CartesianFirstDerivative<molec*, int, C, D> : StaticInterface<CartesianFirstDerivative<molec*, int, C, D> >, CartesianOperation<molec*, int, C, D> {
    using CartesianOperation<molec*, int, C, D>::AIN;
    using CartesianOperation<molec*, int, C, D>::BIN;
    using CartesianOperation<molec*, int, C, D>::CIN;
    using CartesianOperation<molec*, int, C, D>::DIN;
    using CartesianOperation<molec*, int, C, D>::Operation::Value;

    INLINE_MODE
    CartesianFirstDerivative(molec*& mol, const int& coord, const C& q1, const D& q2) :
        CartesianOperation<molec*, int, C, D> (mol, coord, q1, q2, getCartesianFirstDerivative(mol, coord, q1.Value, q2.Value)) {}

    template<unsigned GradLevel, typename T>INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const T& gradient) const {
        CartesianOperation<molec*, int, C, D>::CIN. template Gradient<GradLevel> (gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<molec*, int, C, D>::AIN, CartesianOperation<molec*, int, C, D>::BIN, CartesianOperation<molec*, int, C, D>::CIN.Value, CartesianOperation<molec*, int, C, D>::CIN.Value) * gradient);
        CartesianOperation<molec*, int, C, D>::DIN. template Gradient<GradLevel> (gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<molec*, int, C, D>::AIN, CartesianOperation<molec*, int, C, D>::BIN, CartesianOperation<molec*, int, C, D>::CIN.Value, CartesianOperation<molec*, int, C, D>::DIN.Value) * gradient);
    }

    INLINE_MODE
    void Gradient(real gradData[], unsigned index, const unsigned baseCount, const real& gradient) const {
        CartesianOperation<molec*, int, C, D>::CIN.Gradient(gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<molec*, int, C, D>::AIN, CartesianOperation<molec*, int, C, D>::BIN, CartesianOperation<molec*, int, C, D>::CIN.Value, CartesianOperation<molec*, int, C, D>::CIN.Value) * gradient);
        CartesianOperation<molec*, int, C, D>::DIN.Gradient(gradData, index, baseCount, getCartesianSecondDerivative(CartesianOperation<molec*, int, C, D>::AIN, CartesianOperation<molec*, int, C, D>::BIN, CartesianOperation<molec*, int, C, D>::CIN.Value, CartesianOperation<molec*, int, C, D>::DIN.Value) * gradient);
    }
};

template<typename C, typename D>INLINE_MODE
const CartesianFirstDerivative<molec*, int, C, D> cartesianFirstDerivative(molec*& mol, const int& coord, const StaticInterface<C>& q1, const StaticInterface<D>& q2) {
    return CartesianFirstDerivative<molec*, int, C, D> (mol, coord, q1.Self(), q2.Self());
}

#endif  // CARTESIANFD_H
