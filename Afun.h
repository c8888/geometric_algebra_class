//
// Created by c8888 on 28.07.16.
//

#ifndef GEOMETRIC_ALGEBRA_CLASSES_BFUN_H
#define GEOMETRIC_ALGEBRA_CLASSES_BFUN_H


#include <include/vsr/space/vsr_ega3D_types.h>


class Afun {
public:
    virtual vsr::ega::Biv value(double t) = 0; //THE PSEUDOSCALAR PART SHOULD BE INCLUDED HERE.
    virtual vsr::ega::Rot exactSol(double t, vsr::ega::Rot R0) {return vsr::ega::Rot(0,0,0,0);} //when not specified, exact solution is set to zero!!!

};

class Afun1: public Afun {
public:
    vsr::ega::Biv value(double t);
    vsr::ega::Rot exactSol (double t, vsr::ega::Rot R0);
};

class Afun2: public Afun {
public:
    vsr::ega::Biv value(double t);
    vsr::ega::Rot exactSol (double t, vsr::ega::Rot R0);
private:
    const double B0 = 1;
    const double B1 = 10;
    const double gamma = 1;
    const double w = 1;
    vsr::ega::Pss I = vsr::ega::Pss(1);
};

#endif //GEOMETRIC_ALGEBRA_CLASSES_BFUN_H
