//
// Created by c8888 on 28.07.16.
//

#include "Afun.h"

//FILE INCLUDES DEFINITIONS OF A(t) - BIVECTORS
//CALL CONSECUTIVE FUNCTIONS Afun2, Afun3 AND SO ON
//EACH FUNCTION Afun* CONTAINS ELEMENTS VALUE AND EXACTSOL. EXACTSOL SHOULD BE SET IF THE ANALYTICAL SOLUTION EXISTS.
//IF EXACTSOL() IS NOT SET, IT IS BY DEFUALT Rot(0,0,0,0)

using namespace std;
using namespace vsr::ega;

Biv Afun1::value(double t) {
    //const double gamma = 2.675E8;

    Pss I(1);
    Vec B(0, gamma/2, 0);
    return B*I; //ALREADY INCLUDES THE PSEUDOSCALAR FACTOR
}

Rot Afun1::exactSol(double t, vsr::ega::Rot R0){
    Pss I(1);
    Vec B(0,gamma/2,0);
    Rot Rtmp = Sca(cos((B*I).norm()*t)) + (B*I)*Sca((1/(B*I).norm())*sin((I*B).norm()*t));
    return Rtmp;
}

Biv Afun2::value(double t) {
    return Sca(gamma/2.)*Vec(B1*cos(w*t),B1*sin(w*t),B0)*Pss(1);
}

Rot Afun2::exactSol(double t, vsr::ega::Rot R0){
    Rot Rext1 = Sca(cos(-t*w/2))+Vec(0,0,1)*Pss(1)*Sca(sin(-t*w/2));
    Biv P = Sca(0.5*t)*Vec(gamma*B1,0,(gamma*B0+w))*Pss(1);
    Rot Rext2 = Sca(cos(P.norm()))+P*Sca(sin(P.norm())/P.norm());
    return Rext1*Rext2;
}

