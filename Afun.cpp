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
    const double gamma = 2.675E8;
    Pss I(1);
    Vec B(0, gamma, 0);
    return B*I; //ALREADY INCLUDES THE PSEUDOSCALAR FACTOR
}

Rot Afun1::exactSol(double t, vsr::ega::Rot R0){
    const double gamma = 2.675E8;
    Pss I(1);
    Vec B(0,gamma,0);
    Rot Rtmp = Sca(cos((I*B*t).norm())) + (I*B*t)*(1/(I*B*t).norm())*sin((I*B*t).norm()); //THIS MEANS INITIAL ROTOR SHOULD ALWAYS BE (1,0,0,0) FOR THIS METHOD
    return Rtmp;
}

Biv Afun2::value(double t) {
    const double gamma = 2.675E8;
    Pss I(1);
    Vec B(0, gamma, 0);
    return B*I;
}

Rot Afun2::exactSol(double t, vsr::ega::Rot R0){
    const double gamma = 2.675E8;
    Pss I(1);
    Vec B(0,gamma,0);
    Rot Rtmp = Sca(cos((I*B*t).norm())) + (I*B*t)*(1/(I*B*t).norm())*sin((I*B*t).norm());
    return Rtmp;
}