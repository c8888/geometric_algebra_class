#include <iostream>
#include "EQ.h"
#include "numAlgo.h"
#include "Afun.h"
#include "data.h"

using namespace std;
using namespace vsr::ega;

int main() {
    EulerMethodRotor* euler = new EulerMethodRotor(new Afun2); //how to create new object. It has Afun1 as A(t) and EulerMethodRotor as numerical method.
    EulerMethodRotorRescaling1* emr1 = new EulerMethodRotorRescaling1(new Afun2);
    EulerMethodRotorRescaling2* emr2 = new EulerMethodRotorRescaling2(new Afun2);
    EulerMethodRotorRescaling3* emr3 = new EulerMethodRotorRescaling3(new Afun2);
    RungeKutta4thMethodRotor* runge = new RungeKutta4thMethodRotor(new Afun2);
    RodriguesFormula* rodr = new RodriguesFormula(new Afun2);
    AdamsMulton* Adm = new AdamsMulton(new Afun2);
    EulerMethodConvent* eulerConvent = new EulerMethodConvent(new Afun2);
    EulerMethodConventRescaling* eulerConventRescaling = new EulerMethodConventRescaling(new Afun2);

    //CALCULATE THE ROTOR AND SPIN IN TIME


    int N = 400;
    double dt=1E-12;
    Vec S0(0,0,1);
    Rot R0(1,0,0,0);

/*
    euler->calc(dt, N, S0, R0);
    euler->timing();
    euler->saveTR("TRfirst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    euler->saveTS("TSfirt.dat");

    rodrConst->calc(dt, N, S0, R0);
    rodrConst->timing();
    rodrConst->saveTR("TRrodrConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    rodrConst->saveTS("TSrodrConst.dat");

    AdmConst->calc(dt, N, S0, R0);
    AdmConst->timing();
    AdmConst->saveTR("TRAdamsConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    AdmConst->saveTS("TSAdamsConst.dat");
    */


    //CALCULATE ERROR WITH ERROR MEASUREMENT
    double dtmin = 1E-23;
    double dtmax = 1E-1;
    double logarithm_base = 1.1; //increment od dt for next step, i.e. dt *= logarithm_base

    euler->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    euler->saveError("Error2eulerConst.dat");

    emr1->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr1->saveError("Error2euler1Const.dat");

    emr2->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr2->saveError("Error2euler2Const.dat");


    emr3->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr3->saveError("Error2euler3Const.dat");


    runge->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    runge->saveError("Error2rungeConst.dat");

    rodr->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    rodr->saveError("Error2rodrConst.dat");

    Adm->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    Adm->saveError("Error2adamsConst.dat");

    eulerConvent->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    eulerConvent->saveError("Error2eulerConventConst.dat");

    eulerConventRescaling->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    eulerConventRescaling->saveError("Error2eulerConventRescalingConst.dat");
     

    delete euler;
    delete emr1;
    delete emr2;
    delete emr3;
    delete runge;
    delete rodr;
    delete Adm;
    delete eulerConvent;
    delete eulerConventRescaling;

    return 0;
}