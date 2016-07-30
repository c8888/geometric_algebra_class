#include <iostream>
#include "EQ.h"
#include "numAlgo.h"
#include "Afun.h"
#include "data.h"

using namespace std;
using namespace vsr::ega;

int main() {
    EulerMethodRotor* emr = new EulerMethodRotor(new Afun1); //how to create new object. It has Afun1 as A(t) and EulerMethodRotor as numerical method.
    RungeKutta4thMethodRotor* rkmr = new RungeKutta4thMethodRotor(new Afun2);
    RodriguesFormula* rodrConst = new RodriguesFormula(new Afun1);
    RodriguesFormula* rodrNMR = new RodriguesFormula(new Afun2);
    AdamsMulton* AdmConst = new AdamsMulton(new Afun1);
    AdamsMulton* AdmConst2 = new AdamsMulton(new Afun2);


    //CALCULATE THE ROTOR AND SPIN IN TIME


    int N = 40000;
    double dt=1E-12;
    Vec S0(0,0,1);
    Rot R0(1,0,0,0);

    emr->calc(dt, N, S0, R0);
    emr->timing();
    emr->saveTR("TRfirst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    emr->saveTS("TSfirt.dat");

    rodrConst->calc(dt, N, S0, R0);
    rodrConst->timing();
    rodrConst->saveTR("TRrodrConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    rodrConst->saveTS("TSrodrConst.dat");

    AdmConst->calc(dt, N, S0, R0);
    AdmConst->timing();
    AdmConst->saveTR("TRAdamsConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    AdmConst->saveTS("TSAdamsConst.dat");


    //CALCULATE ERROR WITH ERROR MEASUREMENT METHOD 1 (change in modulus of spin)
    double dtmin = 1E-22;
    double dtmax = 1E-8;
    double logarithm_base = 1.1; //increment od dt for next step, i.e. dt *= logarithm_base

    emr->calcError1(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr->saveError("Error1Euler.dat");

    rodrConst->calcError1(dtmin, dtmax, N, logarithm_base, S0, R0);
    rodrConst->saveError("rodrConstError.dat");

    AdmConst->calcError1(dtmin, dtmax, N, logarithm_base, S0, R0);
    AdmConst->saveError("Error1AdamsConst.dat");


    return 0;
}