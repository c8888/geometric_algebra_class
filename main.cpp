#include <iostream>
#include "EQ.h"
#include "numAlgo.h"
#include "Afun.h"
#include "data.h"


#define cE calcError1
#define af Afun2


using namespace std;
using namespace vsr::ega;

int main() {
    EulerMethodRotor* euler = new EulerMethodRotor(new af); //how to create new object. It has Afun1 as A(t) and EulerMethodRotor as numerical method.
    EulerMethodRotorRescaling1* emr1 = new EulerMethodRotorRescaling1(new af);
    EulerMethodRotorRescaling2* emr2 = new EulerMethodRotorRescaling2(new af);
    EulerMethodRotorRescaling3* emr3 = new EulerMethodRotorRescaling3(new af);
    RungeKutta4thMethodRotor* runge = new RungeKutta4thMethodRotor(new af); //works ok
    RodriguesFormula* rodr = new RodriguesFormula(new af); //works
    AdamsMulton* Adm = new AdamsMulton(new af); //works
    EulerMethodConvent* eulerConvent = new EulerMethodConvent(new af);
    EulerMethodConventRescaling* eulerConventRescaling = new EulerMethodConventRescaling(new af);
    Milne* milne = new Milne(new af);
    MilneCorrected* milnecorr = new MilneCorrected(new af);
    Adams* Ad = new Adams(new af); //works
    ExactSolution* ex = new ExactSolution(new af); //works

    //CALCULATE THE ROTOR AND SPIN IN TIME


    int N = 10000;
    double dt=1E-1;
    Vec S0(0,0,1);
    Rot R0(1,0,0,0);

/*
    euler->calc(dt, N, S0, R0);
    euler->timing();
    euler->saveTR("TRfirst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    euler->saveTS("TSfirt.dat");

    rodr->calc(dt, N, S0, R0);
    rodr->timing();
    rodr->saveTR("TRrodrConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    rodr->saveTS("TSrodrConst.dat");

    Adm->calc(dt, N, S0, R0);
    Adm->timing();
    Adm->saveTR("TRadamsMConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    Adm->saveTS("TSadamsMConst.dat");

    runge->calc(dt, N, S0, R0);
    runge->timing();
    runge->saveTR("TRrungeConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    runge->saveTS("TSrungeConst.dat");

    emr1->calc(dt, N, S0, R0);
    emr1->timing();
    emr1->saveTR("TRemr1Const.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    emr1->saveTS("TSemr1Const.dat");

    eulerConvent->calc(dt, N, S0, R0);
    eulerConvent->timing();
    eulerConvent->saveTR("TReulerConventConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    eulerConvent->saveTS("TSeulerConventConst.dat");

    eulerConvent->calc(dt, N, S0, R0);
    eulerConvent->timing();
    eulerConvent->saveTR("TReulerConventConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    eulerConvent->saveTS("TSeulerConventConst.dat");

    Ad->calc(dt, N, S0, R0);
    Ad->timing();
    Ad->saveTR("TRdamsConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    Ad->saveTS("TSadamsConst.dat");

    ex->calc(dt, N, S0, R0);
    ex->timing();
    ex->saveTR("TRexactConst.dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    ex->saveTS("TSexactConst.dat");
*/


    //CALCULATE ERROR WITH ERROR MEASUREMENT
    double dtmin = 1E-15;
    double dtmax = 10;
    double logarithm_base = 1.1; //increment od dt for next step, i.e. dt *= logarithm_base

    euler->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    euler->saveError("Error2eulerConst.dat");

    emr1->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr1->saveError("Error2euler1Const.dat");

    emr2->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr2->saveError("Error2euler2Const.dat");


    emr3->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr3->saveError("Error2euler3Const.dat");


    runge->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    runge->saveError("Error2rungeConst.dat");

    rodr->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    rodr->saveError("Error2rodrConst.dat");

    Adm->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    Adm->saveError("Error2adamsMConst.dat");

    eulerConvent->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    eulerConvent->saveError("Error2eulerConventConst.dat");

    eulerConventRescaling->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    eulerConventRescaling->saveError("Error2eulerConventRescalingConst.dat");

    milne->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    milne->saveError("Error2milneConst.dat");

    milnecorr->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    milnecorr->saveError("Error2milneCorrectedConst.dat");

    Ad->cE(dtmin, dtmax, N, logarithm_base, S0, R0);
    Ad->saveError("Error2adamsConst.dat");

    ex->calcError2(dtmin, dtmax, N, logarithm_base, S0, R0);
    ex->saveError("Error2exactConst.dat");



    delete euler;
    delete emr1;
    delete emr2;
    delete emr3;
    delete runge;
    delete rodr;
    delete Adm;
    delete eulerConvent;
    delete eulerConventRescaling;
    delete milne;
    delete milnecorr;
    delete Ad;
    delete ex;

    return 0;
}