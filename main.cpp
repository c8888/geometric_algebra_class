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
    //CALCULATE THE ROTOR AND SPIN IN TIME

    Biv B(1,1,1);
    cout<<B.norm()<<endl;
    Rot RRRRR;
    RRRRR = Rot(1,2,3,4);
    Rot YYYY(1,2,1,1);

    int N = 40000;
    double dt=1E-12;
    Vec S0(0,0,1);
    Rot R0(1,0,0,0);

    emr->calc(dt, N, S0, R0);
    emr->timing();
    emr->saveTR(".dat"); //export to file. Filename should be known before compliation. It by far the simplest to type it before with keyboard. Overwrites th file.
    emr->saveTS("TSfirt.dat");


    //CALCULATE ERROR WITH ERROR MEASUREMENT METHOD 1 (change in modulus of spin)
    double dtmin = 1E-22;
    double dtmax = 1E-8;
    double logarithm_base = 1.1; //increment od dt for next step, i.e. dt *= logarithm_base
    emr->calcError1(dtmin, dtmax, N, logarithm_base, S0, R0);
    emr->saveError("Error1Euler.dat");






    return 0;
}