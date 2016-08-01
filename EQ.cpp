//
// Created by c8888 on 28.07.16.
//

#include "EQ.h"

using namespace std;
using namespace vsr::ega;


void EQ::erase() {
    d.erase();
}

void EQ::push_back(std::pair< double, std::pair<vsr::ega::Vec, vsr::ega::Rot> > p) {
    d.push_back(p);
}

void EQ::push_back(std::pair<double, double> p) {
    d.push_back(p);
}

void EQ::set_time(double t) {
    time = t;
}

void EQ::saveTS(const char *filename) {
    d.saveTS(filename);
}
void EQ::saveTR(const char *filename) {
    d.saveTR(filename);
}

void EQ::saveError(const char *filename) {
    d.saveError(filename);
}

void EQ::timing() {
    cout<<"Time: "<<time<<endl<<"Number of iterations: "<<N<<endl;
}

void EQ::set_N(int Niter){
    N = Niter;
}

void EQ::calcError1(double dtmin, double dtmax, int N, double logarithm_base, Vec S0T, Rot R0T){
    if(logarithm_base<=1) __throw_invalid_argument("bad argument logarithm_base. Should be more than 1.");
    erase(); //THIS LINE IS OBLIGATORY!!!!
    double dt = dtmin;

    while(dt<dtmax){
        double diff = 0;
        calc(dt, N, S0T, R0T);
        for(int i = 0; i < d.size(); i++){
            if( abs(1-d.getS(i).norm()) > diff ) diff = abs(1-d.getS(i).norm());
        }
        push_back(pair<double, double> (log10(dt), log10(diff)));
        dt *= logarithm_base;
    }
}

void EQ::calcError2(double dtmin, double dtmax, int N, double logarithm_base, Vec S0T, Rot R0T){
    if(logarithm_base<=1) __throw_invalid_argument("bad argument logarithm_base. Should be more than 1.");
    double dt = dtmin;
    erase(); //THIS LINE IS OBLIGATORY!!!!

    while(dt<dtmax){
        double theta = 0;
        calc(dt, N, S0T, R0T);
        for(unsigned long i=0; i<d.size(); i++){
            Rot Rexact = A->exactSol(d.getT(i), R0T); //getT(i) is the time at which getS(i) is specified
            Vec Sexact = Rexact * S0T * (~Rexact);
            double tmp = acos((d.getS(i)<=Sexact)[0]/(d.getR(i)*(~(d.getR(i))))[0]);
            if(abs(tmp)>abs(theta)) theta = abs(tmp);
        }
        push_back(pair<double, double> (log10(dt), log10(theta)));
        dt *= logarithm_base;
    }
}

void EQ::calcError3(double dtmin, double dtmax, int N, double logarithm_base, Vec S0T, Rot R0T){
    if(logarithm_base<=1) __throw_invalid_argument("bad argument logarithm_base. Should be more than 1.");
    erase(); //THIS LINE IS OBLIGATORY!!!!
    double dt = dtmin;

    while(dt<dtmax){
        double diff = 0;
        calc(dt, N, S0T, R0T);
        for(int i = 0; i < d.size(); i++){
            if( abs(1-(d.getR(i)*(~d.getR(i)))[0]) > diff ) diff = abs(1-(d.getR(i)*(~d.getR(i)))[0]);
        }
        push_back(pair<double, double> (log10(dt), log10(diff)));
        dt *= logarithm_base;
    }
}
void EQ::calcError4(double dtmin, double dtmax, int N, double logarithm_base, Vec S0T, Rot R0T){
    if(logarithm_base<=1) __throw_invalid_argument("bad argument logarithm_base. Should be more than 1.");
    double dt = dtmin;
    erase(); //THIS LINE IS OBLIGATORY!!!!

    while(dt<dtmax){
        double costheta = 1000000000;
        calc(dt, N, S0T, R0T);
        for(unsigned long i=0; i<d.size(); i++){
            Rot Rexact = A->exactSol(d.getT(i), R0T); //getT(i) is the time at which getS(i) is specified
            Vec Sexact = Rexact * S0T * (~Rexact);
            double tmp = abs((d.getS(i)<=Sexact)[0]/(d.getR(i)*(~(d.getR(i))))[0]);
            if(tmp<costheta) costheta = tmp;
        }
        push_back(pair<double, double> (log10(dt), log10(costheta)));
        dt *= logarithm_base;
    }
}