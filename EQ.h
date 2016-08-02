//
// Created by c8888 on 28.07.16.
//

#ifndef GEOMETRIC_ALGEBRA_CLASSES_EQ_H
#define GEOMETRIC_ALGEBRA_CLASSES_EQ_H
#include <iostream>
#include <include/vsr/space/vsr_ega3D_types.h>
#include "vsr.h"
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "data.h"
#include "Afun.h"


class data;
class Bfun;

class EQ {
public:
    //************ CONSTRUCTORS *******************
    EQ(Afun *Af) : A(Af) { }// in a constructor one must specify what is the function A(t)
    //*********************************************

    virtual void calc(double dt, int N,  vsr::ega::Vec S0T, vsr::ega::Rot R0T) = 0; //erase data, fill data with new solution, change time
    void calcError1(double dtmin, double dtmax, int N, double logarithm_base, vsr::ega::Vec S0T, vsr::ega::Rot R0T); //max(abs(1-|s|)), adds pairs of dt and error measurement into data d
    void calcError2(double dtmin, double dtmax, int N, double logarithm_base, vsr::ega::Vec S0T, vsr::ega::Rot R0T); //acos(angle between S and S_exact)
    void calcError3(double dtmin, double dtmax, int N, double logarithm_base, vsr::ega::Vec S0T, vsr::ega::Rot R0T); //max(abs(1-R(~R))), adds pairs of dt and error measurement into data d
    void calcError4(double dtmin, double dtmax, int N, double logarithm_base, vsr::ega::Vec S0T, vsr::ega::Rot R0T); //min(cos(angle between S and S_exact)), adds pairs of dt and error measurement into data d
    void calcError5(double dtmin, double dtmax, int N, double logarithm_base, vsr::ega::Vec S0T, vsr::ega::Rot R0T); //acos(angle between S and S_exact) without dividing by R(~R)
    void timing(); //returns last execution time and N
    void push_back(std::pair< double, std::pair<vsr::ega::Vec, vsr::ega::Rot> > ); //add t,S,R with this
    void push_back(std::pair<double, double>); //add error with this
    void erase(); //erases data
    void set_time(double t);
    void set_N(int N);

    //saving to file
    void saveTS(const char* filename);
    void saveTR(const char* filename);
    void saveError(const char* filename);

    Afun *A;
private:
    data d;
    double time; //store last time of N iterations execution
    int N; //store last number of iterations
};


#endif //GEOMETRIC_ALGEBRA_CLASSES_EQ_H
