//
// Created by c8888 on 28.07.16.
//

#ifndef GEOMETRIC_ALGEBRA_CLASSES_NUMALGO_H
#define GEOMETRIC_ALGEBRA_CLASSES_NUMALGO_H

#include "EQ.h"
#include "Afun.h"


class EulerMethodRotor: public EQ {
public:
    EulerMethodRotor(Afun* Af): EQ(Af){} //always include this line
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T); //erase data, fill data with new solution, change time
    //these two lines followed by definition of void calc() should be in every newly implemented method
};
class RungeKutta4thMethodRotor: public EQ {
public:
    RungeKutta4thMethodRotor(Afun* Af): EQ(Af){}
    auto f(Afun* Af, double t, const vsr::ega::Rot& R){
        return Af->value(t)*R;
    } //simplifies the coding
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);

};

class EulerMethodRotorRescaling1: public EQ {
public:
    EulerMethodRotorRescaling1(Afun* Af): EQ(Af){}
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);

};

class EulerMethodRotorRescaling2: public EQ {
public:
    EulerMethodRotorRescaling2(Afun* Af): EQ(Af){}
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);

};

class EulerMethodRotorRescaling3: public EQ {
public:
    EulerMethodRotorRescaling3(Afun* Af): EQ(Af){}
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

class RodriguesFormula: public EQ {
public:
    RodriguesFormula(Afun* Af): EQ(Af){}
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

class AdamsMulton: public EQ {
public:
    AdamsMulton(Afun* Af): EQ(Af){}
    auto f(Afun* Af, double t, const vsr::ega::Rot& R){
        return Af->value(t)*R;
    } //simplifies the coding
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

class EulerMethodConvent: public EQ {
public:
    EulerMethodConvent(Afun* Af): EQ(Af){}
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

class EulerMethodConventRescaling: public EQ {
public:
    EulerMethodConventRescaling(Afun* Af): EQ(Af){}
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

class EulerHeunMethodRotor: public EQ {
public:
    EulerHeunMethodRotor(Afun* Af): EQ(Af){} //always include this line
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T); //erase data, fill data with new solution, change time
    //these two lines followed by definition of void calc() should be in every newly implemented method
};

class Adams: public EQ {
public:
    Adams(Afun* Af): EQ(Af){}
    auto f(Afun* Af, double t, const vsr::ega::Rot& R){
        return Af->value(t)*R;
    } //simplifies the coding
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

class MilneCorrected: public EQ {
public:
    MilneCorrected(Afun* Af): EQ(Af){}
    auto f(Afun* Af, double t, const vsr::ega::Rot& R){
        return Af->value(t)*R;
    } //simplifies the coding
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

class Milne: public EQ {
public:
    Milne(Afun* Af): EQ(Af){}
    auto f(Afun* Af, double t, const vsr::ega::Rot& R){
        return Af->value(t)*R;
    } //simplifies the coding
    void calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T);
};

#endif //GEOMETRIC_ALGEBRA_CLASSES_NUMALGO_H
