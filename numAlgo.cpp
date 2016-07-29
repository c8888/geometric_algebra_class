//
// Created by c8888 on 28.07.16.
//

#include "numAlgo.h"
#include <ctime>

//PLEASE IMPLEMENT NUMERICAL ALGORITHMS HERE. DO NOT FORGET TO ADD DECLARATIONS TO numAlgo.h FILE!!!
//EVERY IMPLEMENTATION SHOULD MEASURE EXECUTION TIME ON ITS OWN. AT THE END EXECUTE set_time(double t) and set_N(N) TO
//TELL THE PROGRAM HOW FAST THE IMPLEMENTATION IS
//

using namespace std;
using namespace vsr::ega;

void EulerMethodRotor::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Rot R = R0T;

    for(int i = 0; i<N; i++){
        R += (A->value(i*dt)) * R * dt; //A(T) IS THE FUNCTION RETURNING BIVECTOR SO NO NEED TO MULTIPLY IT BY PSEUDOSCALAR. YOU CAN ACCESS IT HERE.
        s = R * S0T * (~R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }

    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}

void RungeKutta4thMethodRotor::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Rot R = R0T;

    Rot k1, k2, k3, k4;

    for(int i = 0; i<N; i++){
        k1 = f(A, i*dt, R) * dt;
        k2 = f(A, i*dt + dt/2, R + k1*0.5) * dt;
        k3 = f(A, i*dt + dt/2, R + k2*0.5) * dt;
        k4 = f(A, i*dt + dt, R + k3) * dt;
        R += (k1 + k2*2 + k3*2 + k4)/6;
        s = R * S0T * (~R);

        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //ADD NEXT POINT
    }

    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE

}

void EulerMethodRotorRescaling1::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Rot R = R0T;

    for(int i = 0; i<N; i++){
        R += (A->value(i*dt)) * R * dt;
        R = R*Sca(1/((R*(~R))[0]));
        s = R * S0T * (~R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }

    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}
void EulerMethodRotorRescaling2::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Rot R = R0T;

    for(int i = 0; i<N; i++){
        R += (A->value(i*dt)) * R * dt;
        s = R * S0T * (!R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }

    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}
void EulerMethodRotorRescaling3::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Rot R = R0T;

    for(int i = 0; i<N; i++){
        R += (A->value(i*dt)) * R * dt;
        R = R*Sca(1/((R*(~R))[0]));
        s = R * S0T * (!R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }

    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}