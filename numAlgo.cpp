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


Rot power(Biv B, int j) {Rot Br=B; for (int i = 1; i<=j; i++) {Br = B*Br;} return Br;}//DEFINING INT POWER OF A BIV WHICH GIVES EITHER BIV OR ROTOR


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

void RodriguesFormula::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Rot R = R0T;

    Biv A1, A2, A3, a1, a2, a3;
    Rot c1, c2;
    Biv Omega;
    double t;



    for(int i=1; i<=N; i++) {
        t = i*dt;
        A1 = A->value(t + (0.5 - sqrt(15)/10.)*dt);
        A2 = A->value(t + dt/2.);
        A3 = A->value(t + (0.5 + sqrt(15)/10.)*dt);
        a1 = A2*dt;
        a2 = Sca(sqrt(15)/3.)*dt*(A3 - A1);
        a3 = Sca(10./3.)*dt*(A3 - Sca(2)*A2 + A1);
        c1 = a1*a2-a2*a1;
        c2 = Sca(-1./60.)*(a1*(Sca(2)*a3+c1) - (Sca(2)*a3+c1)*a1);
        Omega = a1 + a3*Sca(1./12.) + ((Sca(-20)*a1-a3+c1)*(a2+c2)-(a2+c2)*(Sca(-20)*a1-a3+c1))*Sca(1./240.);
        R = (Sca(cos(Omega.norm())) + Omega*Sca(sin(Omega.norm())/Omega.norm()))*R;
        s = R*S0T*(~R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }



    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}

void AdamsMulton::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Pss I =Pss(1);

    Rot R, Rk, Rk_1, Rk_2, Rk_3;

    Rk_3 = R0T;
    Rk_2 = Rk_3 + Sca(dt)*power(A->value(1*dt)*I,1)*Rk_3 + Sca(dt*dt/2.)*power(A->value(1*dt)*I,2)*Rk_3 + Sca(dt*dt*dt/6.)*power(A->value(1*dt)*I,3)*Rk_3 + Sca(dt*dt*dt*dt/24.)*power(A->value(1*dt)*I,4)*Rk_3 + Sca(dt*dt*dt*dt*dt/120.)*power(A->value(1*dt)*I,5)*Rk_3 + Sca(dt*dt*dt*dt*dt*dt/720.)*power(A->value(1*dt)*I,6)*Rk_3 + Sca(dt*dt*dt*dt*dt*dt*dt/5040.)*power(A->value(1*dt)*I,7)*Rk_3;
    Rk_1 = Rk_2 + Sca(dt)*power(A->value(2*dt)*I,1)*Rk_2 + Sca(dt*dt/2.)*power(A->value(2*dt)*I,2)*Rk_2 + Sca(dt*dt*dt/6.)*power(A->value(2*dt)*I,3)*Rk_2 + Sca(dt*dt*dt*dt/24.)*power(A->value(2*dt)*I,4)*Rk_2 + Sca(dt*dt*dt*dt*dt/120.)*power(A->value(2*dt)*I,5)*Rk_2 + Sca(dt*dt*dt*dt*dt*dt/720.)*power(A->value(2*dt)*I,6)*Rk_2 + Sca(dt*dt*dt*dt*dt*dt*dt/5040.)*power(A->value(2*dt)*I,7)*Rk_2;
    Rk = Rk_1 + Sca(dt)*power(A->value(3*dt)*I,1)*Rk_1 + Sca(dt*dt/2.)*power(A->value(3*dt)*I,2)*Rk_1 + Sca(dt*dt*dt/6.)*power(A->value(3*dt)*I,3)*Rk_1 + Sca(dt*dt*dt*dt/24.)*power(A->value(3*dt)*I,4)*Rk_1 + Sca(dt*dt*dt*dt*dt/120.)*power(A->value(3*dt)*I,5)*Rk_1 + Sca(dt*dt*dt*dt*dt*dt/720.)*power(A->value(3*dt)*I,6)*Rk_1 + Sca(dt*dt*dt*dt*dt*dt*dt/5040.)*power(A->value(3*dt)*I,7)*Rk_1;

    s = Rk_3*S0T*(~Rk_3);
    push_back(make_pair< double, pair<Vec, Rot> >(0*dt, pair<Vec, Rot>(s, Rk_3))); //EXAMPLE OF HOW TO ADD NEXT POINT

    s = Rk_2*S0T*(~Rk_2);
    push_back(make_pair< double, pair<Vec, Rot> >(1*dt, pair<Vec, Rot>(s, Rk_2))); //EXAMPLE OF HOW TO ADD NEXT POINT

    s = Rk_1*S0T*(~Rk_1);
    push_back(make_pair< double, pair<Vec, Rot> >(2*dt, pair<Vec, Rot>(s, Rk_1))); //EXAMPLE OF HOW TO ADD NEXT POINT

    s = Rk*S0T*(~Rk);
    push_back(make_pair< double, pair<Vec, Rot> >(3*dt, pair<Vec, Rot>(s, Rk))); //EXAMPLE OF HOW TO ADD NEXT POINT


    for(int i=4; i<N; i++) {
        R = Rk + Sca(dt/24.)*(Sca(55.)*A->value((i-1)*dt)*Rk - Sca(59.)*A->value((i-2)*dt)*Rk_1 + Sca(37.)*A->value((i-3)*dt)*Rk_2 - Sca(9.)*A->value((i-4)*dt)*Rk_3);
        R = Rk + Sca(dt/24.)*(Sca(9.)*A->value((i)*dt)*R + Sca(19.)*A->value((i-1)*dt)*Rk - Sca(5.)*A->value((i-2)*dt)*Rk_1 + A->value((i-3)*dt)*Rk_2);

        Rk_3 = Rk_2;
        Rk_2 = Rk_1;
        Rk_1 = Rk;

        s = R*S0T*(~R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }



    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}

void EulerMethodConvent::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;

    for(int i = 0; i<N; i++){
        Vec Bfield(A->value(i*dt)[2], -A->value(i*dt)[1],  A->value(i*dt)[0]);
        s += (s^Bfield).dual() * dt * 2;
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R0T)));
    }


    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}

void EulerMethodConventRescaling::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;

    for(int i = 0; i<N; i++){
        Vec Bfield(A->value(i*dt)[2], -A->value(i*dt)[1],  A->value(i*dt)[0]);
        s += (s^Bfield).dual() * dt * 2;
        s = s*(1/s.norm());
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R0T)));
    }

    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}

void EulerHeunMethodRotor::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Rot R = R0T;
    Rot R0;

    for(int i = 1; i<=N; i++){
        R0=R;
        R += (A->value(i*dt)) * R * dt; //A(T) IS THE FUNCTION RETURNING BIVECTOR SO NO NEED TO MULTIPLY IT BY PSEUDOSCALAR. YOU CAN ACCESS IT HERE.
        R = R0 + Sca(dt/2.)*(A->value(i*dt)*R0+A->value((i+1)*dt)*R);
        s = R * S0T * (~R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }

    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}

void Adams::calc(double dt, int N, vsr::ega::Vec S0T, vsr::ega::Rot R0T) {
    std::clock_t c_start = std::clock(); //ADD THIS LINE TO EVERY FUNCTION HERE
    erase(); //THIS LINE IS OBLIGATORY FOR EVERY METHOD!!

    Vec s = S0T;
    Pss I =Pss(1);

    Rot R, Rk, Rk_1, Rk_2, Rk_3;

    Rk_3 = R0T;
    Rk_2 = Rk_3 + Sca(dt)*power(A->value(1*dt)*I,1)*Rk_3 + Sca(dt*dt/2.)*power(A->value(1*dt)*I,2)*Rk_3 + Sca(dt*dt*dt/6.)*power(A->value(1*dt)*I,3)*Rk_3 + Sca(dt*dt*dt*dt/24.)*power(A->value(1*dt)*I,4)*Rk_3 + Sca(dt*dt*dt*dt*dt/120.)*power(A->value(1*dt)*I,5)*Rk_3 + Sca(dt*dt*dt*dt*dt*dt/720.)*power(A->value(1*dt)*I,6)*Rk_3 + Sca(dt*dt*dt*dt*dt*dt*dt/5040.)*power(A->value(1*dt)*I,7)*Rk_3;
    Rk_1 = Rk_2 + Sca(dt)*power(A->value(2*dt)*I,1)*Rk_2 + Sca(dt*dt/2.)*power(A->value(2*dt)*I,2)*Rk_2 + Sca(dt*dt*dt/6.)*power(A->value(2*dt)*I,3)*Rk_2 + Sca(dt*dt*dt*dt/24.)*power(A->value(2*dt)*I,4)*Rk_2 + Sca(dt*dt*dt*dt*dt/120.)*power(A->value(2*dt)*I,5)*Rk_2 + Sca(dt*dt*dt*dt*dt*dt/720.)*power(A->value(2*dt)*I,6)*Rk_2 + Sca(dt*dt*dt*dt*dt*dt*dt/5040.)*power(A->value(2*dt)*I,7)*Rk_2;
    Rk = Rk_1 + Sca(dt)*power(A->value(3*dt)*I,1)*Rk_1 + Sca(dt*dt/2.)*power(A->value(3*dt)*I,2)*Rk_1 + Sca(dt*dt*dt/6.)*power(A->value(3*dt)*I,3)*Rk_1 + Sca(dt*dt*dt*dt/24.)*power(A->value(3*dt)*I,4)*Rk_1 + Sca(dt*dt*dt*dt*dt/120.)*power(A->value(3*dt)*I,5)*Rk_1 + Sca(dt*dt*dt*dt*dt*dt/720.)*power(A->value(3*dt)*I,6)*Rk_1 + Sca(dt*dt*dt*dt*dt*dt*dt/5040.)*power(A->value(3*dt)*I,7)*Rk_1;

    s = Rk_3*S0T*(~Rk_3);
    push_back(make_pair< double, pair<Vec, Rot> >(0*dt, pair<Vec, Rot>(s, Rk_3))); //EXAMPLE OF HOW TO ADD NEXT POINT

    s = Rk_2*S0T*(~Rk_2);
    push_back(make_pair< double, pair<Vec, Rot> >(1*dt, pair<Vec, Rot>(s, Rk_2))); //EXAMPLE OF HOW TO ADD NEXT POINT

    s = Rk_1*S0T*(~Rk_1);
    push_back(make_pair< double, pair<Vec, Rot> >(2*dt, pair<Vec, Rot>(s, Rk_1))); //EXAMPLE OF HOW TO ADD NEXT POINT

    s = Rk*S0T*(~Rk);
    push_back(make_pair< double, pair<Vec, Rot> >(3*dt, pair<Vec, Rot>(s, Rk))); //EXAMPLE OF HOW TO ADD NEXT POINT


    for(int i=4; i<N; i++) {
        R = Rk + Sca(dt/24.)*(Sca(55.)*A->value((i-1)*dt)*Rk - Sca(59.)*A->value((i-2)*dt)*Rk_1 + Sca(37.)*A->value((i-3)*dt)*Rk_2 - Sca(9.)*A->value((i-4)*dt)*Rk_3);

        Rk_3 = Rk_2;
        Rk_2 = Rk_1;
        Rk_1 = Rk;

        s = R*S0T*(~R);
        push_back(make_pair< double, pair<Vec, Rot> >(i*dt, pair<Vec, Rot>(s, R))); //EXAMPLE OF HOW TO ADD NEXT POINT
    }



    std::clock_t c_end = std::clock();
    set_time(1000.0 * (c_end-c_start) / CLOCKS_PER_SEC); //set time of executions in miliseconds
    set_N(N); //ADD THESE THREE LINES TO EVERY FUNCTION HERE
}