//
// Created by c8888 on 28.07.16.
//

#include <fstream>
#include "data.h"

using namespace std;
using namespace vsr::ega;

void Sdata::erase(){
    V.erase(V.begin(), V.end());
}

void Rdata::erase(){
    R.erase(R.begin(), R.end());
}
void Tdata::erase(){
    T.erase(T.begin(), T.end());
}
void Edata::erase(){
    E.erase(E.begin(), E.end());
}

void data::erase(){
    S.erase();
    R.erase();
    T.erase();
}
void data::eraseE(){
    E.erase();
}

void Sdata::push_back(Vec SS) {
    V.push_back(SS);
}

void Rdata::push_back(Rot RR) {
    R.push_back(RR);
}
void Tdata::push_back(double TT) {
    T.push_back(TT);
}

void data::push_back(pair<double, pair<Vec,Rot> > p) {
    T.push_back(p.first);
    S.push_back(p.second.first);
    R.push_back(p.second.second);
}

Vec Sdata::at(size_t i) {
    return V.at(i);
}

Rot Rdata::at(size_t i) {
    return R.at(i);
}

double Tdata::at(size_t i) {
    return T.at(i);
}
size_t Sdata::size() {
    return V.size();
}

size_t Rdata::size() {
    return R.size();
}

size_t Tdata::size() {
    return T.size();
}

size_t data::size() {
    if(S.size() == R.size() && R.size() == T.size()) return S.size();
    else __throw_length_error("S is not the same size as R or T! Terminating.");
}

void data::saveTS(const char* filename) {
    ofstream of;
    of.open(filename);
    if(!of) __throw_ios_failure("Could not create the file. Terminating.");
    for(size_t i = 0; i<S.size(); i++) {
        of << T.at(i) << '\t' << S.at(i) << endl; //operator << can be applied to Vec. It is overloaded in the library.
    }
}
void data::saveTR(const char* filename) {
    ofstream of;
    of.open(filename);
    if(!of) __throw_ios_failure("Could not create the file. Terminating.");
    for(size_t i = 0; i<R.size(); i++) {
        of << T.at(i) << '\t' << R.at(i) << endl; //operator << can be applied to Rot. It is overloaded in the library.
    }
}

pair<double, double> Edata::at(size_t i){
    return E.at(i);
}
void Edata::push_back(pair<double, double> p){
    E.push_back(p);
}
size_t Edata::size(){
    return E.size();
}

void data::saveError(const char *filename) {
    ofstream of;
    of.open(filename);
    if(!of) __throw_ios_failure("Could not create the file. Terminating.");
    for(size_t i = 0; i<E.size(); i++) {
        of << E.at(i).first << '\t' << E.at(i).second << endl; //operator << can be applied to Rot. It is overloaded in the library.
    }
}

void data::push_back(pair<double , double> errorAtdt){
    E.push_back(errorAtdt);
}

Vec data::getS(size_t i) {
    return S.at(i);
}
Rot data::getR(size_t i) {
    return R.at(i);
}
double data::getT(size_t i) {
    return T.at(i);
}
pair<double, double> data::getE(size_t i) {
    return E.at(i);
};
