//
// Created by c8888 on 28.07.16.
//

#ifndef GEOMETRIC_ALGEBRA_CLASSES_DATA_H
#define GEOMETRIC_ALGEBRA_CLASSES_DATA_H


#include <include/vsr/space/vsr_ega3D_types.h>
#include <vector>


class Sdata {
public:
    vsr::ega::Vec at(size_t i); //return ith element
    void push_back(vsr::ega::Vec);
    size_t size();
    void erase();
private:
    std::vector<vsr::ega::Vec> V;
};
class Rdata {
public:
    vsr::ega::Rot at(size_t i); //return ith element
    void push_back(vsr::ega::Rot);
    size_t size();
    void erase();
private:
    std::vector<vsr::ega::Rot> R;
};

class Tdata {
public:
    double at(size_t i); //return ith element
    void push_back(double t);
    size_t size();

    void erase();
private:
    std::vector<double> T;
};

class Edata {
public:
    std::pair<double, double> at(size_t i); //return ith pair of dt and error
    void push_back(std::pair<double, double>);
    size_t size();
    void erase();

private:
    std::vector<std::pair<double, double> > E;
};

class data {
public:
    data() {}
    void erase(); //erase all data except for E
    void eraseE(); //erase E container
    void push_back(std::pair< double, std::pair<vsr::ega::Vec, vsr::ega::Rot> >); //need to create pair of time and pair of S and R to push it back to data
    void push_back(std::pair<double, double> errorAtdt);
    vsr::ega::Vec getS(size_t i);
    vsr::ega::Rot getR(size_t i);
    double getT(size_t i);
    std::pair<double, double> getE(size_t i);

    size_t size();
    void saveTS(const char* filename); //saves (t, S(t)) fo file "filename"
    void saveTR(const char* filename); //the same for Rotor
    void saveError(const char* filename); //and for Error
private:
    Sdata S;
    Rdata R;
    Tdata T;
    Edata E; //stores pairs of dt and error
};




#endif //GEOMETRIC_ALGEBRA_CLASSES_DATA_H
