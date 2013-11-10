#ifndef BASISSET_H
#define BASISSET_H

#include <iostream>
#include <armadillo>

#include "contractedGTO/contractedGTO.h"

using namespace arma;
using namespace std;

class BasisSet
{
public:
    BasisSet();

    rowvec corePosition() const;
    void setCorePosition(const rowvec &corePosition);
    const ContractedGTO &getContracted(const int c) const;
    int getNumContracted() const;

    virtual int getAngularMomentum() const = 0;

protected:
    vector<ContractedGTO> m_contractedGTOs;
    rowvec m_corePosition;
    int m_angularMomentum;

};

#endif // BASISSET_H
