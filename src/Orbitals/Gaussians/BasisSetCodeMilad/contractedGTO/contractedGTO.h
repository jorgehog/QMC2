#ifndef CONTRACTEDGTO_H
#define CONTRACTEDGTO_H

#include <iostream>
#include <armadillo>
#include "../primitiveGTO/primitiveGTO.h"

using namespace arma;
using namespace std;


class ContractedGTO
{
public:
    ContractedGTO();

    void addPrimitive(PrimitiveGTO primitiveGTO);
    int getNumPrimitives() const;
    const PrimitiveGTO &getPrimitive(const int p) const;

private:
    vector<PrimitiveGTO> m_primitivesGTOs;
};


#endif // CONTRACTEDGTO_H
