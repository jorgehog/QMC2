#ifndef O_431G_H
#define O_431G_H

#include "splitvalence.h"

class O_431G : public SplitValence
{
public:
    O_431G();

    int getAngularMomentum() const;

private:
    ContractedGTO m_contractedGTO_1;
    ContractedGTO m_contractedGTO_2;
    ContractedGTO m_contractedGTO_3;
    ContractedGTO m_contractedGTO_4;
    ContractedGTO m_contractedGTO_5;
    ContractedGTO m_contractedGTO_6;
    ContractedGTO m_contractedGTO_7;
    ContractedGTO m_contractedGTO_8;
    ContractedGTO m_contractedGTO_9;
};

#endif // O_431G_H
