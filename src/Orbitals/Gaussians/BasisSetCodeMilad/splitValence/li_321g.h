#ifndef LI_321G_H
#define LI_321G_H


#include "splitvalence.h"

class Li_321G : public SplitValence
{
public:
    Li_321G();
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

#endif // LI_321G_H
