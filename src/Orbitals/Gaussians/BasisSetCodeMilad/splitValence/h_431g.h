#ifndef H_431G_H
#define H_431G_H

#include "splitvalence.h"

class H_431G : public SplitValence
{
public:
    H_431G();
    int getAngularMomentum() const;

private:
    ContractedGTO m_contractedGTO_1;
    ContractedGTO m_contractedGTO_2;
};

#endif // H_431G_H
