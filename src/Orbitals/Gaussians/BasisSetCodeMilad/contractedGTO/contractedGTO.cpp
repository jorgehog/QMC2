#include "contractedGTO.h"

ContractedGTO::ContractedGTO()
{
}

void ContractedGTO::addPrimitive(PrimitiveGTO primitiveGTO)
{
    m_primitivesGTOs.push_back(primitiveGTO);
}

const PrimitiveGTO &ContractedGTO::getPrimitive(const int p) const
{
    return m_primitivesGTOs.at(p);
}

int ContractedGTO::getNumPrimitives() const
{
    return m_primitivesGTOs.size();
}
