#include "hforbitals.h"


using namespace QMC2;

HFOrbitals::HFOrbitals(hf::System *system) :
    Orbitals(n_p = system->getNumOfElectrons(), 3),
    m_system(system)
{

}

double HFOrbitals::phi(const Walker *walker, int particle, int q_num)
{
    const rowvec & ri = walker->r.row(particle);

    return 0;
}

double HFOrbitals::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    return num_diff(walker, particle, q_num, d);
}

double HFOrbitals::lapl_phi(const Walker *walker, int particle, int q_num)
{
    return num_ddiff(walker, particle, q_num);
}


//double  Gab = 0.0;

//const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
//const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));

//const rowvec &corePositionA = coreA->corePosition();

//double Xa = x - corePositionA(0);
//double Ya = y - corePositionA(1);
//double Za = z - corePositionA(2);

//double Ra = Xa * Xa + Ya * Ya + Za * Za;

//for(int i = 0; i < contractedA.getNumPrimitives(); i++){
//    const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);
//    Gab +=  primitiveA.weight()
//            * pow(Xa, primitiveA.xPower())
//            * pow(Ya, primitiveA.yPower())
//            * pow(Za, primitiveA.zPower())
//            * exp(-primitiveA.exponent()*Ra);

//}

//return Gab;
