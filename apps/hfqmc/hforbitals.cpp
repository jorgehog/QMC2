#include "hforbitals.h"


using namespace QMC2;

HFOrbitals::HFOrbitals(const uint n_p, vector<const hf::ContractedGTO *> contractedGTOs) :
    Orbitals(n_p, 3),
    m_contractedGTOs(contractedGTOs)
{

}

double HFOrbitals::evalContracted(uint k, const rowvec & r)
{
    double  Gab = 0.0;

    const hf::ContractedGTO &contracted = *m_contractedGTOs.at(k);

    const rowvec &corePosition = contracted.center();

    double X = r(0) - corePosition(0);
    double Y = r(1) - corePosition(1);
    double Z = r(2) - corePosition(2);

    double R2 = X * X + Y * Y + Z * Z;

    for(const hf::PrimitiveGTO &primitive : contracted.primitivesGTOs())
    {
        Gab +=  primitive.weight()
                * std::pow(X, primitive.xPower()) //using std::X to avoid calling arma::X
                * std::pow(Y, primitive.yPower())
                * std::pow(Z, primitive.zPower())
                * std::exp(-primitive.exponent()*R2);
    }

    return Gab;
}

double HFOrbitals::phi(const Walker *walker, int particle, int q_num)
{

    const rowvec & ri = walker->r.row(particle);



    return 0.3285*evalContracted(0, ri) +
           0.2694*evalContracted(1, ri) +
           0.3285*evalContracted(2, ri) +
           0.2694*evalContracted(3, ri);
//    return 0.3285*evalContracted(0, ri) + 0.1207*evalContracted(1, ri) + 0.7620*evalContracted(2, ri) + 1.1308*evalContracted(3, ri);

}

double HFOrbitals::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    return num_diff(walker, particle, q_num, d);
}

double HFOrbitals::lapl_phi(const Walker *walker, int particle, int q_num)
{
    return num_ddiff(walker, particle, q_num);
}

