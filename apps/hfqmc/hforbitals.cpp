#include "hforbitals.h"


using namespace QMC2;

HFOrbitals::HFOrbitals(const uint n_p, vector<hf::ContractedGTO*> contractedGTOs) :
    Orbitals(n_p, 3),
    m_contractedGTOs(contractedGTOs)
{

}

double HFOrbitals::phi(const Walker *walker, int particle, int q_num)
{
    const rowvec & ri = walker->r.row(particle);

    double  Gab = 0.0;

    const hf::ContractedGTO &contracted = m_contractedGTOs.at(q_num);

    const rowvec &corePosition = contracted.center();

    double X = ri(0) - corePosition(0);
    double Y = ri(1) - corePosition(1);
    double Z = ri(2) - corePosition(2);

    double R2 = X * X + Y * Y + Z * Z;

    for(const hf::PrimitiveGTO &primitive : contracted.getPrimitives())
    {
        Gab +=  primitive.weight()
                * std::pow(X, primitive.xPower()) //using std::X to avoid calling arma::X
                * std::pow(Y, primitive.yPower())
                * std::pow(Z, primitive.zPower())
                * std::exp(-primitive.exponent()*R2);
    }

    return Gab;

}

double HFOrbitals::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    return num_diff(walker, particle, q_num, d);
}

double HFOrbitals::lapl_phi(const Walker *walker, int particle, int q_num)
{
    return num_ddiff(walker, particle, q_num);
}

