#include "hforbitals.h"


using namespace QMC2;

HFOrbitals::HFOrbitals(const uint n_p, vector<const hf::ContractedGTO *> contractedGTOs, const mat& coeffs) :
    Orbitals(n_p, 3),
    m_contractedGTOs(contractedGTOs),
    m_coeffs(coeffs)
{

}

double HFOrbitals::phi(const Walker *walker, int particle, int q_num)
{

    double value = 0;
    const rowvec & ri = walker->r.row(particle);

    for(int i = 0; i < signed(m_coeffs.n_rows); i++){
        value += m_coeffs(i, q_num) * evalContracted(i ,ri);
    }

    return value;
}


double HFOrbitals::del_phi(const Walker *walker, int particle, int q_num, int d)
{
    double value = 0;
    const rowvec & ri = walker->r.row(particle);

    for(int i = 0; i < signed(m_coeffs.n_rows); i++){
        value += m_coeffs(i, q_num) * evaldelContracted(i, ri, d);
    }

    return value;
//    return num_diff(walker, particle, q_num, d);
}


double HFOrbitals::lapl_phi(const Walker *walker, int particle, int q_num)
{

    double value = 0;
    const rowvec & ri = walker->r.row(particle);

    for(int i = 0; i < signed(m_coeffs.n_rows); i++){
        value += m_coeffs(i, q_num) * evallaplContracted(i ,ri);
    }


//    cout << value << "       " << num_ddiff(walker, particle, q_num) << endl;
//    sleep(5);

    return value;
//    return num_ddiff(walker, particle, q_num);
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

    for(const hf::PrimitiveGTO &primitive : contracted.primitiveGTOs())
    {
        Gab +=  primitive.weight()
                * std::pow(X, primitive.xPower()) //using std::X to avoid calling arma::X
                * std::pow(Y, primitive.yPower())
                * std::pow(Z, primitive.zPower())
                * std::exp(-primitive.exponent()*R2);
    }

    return Gab;
}


double HFOrbitals::evallaplContracted(uint k, const rowvec & r)
{
    double  ddGab = 0.0;

    const hf::ContractedGTO &contracted = *m_contractedGTOs.at(k);

    const rowvec &corePosition = contracted.center();

    double X = r(0) - corePosition(0);
    double Y = r(1) - corePosition(1);
    double Z = r(2) - corePosition(2);

    double R2 = X * X + Y * Y + Z * Z;

    for(const hf::PrimitiveGTO &primitive : contracted.primitiveGTOs())
    {

        double G =  primitive.weight()
                * std::pow(X, primitive.xPower()) //using std::X to avoid calling arma::X
                * std::pow(Y, primitive.yPower())
                * std::pow(Z, primitive.zPower())
                * std::exp(-primitive.exponent()*R2);


        double Gp2 = G * 4.0 * primitive.exponent() * primitive.exponent() * R2;

        double Gm2 = G * (
                      primitive.xPower() * (primitive.xPower() - 1.) / (X*X)
                    + primitive.yPower() * (primitive.yPower() - 1.) / (Y*Y)
                    + primitive.zPower() * (primitive.zPower() - 1.) / (Z*Z)
                    );

        ddGab += -G * 2.0 * primitive.exponent()
                *  ( (1 + 2.0 * primitive.xPower() )
                   + (1 + 2.0 * primitive.yPower() )
                   + (1 + 2.0 * primitive.zPower() ))
                + Gp2 + Gm2;

    }

    return ddGab;
}


double HFOrbitals::evaldelContracted(uint k, const rowvec & r, int d)
{
    double  dGab = 0.0;

    const hf::ContractedGTO &contracted = *m_contractedGTOs.at(k);

    const rowvec &corePosition = contracted.center();

    double X = r(0) - corePosition(0);
    double Y = r(1) - corePosition(1);
    double Z = r(2) - corePosition(2);

    double R2 = X * X + Y * Y + Z * Z;

    for(const hf::PrimitiveGTO &primitive : contracted.primitiveGTOs())
    {
        double expTerm = std::exp(-primitive.exponent()*R2);

        double G =  primitive.weight()
                * std::pow(X, primitive.xPower()) //using std::X to avoid calling arma::X
                * std::pow(Y, primitive.yPower())
                * std::pow(Z, primitive.zPower())
                *  expTerm;


        double Gp = G * (r(d) - corePosition(d));
        double Gm = G / (r(d) - corePosition(d));
        dGab += primitive.powers()[d] * Gm  - 2.0 * primitive.exponent() * Gp;

    }

    return dGab;
}



