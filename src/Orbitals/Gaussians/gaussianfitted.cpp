#include "gaussianfitted.h"
#include "../../Walker/Walker.h"
#include "../../BasisFunctions/expandedbasisfunctions.h"
#include "../../BasisFunctions/Gaussians/gaussians.h"


GaussianFitted::GaussianFitted(int n_p, int dim) :
    Orbitals(n_p, dim),
    CURRENT(0)
{
    name = "gaussianFitted";
}

void GaussianFitted::getGaussianFromPGTO(const PrimitiveGTO &PGTO)
{
    double alpha = PGTO.exponent();
    arma::rowvec powers = PGTO.powers();

    double* expFactor = NULL;
    for (unsigned int i = 0; i < uniqueAlphas.size(); ++i) {
        if (uniqueAlphas.at(i) == alpha) {
            expFactor = expFactors.at(i);
        }
    }

    if (expFactor == NULL) {
        uniqueAlphas.push_back(alpha);

        expFactor = new double();
        expFactors.push_back(expFactor);
    }

    if (expFactor == NULL) {
        std::cout << "SOMETHING IS WRONG EXP FACTOR SHOULD NEVER BE NULL..." << std::endl;
    }

    locateGaussianFromPowers(powers, alpha, expFactor);

}

void GaussianFitted::locateGaussianFromPowers(const arma::rowvec &powers, double alpha, double *expFactor)
{
    gaussians* phi;
    gaussians* del_phi_x;
    gaussians* del_phi_y;
    gaussians* del_phi_z;
    gaussians* lapl_phi;

    //DO ALL THE MAGIC HERE...
    phi = new gaussians_000(expFactor, alpha);
    del_phi_x = new gaussians_000(expFactor, alpha);
    del_phi_y = new gaussians_000(expFactor, alpha);
    del_phi_z = new gaussians_000(expFactor, alpha);
    lapl_phi = new gaussians_000(expFactor, alpha);
    //

    primitivesPhi.push_back(phi);
    primitivesDelX.push_back(del_phi_x);
    primitivesDelY.push_back(del_phi_y);
    primitivesDelZ.push_back(del_phi_z);
    primitivesLapl.push_back(lapl_phi);


}

void GaussianFitted::addGaussianFitFromCGTOs(const ContractedGTO &CGTO)
{

    primitivesPhi.clear();
    primitivesDelX.clear();
    primitivesDelY.clear();
    primitivesDelZ.clear();
    primitivesLapl.clear();

    coeffs.reset();
    coeffs.set_size(CGTO.getNumPrimitives());

    for (int k = 0; k < CGTO.getNumPrimitives(); ++k) {

        const PrimitiveGTO & PGTO = CGTO.getPrimitive(k);

        coeffs(k) = PGTO.weight();
        getGaussianFromPGTO(PGTO);

    }

    basis_functions[CURRENT]         = new expandedBasisFunctions(primitivesPhi, coeffs);
    dell_basis_functions[0][CURRENT] = new expandedBasisFunctions(primitivesDelX, coeffs);
    dell_basis_functions[1][CURRENT] = new expandedBasisFunctions(primitivesDelY, coeffs);
    dell_basis_functions[2][CURRENT] = new expandedBasisFunctions(primitivesDelZ, coeffs);
    lapl_basis_functions[CURRENT]    = new expandedBasisFunctions(primitivesLapl, coeffs);

    CURRENT++;
}


void GaussianFitted::set_qnum_indie_terms(Walker *walker, int i)
{
    for (unsigned int k = 0; k < expFactors.size(); ++k) {
        *(expFactors.at(k)) = exp(-uniqueAlphas.at(k)*walker->get_r_i2(i));
    }
}




