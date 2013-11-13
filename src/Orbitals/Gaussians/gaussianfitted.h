#ifndef GAUSSIANFITTED_H
#define GAUSSIANFITTED_H

#include "../Orbitals.h"

#include "BasisSetCodeMilad/contractedGTO/contractedGTO.h"

#include "BasisSetCodeMilad/splitValence/splitvalence.h"

#include <vector>

class gaussians;

class GaussianFitted : public Orbitals
{
public:
    GaussianFitted(int n_p, int dim, SplitValence * basis);

private:
    int CURRENT;

    std::vector<double*> expFactors;
    std::vector<double> uniqueAlphas;

    std::vector<BasisFunctions*> primitivesPhi;
    std::vector<BasisFunctions*> primitivesDelX;
    std::vector<BasisFunctions*> primitivesDelY;
    std::vector<BasisFunctions*> primitivesDelZ;
    std::vector<BasisFunctions*> primitivesLapl;

    arma::vec coeffs;

    void getGaussianFromPGTO(const PrimitiveGTO & PGTO);
    void locateGaussianFromPowers(const arma::rowvec & powers, double alpha, double * expFactor);
    void addGaussianFitFromCGTOs(const ContractedGTO & CGTO);

    // Orbitals interface
public:
    void set_qnum_indie_terms(Walker *walker, int i);



};

#endif // GAUSSIANFITTED_H
