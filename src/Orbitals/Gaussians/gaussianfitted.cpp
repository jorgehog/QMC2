#include "gaussianfitted.h"
#include "../../Walker/Walker.h"
#include "../../BasisFunctions/expandedbasisfunctions.h"
#include "../../BasisFunctions/Gaussians/gaussians.h"


GaussianFitted::GaussianFitted(int n_p, int dim, SplitValence * basis) :
    Orbitals(n_p, dim),
    CURRENT(0)
{
    name = "gaussianFitted";


    for (int i = 0; i < basis->getNumContracted(); ++i) {
        addGaussianFitFromCGTOs(basis->getContracted(i));
    }

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

void GaussianFitted::locateGaussianFromPowers(const arma::rowvec &powers, double alpha, double *expFactor)
{
    gaussians* phi;
    gaussians* del_phi_x;
    gaussians* del_phi_y;
    gaussians* del_phi_z;
    gaussians* lapl_phi;

    //DO ALL THE MAGIC HERE...
    int px = powers(0);
    int py = powers(1);
    int pz = powers(2);

    if ((px == 0) && (py == 0) && (pz == 0)) {
        phi = new gaussians_000(expFactor, alpha);
        del_phi_x = new dell_gaussians_000_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_000_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_000_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_000(expFactor, alpha);
    }
    if ((px == 0) && (py == 0) && (pz == 1)) {
        phi = new gaussians_001(expFactor, alpha);
        del_phi_x = new dell_gaussians_001_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_001_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_001_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_001(expFactor, alpha);
    }
    if ((px == 0) && (py == 1) && (pz == 0)) {
        phi = new gaussians_010(expFactor, alpha);
        del_phi_x = new dell_gaussians_010_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_010_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_010_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_010(expFactor, alpha);
    }
    if ((px == 1) && (py == 0) && (pz == 0)) {
        phi = new gaussians_100(expFactor, alpha);
        del_phi_x = new dell_gaussians_100_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_100_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_100_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_100(expFactor, alpha);
    }
    if ((px == 0) && (py == 0) && (pz == 2)) {
        phi = new gaussians_002(expFactor, alpha);
        del_phi_x = new dell_gaussians_002_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_002_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_002_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_002(expFactor, alpha);
    }
    if ((px == 0) && (py == 1) && (pz == 1)) {
        phi = new gaussians_011(expFactor, alpha);
        del_phi_x = new dell_gaussians_011_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_011_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_011_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_011(expFactor, alpha);
    }
    if ((px == 0) && (py == 2) && (pz == 0)) {
        phi = new gaussians_020(expFactor, alpha);
        del_phi_x = new dell_gaussians_020_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_020_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_020_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_020(expFactor, alpha);
    }
    if ((px == 1) && (py == 0) && (pz == 1)) {
        phi = new gaussians_101(expFactor, alpha);
        del_phi_x = new dell_gaussians_101_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_101_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_101_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_101(expFactor, alpha);
    }
    if ((px == 1) && (py == 1) && (pz == 0)) {
        phi = new gaussians_110(expFactor, alpha);
        del_phi_x = new dell_gaussians_110_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_110_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_110_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_110(expFactor, alpha);
    }
    if ((px == 2) && (py == 0) && (pz == 0)) {
        phi = new gaussians_200(expFactor, alpha);
        del_phi_x = new dell_gaussians_200_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_200_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_200_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_200(expFactor, alpha);
    }
    if ((px == 0) && (py == 0) && (pz == 3)) {
        phi = new gaussians_003(expFactor, alpha);
        del_phi_x = new dell_gaussians_003_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_003_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_003_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_003(expFactor, alpha);
    }
    if ((px == 0) && (py == 1) && (pz == 2)) {
        phi = new gaussians_012(expFactor, alpha);
        del_phi_x = new dell_gaussians_012_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_012_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_012_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_012(expFactor, alpha);
    }
    if ((px == 0) && (py == 2) && (pz == 1)) {
        phi = new gaussians_021(expFactor, alpha);
        del_phi_x = new dell_gaussians_021_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_021_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_021_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_021(expFactor, alpha);
    }
    if ((px == 0) && (py == 3) && (pz == 0)) {
        phi = new gaussians_030(expFactor, alpha);
        del_phi_x = new dell_gaussians_030_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_030_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_030_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_030(expFactor, alpha);
    }
    if ((px == 1) && (py == 0) && (pz == 2)) {
        phi = new gaussians_102(expFactor, alpha);
        del_phi_x = new dell_gaussians_102_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_102_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_102_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_102(expFactor, alpha);
    }
    if ((px == 1) && (py == 1) && (pz == 1)) {
        phi = new gaussians_111(expFactor, alpha);
        del_phi_x = new dell_gaussians_111_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_111_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_111_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_111(expFactor, alpha);
    }
    if ((px == 1) && (py == 2) && (pz == 0)) {
        phi = new gaussians_120(expFactor, alpha);
        del_phi_x = new dell_gaussians_120_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_120_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_120_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_120(expFactor, alpha);
    }
    if ((px == 2) && (py == 0) && (pz == 1)) {
        phi = new gaussians_201(expFactor, alpha);
        del_phi_x = new dell_gaussians_201_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_201_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_201_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_201(expFactor, alpha);
    }
    if ((px == 2) && (py == 1) && (pz == 0)) {
        phi = new gaussians_210(expFactor, alpha);
        del_phi_x = new dell_gaussians_210_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_210_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_210_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_210(expFactor, alpha);
    }
    if ((px == 3) && (py == 0) && (pz == 0)) {
        phi = new gaussians_300(expFactor, alpha);
        del_phi_x = new dell_gaussians_300_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_300_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_300_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_300(expFactor, alpha);
    }
    if ((px == 0) && (py == 0) && (pz == 4)) {
        phi = new gaussians_004(expFactor, alpha);
        del_phi_x = new dell_gaussians_004_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_004_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_004_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_004(expFactor, alpha);
    }
    if ((px == 0) && (py == 1) && (pz == 3)) {
        phi = new gaussians_013(expFactor, alpha);
        del_phi_x = new dell_gaussians_013_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_013_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_013_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_013(expFactor, alpha);
    }
    if ((px == 0) && (py == 2) && (pz == 2)) {
        phi = new gaussians_022(expFactor, alpha);
        del_phi_x = new dell_gaussians_022_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_022_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_022_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_022(expFactor, alpha);
    }
    if ((px == 0) && (py == 3) && (pz == 1)) {
        phi = new gaussians_031(expFactor, alpha);
        del_phi_x = new dell_gaussians_031_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_031_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_031_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_031(expFactor, alpha);
    }
    if ((px == 0) && (py == 4) && (pz == 0)) {
        phi = new gaussians_040(expFactor, alpha);
        del_phi_x = new dell_gaussians_040_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_040_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_040_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_040(expFactor, alpha);
    }
    if ((px == 1) && (py == 0) && (pz == 3)) {
        phi = new gaussians_103(expFactor, alpha);
        del_phi_x = new dell_gaussians_103_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_103_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_103_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_103(expFactor, alpha);
    }
    if ((px == 1) && (py == 1) && (pz == 2)) {
        phi = new gaussians_112(expFactor, alpha);
        del_phi_x = new dell_gaussians_112_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_112_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_112_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_112(expFactor, alpha);
    }
    if ((px == 1) && (py == 2) && (pz == 1)) {
        phi = new gaussians_121(expFactor, alpha);
        del_phi_x = new dell_gaussians_121_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_121_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_121_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_121(expFactor, alpha);
    }
    if ((px == 1) && (py == 3) && (pz == 0)) {
        phi = new gaussians_130(expFactor, alpha);
        del_phi_x = new dell_gaussians_130_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_130_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_130_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_130(expFactor, alpha);
    }
    if ((px == 2) && (py == 0) && (pz == 2)) {
        phi = new gaussians_202(expFactor, alpha);
        del_phi_x = new dell_gaussians_202_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_202_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_202_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_202(expFactor, alpha);
    }
    if ((px == 2) && (py == 1) && (pz == 1)) {
        phi = new gaussians_211(expFactor, alpha);
        del_phi_x = new dell_gaussians_211_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_211_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_211_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_211(expFactor, alpha);
    }
    if ((px == 2) && (py == 2) && (pz == 0)) {
        phi = new gaussians_220(expFactor, alpha);
        del_phi_x = new dell_gaussians_220_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_220_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_220_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_220(expFactor, alpha);
    }
    if ((px == 3) && (py == 0) && (pz == 1)) {
        phi = new gaussians_301(expFactor, alpha);
        del_phi_x = new dell_gaussians_301_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_301_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_301_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_301(expFactor, alpha);
    }
    if ((px == 3) && (py == 1) && (pz == 0)) {
        phi = new gaussians_310(expFactor, alpha);
        del_phi_x = new dell_gaussians_310_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_310_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_310_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_310(expFactor, alpha);
    }
    if ((px == 4) && (py == 0) && (pz == 0)) {
        phi = new gaussians_400(expFactor, alpha);
        del_phi_x = new dell_gaussians_400_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_400_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_400_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_400(expFactor, alpha);
    }
    if ((px == 0) && (py == 0) && (pz == 5)) {
        phi = new gaussians_005(expFactor, alpha);
        del_phi_x = new dell_gaussians_005_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_005_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_005_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_005(expFactor, alpha);
    }
    if ((px == 0) && (py == 1) && (pz == 4)) {
        phi = new gaussians_014(expFactor, alpha);
        del_phi_x = new dell_gaussians_014_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_014_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_014_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_014(expFactor, alpha);
    }
    if ((px == 0) && (py == 2) && (pz == 3)) {
        phi = new gaussians_023(expFactor, alpha);
        del_phi_x = new dell_gaussians_023_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_023_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_023_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_023(expFactor, alpha);
    }
    if ((px == 0) && (py == 3) && (pz == 2)) {
        phi = new gaussians_032(expFactor, alpha);
        del_phi_x = new dell_gaussians_032_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_032_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_032_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_032(expFactor, alpha);
    }
    if ((px == 0) && (py == 4) && (pz == 1)) {
        phi = new gaussians_041(expFactor, alpha);
        del_phi_x = new dell_gaussians_041_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_041_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_041_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_041(expFactor, alpha);
    }
    if ((px == 0) && (py == 5) && (pz == 0)) {
        phi = new gaussians_050(expFactor, alpha);
        del_phi_x = new dell_gaussians_050_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_050_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_050_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_050(expFactor, alpha);
    }
    if ((px == 1) && (py == 0) && (pz == 4)) {
        phi = new gaussians_104(expFactor, alpha);
        del_phi_x = new dell_gaussians_104_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_104_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_104_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_104(expFactor, alpha);
    }
    if ((px == 1) && (py == 1) && (pz == 3)) {
        phi = new gaussians_113(expFactor, alpha);
        del_phi_x = new dell_gaussians_113_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_113_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_113_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_113(expFactor, alpha);
    }
    if ((px == 1) && (py == 2) && (pz == 2)) {
        phi = new gaussians_122(expFactor, alpha);
        del_phi_x = new dell_gaussians_122_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_122_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_122_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_122(expFactor, alpha);
    }
    if ((px == 1) && (py == 3) && (pz == 1)) {
        phi = new gaussians_131(expFactor, alpha);
        del_phi_x = new dell_gaussians_131_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_131_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_131_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_131(expFactor, alpha);
    }
    if ((px == 1) && (py == 4) && (pz == 0)) {
        phi = new gaussians_140(expFactor, alpha);
        del_phi_x = new dell_gaussians_140_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_140_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_140_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_140(expFactor, alpha);
    }
    if ((px == 2) && (py == 0) && (pz == 3)) {
        phi = new gaussians_203(expFactor, alpha);
        del_phi_x = new dell_gaussians_203_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_203_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_203_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_203(expFactor, alpha);
    }
    if ((px == 2) && (py == 1) && (pz == 2)) {
        phi = new gaussians_212(expFactor, alpha);
        del_phi_x = new dell_gaussians_212_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_212_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_212_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_212(expFactor, alpha);
    }
    if ((px == 2) && (py == 2) && (pz == 1)) {
        phi = new gaussians_221(expFactor, alpha);
        del_phi_x = new dell_gaussians_221_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_221_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_221_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_221(expFactor, alpha);
    }
    if ((px == 2) && (py == 3) && (pz == 0)) {
        phi = new gaussians_230(expFactor, alpha);
        del_phi_x = new dell_gaussians_230_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_230_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_230_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_230(expFactor, alpha);
    }
    if ((px == 3) && (py == 0) && (pz == 2)) {
        phi = new gaussians_302(expFactor, alpha);
        del_phi_x = new dell_gaussians_302_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_302_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_302_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_302(expFactor, alpha);
    }
    if ((px == 3) && (py == 1) && (pz == 1)) {
        phi = new gaussians_311(expFactor, alpha);
        del_phi_x = new dell_gaussians_311_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_311_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_311_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_311(expFactor, alpha);
    }
    if ((px == 3) && (py == 2) && (pz == 0)) {
        phi = new gaussians_320(expFactor, alpha);
        del_phi_x = new dell_gaussians_320_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_320_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_320_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_320(expFactor, alpha);
    }
    if ((px == 4) && (py == 0) && (pz == 1)) {
        phi = new gaussians_401(expFactor, alpha);
        del_phi_x = new dell_gaussians_401_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_401_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_401_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_401(expFactor, alpha);
    }
    if ((px == 4) && (py == 1) && (pz == 0)) {
        phi = new gaussians_410(expFactor, alpha);
        del_phi_x = new dell_gaussians_410_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_410_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_410_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_410(expFactor, alpha);
    }
    if ((px == 5) && (py == 0) && (pz == 0)) {
        phi = new gaussians_500(expFactor, alpha);
        del_phi_x = new dell_gaussians_500_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_500_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_500_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_500(expFactor, alpha);
    }
    if ((px == 0) && (py == 0) && (pz == 6)) {
        phi = new gaussians_006(expFactor, alpha);
        del_phi_x = new dell_gaussians_006_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_006_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_006_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_006(expFactor, alpha);
    }
    if ((px == 0) && (py == 1) && (pz == 5)) {
        phi = new gaussians_015(expFactor, alpha);
        del_phi_x = new dell_gaussians_015_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_015_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_015_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_015(expFactor, alpha);
    }
    if ((px == 0) && (py == 2) && (pz == 4)) {
        phi = new gaussians_024(expFactor, alpha);
        del_phi_x = new dell_gaussians_024_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_024_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_024_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_024(expFactor, alpha);
    }
    if ((px == 0) && (py == 3) && (pz == 3)) {
        phi = new gaussians_033(expFactor, alpha);
        del_phi_x = new dell_gaussians_033_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_033_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_033_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_033(expFactor, alpha);
    }
    if ((px == 0) && (py == 4) && (pz == 2)) {
        phi = new gaussians_042(expFactor, alpha);
        del_phi_x = new dell_gaussians_042_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_042_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_042_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_042(expFactor, alpha);
    }
    if ((px == 0) && (py == 5) && (pz == 1)) {
        phi = new gaussians_051(expFactor, alpha);
        del_phi_x = new dell_gaussians_051_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_051_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_051_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_051(expFactor, alpha);
    }
    if ((px == 0) && (py == 6) && (pz == 0)) {
        phi = new gaussians_060(expFactor, alpha);
        del_phi_x = new dell_gaussians_060_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_060_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_060_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_060(expFactor, alpha);
    }
    if ((px == 1) && (py == 0) && (pz == 5)) {
        phi = new gaussians_105(expFactor, alpha);
        del_phi_x = new dell_gaussians_105_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_105_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_105_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_105(expFactor, alpha);
    }
    if ((px == 1) && (py == 1) && (pz == 4)) {
        phi = new gaussians_114(expFactor, alpha);
        del_phi_x = new dell_gaussians_114_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_114_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_114_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_114(expFactor, alpha);
    }
    if ((px == 1) && (py == 2) && (pz == 3)) {
        phi = new gaussians_123(expFactor, alpha);
        del_phi_x = new dell_gaussians_123_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_123_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_123_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_123(expFactor, alpha);
    }
    if ((px == 1) && (py == 3) && (pz == 2)) {
        phi = new gaussians_132(expFactor, alpha);
        del_phi_x = new dell_gaussians_132_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_132_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_132_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_132(expFactor, alpha);
    }
    if ((px == 1) && (py == 4) && (pz == 1)) {
        phi = new gaussians_141(expFactor, alpha);
        del_phi_x = new dell_gaussians_141_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_141_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_141_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_141(expFactor, alpha);
    }
    if ((px == 1) && (py == 5) && (pz == 0)) {
        phi = new gaussians_150(expFactor, alpha);
        del_phi_x = new dell_gaussians_150_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_150_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_150_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_150(expFactor, alpha);
    }
    if ((px == 2) && (py == 0) && (pz == 4)) {
        phi = new gaussians_204(expFactor, alpha);
        del_phi_x = new dell_gaussians_204_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_204_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_204_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_204(expFactor, alpha);
    }
    if ((px == 2) && (py == 1) && (pz == 3)) {
        phi = new gaussians_213(expFactor, alpha);
        del_phi_x = new dell_gaussians_213_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_213_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_213_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_213(expFactor, alpha);
    }
    if ((px == 2) && (py == 2) && (pz == 2)) {
        phi = new gaussians_222(expFactor, alpha);
        del_phi_x = new dell_gaussians_222_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_222_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_222_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_222(expFactor, alpha);
    }
    if ((px == 2) && (py == 3) && (pz == 1)) {
        phi = new gaussians_231(expFactor, alpha);
        del_phi_x = new dell_gaussians_231_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_231_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_231_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_231(expFactor, alpha);
    }
    if ((px == 2) && (py == 4) && (pz == 0)) {
        phi = new gaussians_240(expFactor, alpha);
        del_phi_x = new dell_gaussians_240_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_240_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_240_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_240(expFactor, alpha);
    }
    if ((px == 3) && (py == 0) && (pz == 3)) {
        phi = new gaussians_303(expFactor, alpha);
        del_phi_x = new dell_gaussians_303_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_303_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_303_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_303(expFactor, alpha);
    }
    if ((px == 3) && (py == 1) && (pz == 2)) {
        phi = new gaussians_312(expFactor, alpha);
        del_phi_x = new dell_gaussians_312_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_312_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_312_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_312(expFactor, alpha);
    }
    if ((px == 3) && (py == 2) && (pz == 1)) {
        phi = new gaussians_321(expFactor, alpha);
        del_phi_x = new dell_gaussians_321_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_321_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_321_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_321(expFactor, alpha);
    }
    if ((px == 3) && (py == 3) && (pz == 0)) {
        phi = new gaussians_330(expFactor, alpha);
        del_phi_x = new dell_gaussians_330_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_330_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_330_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_330(expFactor, alpha);
    }
    if ((px == 4) && (py == 0) && (pz == 2)) {
        phi = new gaussians_402(expFactor, alpha);
        del_phi_x = new dell_gaussians_402_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_402_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_402_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_402(expFactor, alpha);
    }
    if ((px == 4) && (py == 1) && (pz == 1)) {
        phi = new gaussians_411(expFactor, alpha);
        del_phi_x = new dell_gaussians_411_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_411_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_411_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_411(expFactor, alpha);
    }
    if ((px == 4) && (py == 2) && (pz == 0)) {
        phi = new gaussians_420(expFactor, alpha);
        del_phi_x = new dell_gaussians_420_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_420_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_420_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_420(expFactor, alpha);
    }
    if ((px == 5) && (py == 0) && (pz == 1)) {
        phi = new gaussians_501(expFactor, alpha);
        del_phi_x = new dell_gaussians_501_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_501_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_501_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_501(expFactor, alpha);
    }
    if ((px == 5) && (py == 1) && (pz == 0)) {
        phi = new gaussians_510(expFactor, alpha);
        del_phi_x = new dell_gaussians_510_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_510_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_510_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_510(expFactor, alpha);
    }
    if ((px == 6) && (py == 0) && (pz == 0)) {
        phi = new gaussians_600(expFactor, alpha);
        del_phi_x = new dell_gaussians_600_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_600_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_600_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_600(expFactor, alpha);
    }
    if ((px == 0) && (py == 0) && (pz == 7)) {
        phi = new gaussians_007(expFactor, alpha);
        del_phi_x = new dell_gaussians_007_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_007_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_007_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_007(expFactor, alpha);
    }
    if ((px == 0) && (py == 1) && (pz == 6)) {
        phi = new gaussians_016(expFactor, alpha);
        del_phi_x = new dell_gaussians_016_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_016_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_016_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_016(expFactor, alpha);
    }
    if ((px == 0) && (py == 2) && (pz == 5)) {
        phi = new gaussians_025(expFactor, alpha);
        del_phi_x = new dell_gaussians_025_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_025_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_025_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_025(expFactor, alpha);
    }
    if ((px == 0) && (py == 3) && (pz == 4)) {
        phi = new gaussians_034(expFactor, alpha);
        del_phi_x = new dell_gaussians_034_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_034_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_034_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_034(expFactor, alpha);
    }
    if ((px == 0) && (py == 4) && (pz == 3)) {
        phi = new gaussians_043(expFactor, alpha);
        del_phi_x = new dell_gaussians_043_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_043_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_043_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_043(expFactor, alpha);
    }
    if ((px == 0) && (py == 5) && (pz == 2)) {
        phi = new gaussians_052(expFactor, alpha);
        del_phi_x = new dell_gaussians_052_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_052_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_052_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_052(expFactor, alpha);
    }
    if ((px == 0) && (py == 6) && (pz == 1)) {
        phi = new gaussians_061(expFactor, alpha);
        del_phi_x = new dell_gaussians_061_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_061_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_061_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_061(expFactor, alpha);
    }
    if ((px == 0) && (py == 7) && (pz == 0)) {
        phi = new gaussians_070(expFactor, alpha);
        del_phi_x = new dell_gaussians_070_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_070_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_070_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_070(expFactor, alpha);
    }
    if ((px == 1) && (py == 0) && (pz == 6)) {
        phi = new gaussians_106(expFactor, alpha);
        del_phi_x = new dell_gaussians_106_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_106_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_106_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_106(expFactor, alpha);
    }
    if ((px == 1) && (py == 1) && (pz == 5)) {
        phi = new gaussians_115(expFactor, alpha);
        del_phi_x = new dell_gaussians_115_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_115_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_115_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_115(expFactor, alpha);
    }
    if ((px == 1) && (py == 2) && (pz == 4)) {
        phi = new gaussians_124(expFactor, alpha);
        del_phi_x = new dell_gaussians_124_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_124_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_124_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_124(expFactor, alpha);
    }
    if ((px == 1) && (py == 3) && (pz == 3)) {
        phi = new gaussians_133(expFactor, alpha);
        del_phi_x = new dell_gaussians_133_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_133_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_133_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_133(expFactor, alpha);
    }
    if ((px == 1) && (py == 4) && (pz == 2)) {
        phi = new gaussians_142(expFactor, alpha);
        del_phi_x = new dell_gaussians_142_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_142_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_142_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_142(expFactor, alpha);
    }
    if ((px == 1) && (py == 5) && (pz == 1)) {
        phi = new gaussians_151(expFactor, alpha);
        del_phi_x = new dell_gaussians_151_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_151_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_151_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_151(expFactor, alpha);
    }
    if ((px == 1) && (py == 6) && (pz == 0)) {
        phi = new gaussians_160(expFactor, alpha);
        del_phi_x = new dell_gaussians_160_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_160_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_160_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_160(expFactor, alpha);
    }
    if ((px == 2) && (py == 0) && (pz == 5)) {
        phi = new gaussians_205(expFactor, alpha);
        del_phi_x = new dell_gaussians_205_x(expFactor, alpha);
        del_phi_y = new dell_gaussians_205_y(expFactor, alpha);
        del_phi_z = new dell_gaussians_205_z(expFactor, alpha);
        lapl_phi = new lapl_gaussians_205(expFactor, alpha);
    }

    //

    primitivesPhi.push_back(phi);
    primitivesDelX.push_back(del_phi_x);
    primitivesDelY.push_back(del_phi_y);
    primitivesDelZ.push_back(del_phi_z);
    primitivesLapl.push_back(lapl_phi);


}






