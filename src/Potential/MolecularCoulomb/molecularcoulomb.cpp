#include "molecularcoulomb.h"

#include "../../structs.h"
#include "../../Orbitals/NBodyTransform/nbodytransform.h"
#include "../../Walker/Walker.h"


using namespace QMC2;

MolecularCoulomb::MolecularCoulomb(GeneralParams &gP) :
    Potential(gP.n_p, gP.dim)
{
    name = "MolecularCoulomb";
}

double MolecularCoulomb::get_pot_E(const Walker *walker) const
{

    (void) walker;

    double potE = 0;

    double potE_i1;
    double potE_i2;

    for (int k = 0; k < nBodyOrbitals->N; ++k) {

        potE_i1 = 0;


        //Nucleus - electron interaction
        for (int i = 0; i < n_p; ++i) {
//            potE_i1 -= 1.0/arma::norm(walker->r.row(i) - nBodyOrbitals->origins.at(k), 2);
            potE_i1 -= 1.0/nBodyOrbitals->nuclei_walkers.at(k)->get_r_i(i);
        }

        potE_i2 = 0;

        //Nucleus - nucleus interaction
        for (int i = k + 1; i < nBodyOrbitals->N; ++i) {
            potE_i2 += nBodyOrbitals->populations(i)/nBodyOrbitals->r_rel_nuclei(k, i);
        }

        potE += (potE_i1 + potE_i2)*nBodyOrbitals->populations(k);

    }

    return potE;

}
