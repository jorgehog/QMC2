#include "molecularcoulomb.h"

#include "../../structs.h"
#include "../../Orbitals/NBodyTransform/nbodytransform.h"
#include "../../Walker/Walker.h"


using namespace QMC2;

MolecularCoulomb::MolecularCoulomb(GeneralParams &gP, const vector<const rowvec*> & corePositions, const rowvec & coreCharge) :
    Potential(gP.n_p, gP.dim),
    NCores(coreCharge.n_elem),
    corePositions(corePositions),
    coreCharge(coreCharge)

{
    name = "MolecularCoulomb";

    makeRRelNucleiMatrix();
}

double MolecularCoulomb::get_pot_E(const Walker *walker) const
{

    (void) walker;

    double potE = 0;

    double potE_i1;
    double potE_i2;


    for (int k = 0; k < (int)NCores; ++k) {

        potE_i1 = 0;


        //Nucleus - electron interaction
        for (int i = 0; i < n_p; ++i) {
            potE_i1 -= 1.0/arma::norm(walker->r.row(i) - *corePositions.at(k), 2);
        }

        potE_i2 = 0;

        //Nucleus - nucleus interaction
        for (int i = k + 1; i < NCores; ++i) {
            potE_i2 += coreCharge(i)/r_rel_nuclei(k, i);
        }

        potE += (potE_i1 + potE_i2)*coreCharge(k);

    }

    return potE;

}

void MolecularCoulomb::makeRRelNucleiMatrix()
{
    r_rel_nuclei.set_size(NCores, NCores);

    for (int i = 0; i < NCores - 1; i++) {
        for (int j = i + 1; j < NCores; j++) {
            r_rel_nuclei(i, j) = r_rel_nuclei(j, i) = arma::norm(*corePositions.at(i) - *corePositions.at(j), 2);
        }
    }


//    cout << r_rel_nuclei << endl;
//    exit(1);

}
