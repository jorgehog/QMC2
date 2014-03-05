#pragma once


#include "Orbitals.h"

#include "../structs.h"

#include "ExpandedBasis/ExpandedBasis.h"
#include "hydrogenicOrbitals/hydrogenicOrbitals.h"


namespace QMC2
{


enum FACTORYTYPE {EXPANDED,
                  ATOMS};


class OrbitalsFactory {
public:
    OrbitalsFactory(FACTORYTYPE type) :
        type(type)
    {

    }

    Orbitals* create(GeneralParams & gP, VariationalParams & vP) {

        Orbitals* basis;

        switch (type) {
        case EXPANDED:
        {
            ExpandedBasis* expbasis = new ExpandedBasis();
            expbasis->setBasis(basisForExpanded->create(gP, vP));
            expbasis->setCoeffs(C);

            basis = expbasis;

            break;
        }

        case ATOMS:
        {
            basis = new hydrogenicOrbitals(gP, vP);
            break;
        }

        default:

	    basis = NULL;
            std::cout << "UNREGISTERED ORBITAL MET" << std::endl;
            break;
        }

        return basis;
    }

    FACTORYTYPE type;

    OrbitalsFactory* basisForExpanded;
    arma::mat C;

};

}
