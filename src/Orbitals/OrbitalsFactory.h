#pragma once


#include "Orbitals.h"

#include "../structs.h"

#include "ExpandedBasis/ExpandedBasis.h"
#include "Gaussians/oxygen3-21G/oxygen3_21g.h"
#include "hydrogenicOrbitals/hydrogenicOrbitals.h"

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
            ExpandedBasis* expbasis = new ExpandedBasis(gP);
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
