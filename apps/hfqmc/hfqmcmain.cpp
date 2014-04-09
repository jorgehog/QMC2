#include <QMC2.h>
#include <hf.h>
#include "hforbitals.h"

#include <armadillo>

using namespace arma;
using namespace hf;

ElectronicSystem *setupSystem(string name);

int main(int argc, char** argv)
{
    boost::mpi::environment env;
    boost::mpi::communicator world;
    boost::mpi::timer timer;
    timer.restart();

    (void) argc;
    (void) argv;

    //system:
    string chemicalSystem = "H2";

    if(world.rank() == 0){
            cout << "---------------------------Hartree-Fock------------------------------"  << endl;
            cout << "system:    " << chemicalSystem << endl;
        }

    //Setup system:
    ElectronicSystem *system = setupSystem(chemicalSystem);

    //Choose method:
    HFsolver* solver;
    solver = new RHF(system);

    //run solver:
    solver->runSolver();

    //get expansion coefficients:
    //usage: - row k contains the coeffs of orbital k
    //       - the first N/2 cols are the coeffs of the occupied orbitals
    //       - the last  N/2 cols are the coeffs of virtual orbitals
    const mat* allCoeffs = solver->expansionCoefficients().at(0);
    const mat coeffs = allCoeffs->cols(0, system->nElectrons()/2 - 1);
    //----------------------------------------------------------------------------------------------//
    vector<const rowvec*> corePositions;
    rowvec coreCharges = zeros<rowvec>(system->nAtoms());

    int i = 0;
    for (const Atom* atom : system->atoms())
    {
        corePositions.push_back(&atom->corePosition());
        coreCharges[i] = atom->coreCharge();
        i++;
    }


    QMC2::ParParams parParams;
    int n_nodes, node;
    MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);

    parParams.n_nodes = n_nodes;
    parParams.node = node;
    parParams.parallel = (parParams.n_nodes > 1);
    parParams.is_master = (parParams.node == 0);


    QMC2::GeneralParams gP;
    QMC2::SystemObjects sO;
    QMC2::MinimizerParams mP;
    QMC2::VMCparams vmcP;
    QMC2::VariationalParams vP;
    QMC2::DMCparams dP;

    dP.n_b = 100;
    dP.n_w = 1000;
    dP.therm = 100;

    vmcP.n_c = 1E6;
    mP.n_c_SGD = 2000;
    mP.n_c = 200;

    QMC2::scaleWithProcs(parParams, gP, mP, vmcP, dP);

    gP.n_p = system->nElectrons();
    gP.dim = 3;

    QMC2::HFOrbitals hfOrbitals(gP.n_p, system->basisFunctions(), coeffs);
//    QMC2::ExpandedBasis expHFBasis(&hfOrbitals, coeffs.t());
//    sO.SP_basis = &expHFBasis;
    sO.SP_basis = &hfOrbitals;

    sO.system = new QMC2::Fermions(gP, sO.SP_basis);


    QMC2::MolecularCoulomb MCO(gP, corePositions, coreCharges);
    QMC2::Coulomb COL(gP);

    sO.system->add_potential(&MCO);
    sO.system->add_potential(&COL);

    vP.beta = 0.53;
    sO.jastrow = new QMC2::Pade_Jastrow(gP, vP);


    QMC2::VMC vmc(gP, vmcP, sO, parParams, dP.n_w);
    QMC2::DMC dmc(gP, dP, sO, parParams);

    mP.alpha.reset();
    QMC2::ASGD asgd(&vmc, mP, parParams, gP.runpath);

    vmc.set_error_estimator(new QMC2::SimpleVar(parParams));
    dmc.set_error_estimator(new QMC2::SimpleVar(parParams));

//    asgd.minimize();
    vmc.run_method();
    vmc.dump_subsamples();

    dmc.initFromVMC(&vmc);

    dmc.run_method();

    MPI_Finalize();

    return 0;
}


ElectronicSystem* setupSystem(string name)
{
    vector<Atom *> atoms;

    if(name =="H2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { -0.69, 0, 0 }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { 0.69, 0, 0 }));

    }else if(name =="Li2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", {-2.5255, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", { 2.5255, 0.0, 0.0}));

    }else if(name =="O2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-1.14, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 1.14, 0.0, 0.0}));

    }else if(name =="H2O"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {1.797, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -1.797*cos((180-104.45) *M_PI/180.0),
                                                                              1.797*sin((180-104.45) *M_PI/180.0), 0.0}));
    }else if(name =="CO2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.0, 0.0, 0.0}));

    }else if(name =="CH4"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {2.043/sqrt(3), 2.043/sqrt(3), 2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-2.043/sqrt(3), -2.043/sqrt(3), 2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {2.043/sqrt(3), -2.043/sqrt(3), -2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-2.043/sqrt(3), 2.043/sqrt(3), -2.043/sqrt(3)}));

    }else if(name =="SiO4"){
        double D = 4.9;
        double T = 2.0*D /sqrt(3);
        atoms.push_back(new Atom("infiles/turbomole/atom_14_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {D/sqrt(3), D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-D/sqrt(3), -D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{D/sqrt(3), -D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-D/sqrt(3), D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_14_basis_3-21G.tm", {T, -T, -T}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {T - D/sqrt(3), -T -D/sqrt(3), -T - D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{T + D/sqrt(3), -T -D/sqrt(3), -T + D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{T + D/sqrt(3), -T -D/sqrt(3), -T + D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {T + D/sqrt(3), -T +D/sqrt(3), -T -D/sqrt(3)}));

    }else if(name =="Fe2S2"){
        double FeFe = 4.556129358;
        double SS   = 6.908838214;
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31G.tm", {0.0,  FeFe*0.5, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31G.tm", {0.0, -FeFe*0.5, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_6-31G.tm", {-SS*0.5, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_6-31G.tm", {SS*0.5, 0.0, 0.0}));

    }else if(name =="FeS2"){
        double FeS = 3.802128689;
        double SFeS = 114.0;
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31Gd.tm", {0.0, 0.0 , 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21Gd.tm", {FeS, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21Gd.tm",  {-FeS*cos((180.0-SFeS) *M_PI/180.0),
                                                                               FeS*sin((180.0-SFeS) *M_PI/180.0), 0.0}));

    }else if(name =="benzene"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        ,  2.63805748,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872,  1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872, -1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        , -2.63805748,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872, -1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872,  1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        ,  4.68463073,  0.         }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 ,  2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 , -2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        , -4.68463073,   0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 , -2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 ,  2.34326023,  0.        }));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

    return system;

}

