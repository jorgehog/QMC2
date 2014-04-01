#include <QMC2.h>
#include <hf.h>
#include "hforbitals.h"

#include <armadillo>

using namespace arma;
using namespace hf;

int main(int argc, char** argv)
{
    boost::mpi::environment env;
    boost::mpi::communicator world;
    boost::mpi::timer timer;
    timer.restart();

    (void) argc;
    (void) argv;

    vector<Atom *> atoms;

    atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { -0.69, 0, 0 }));
    atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { 0.69, 0, 0 }));
    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

    //Choose method:
//    HFsolver* solver;
//    solver = new RHF(system);

//    solver->runSolver();

    mat coeffs = "0.3285   0.1207   0.7620   1.1308";

    // get the number of processes, and the id of this process




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

    dP.n_b = 20;
    dP.n_w = 1000;

    vmcP.n_c = 1E7;
    mP.n_c_SGD = 2000;
    mP.n_c = 200;

    QMC2::scaleWithProcs(parParams, gP, mP, vmcP, dP);

    gP.n_p = system->nElectrons();
    gP.dim = 3;

    QMC2::HFOrbitals hfOrbitals(gP.n_p, system->basisFunctions());
//    QMC2::ExpandedBasis expHFBasis(&hfOrbitals, coeffs.t());
//    sO.SP_basis = &expHFBasis;
    sO.SP_basis = &hfOrbitals;

    sO.system = new QMC2::Fermions(gP, sO.SP_basis);

    vector<const rowvec*> corePositions;

    for (const Atom* atom : system->atoms())
    {
        corePositions.push_back(&atom->corePosition());
    }

    rowvec coreCharges = {1,1};

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

//    asgd.minimize();
    vmc.run_method();
    vmc.dump_subsamples();

    dmc.initFromVMC(&vmc);

    dmc.run_method();



    MPI_Finalize();

    return 0;

}


