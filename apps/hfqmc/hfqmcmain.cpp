#include <QMC2.h>

#include "hforbitals.h"

#include "HF/src/libs/hfSolver/hfsolver.h"
#include "HF/src/libs/basisSet/basisset.h"
#include "HF/src/libs/system/system.h"

#include <armadillo>

using namespace arma;

void runQMC(mat & coeffs, System * system);
void runHF(mat &coeffs, System *system);

int main(int argc, char** argv)
{

    (void) argc;
    (void) argv;

    System * system;
    mat coeffs;

    MPI::Init();

    // get the number of processes, and the id of this process
    int rank = MPI::COMM_WORLD.Get_rank();
    int nProcs = MPI::COMM_WORLD.Get_size();

    clock_t begin = clock();
    int nElectrons;
    rowvec coreCharges,coreMass, A, B, C, D, E;


    BasisSet *basisCoreA;
    BasisSet *basisCoreB;
    BasisSet *basisCoreC;
    BasisSet *basisCoreD;
    BasisSet *basisCoreE;


    int m_case = 1;

    if(m_case == 1){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-0.5, 0.0, 0.0};
        B = {0.5, 0.0, 0.0};
        coreCharges = {1 , 1};
        coreMass = {1 , 1};
        basisCoreA = new BasisSet("infiles/turbomole/H_Qzeta");
        basisCoreB = new BasisSet("infiles/turbomole/H_Qzeta");

    }else if(m_case==2){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-2.5, 0.0, 0.0};
        B = {2.5, 0.0, 0.0};
        coreCharges = {1 , 1};
        coreMass = {1 , 1};
        basisCoreA = new BasisSet("infiles/turbomole/H_3-21G");
        basisCoreB = new BasisSet("infiles/turbomole/H_3-21G");

    }else if(m_case==3){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-0.69, 0.0, 0.0};
        B = {0.69, 0.0, 0.0};
        coreCharges = {1 , 1};
        coreMass = {1 , 1};
        basisCoreA = new BasisSet("infiles/turbomole/H_4-31G");
        basisCoreB = new BasisSet("infiles/turbomole/H_4-31G");

    }else if(m_case==4){
        //Lithium molecule
        nElectrons = 6;
        A = {-2.5255, 0.0, 0.0};
        B = { 2.5255, 0.0, 0.0};
        coreCharges = {3 , 3};
        coreMass = {7 , 7};
        basisCoreA = new BasisSet("infiles/turbomole/Li_3-21G");
        basisCoreB = new BasisSet("infiles/turbomole/Li_3-21G");

    }else if(m_case==5){
        //Oxygen molecule
        nElectrons = 8;
        A = {-1.14, 0.0, 0.0};
        B = { 1.14, 0.0, 0.0};
        coreMass = {16 , 16};
        coreCharges = {8 , 8};
        basisCoreA = new BasisSet("infiles/turbomole/O_4-31G");
        basisCoreB = new BasisSet("infiles/turbomole/O_4-31G");

    }else if(m_case==6){
        double l = 1.797;
        double x = l*cos((180-104.45) *M_PI/180.0);
        double y = l*sin((180-104.45) *M_PI/180.0);

        //Water molecule
        nElectrons = 10;
        A = {l, 0.0, 0.0};
        B = { -x, y, 0.0};
        C = { 0.0, 0.0, 0.0};

        coreCharges = {1 , 1, 8};
        coreMass = {1 , 1, 16};

        basisCoreA = new BasisSet("infiles/turbomole/H_3-21G");
        basisCoreB = new BasisSet("infiles/turbomole/H_3-21G");
        basisCoreC = new BasisSet("infiles/turbomole/O_3-21G");


    }else if(m_case==7){
        //Carbon dioxide
        nElectrons = 22;
        A = {-2.2, 0.0, 0.0};
        B = { 2.2, 0.0, 0.0};
        C = { 0.0, 0.0, 0.0};
        coreCharges = {8 , 8, 6};
        coreMass = {16 , 16, 12};
        basisCoreA = new BasisSet("infiles/turbomole/O_3-21G");
        basisCoreB = new BasisSet("infiles/turbomole/O_3-21G");
        basisCoreC = new BasisSet("infiles/turbomole/C_3-21G");

    }else if(m_case==8){
        //SiO4
        nElectrons = 46;
        A = {0.0, 0.0, 0.0};
        B = {1.0, 1.0, 1.0};
        C = {-1.0, -1.0, 1.0};
        D = {1.0, -1.0, -1.0};
        E = {-1.0, 1.0, -1.0};

        B *=4.9/sqrt(3);C *=4.9/sqrt(3);D *=4.9/sqrt(3);E *=4.9/sqrt(3);

        coreCharges = {14, 8 , 8, 8, 8};
        coreMass = {28 , 16, 16, 16, 16};

        basisCoreA = new BasisSet("infiles/turbomole/Si_3-21G");
        basisCoreB = new BasisSet("infiles/turbomole/O_3-21G");
        basisCoreC = new BasisSet("infiles/turbomole/O_3-21G");
        basisCoreD = new BasisSet("infiles/turbomole/O_3-21G");
        basisCoreE = new BasisSet("infiles/turbomole/O_3-21G");

    }else if(m_case==9){
        //CH4
        nElectrons = 10;
        A = {0.0, 0.0, 0.0};
        B = {1.0, 1.0, 1.0};
        C = {-1.0, -1.0, 1.0};
        D = {1.0, -1.0, -1.0};
        E = {-1.0, 1.0, -1.0};

        B *=2.043/sqrt(3);C *=2.043/sqrt(3);D *=2.043/sqrt(3);E *=2.043/sqrt(3);

        coreCharges = {6, 1 , 1, 1, 1};
        coreMass = {6 , 1, 1, 1, 1};

        basisCoreA = new BasisSet("infiles/turbomole/C_4-31G");
        basisCoreB = new BasisSet("infiles/turbomole/H_4-31G");
        basisCoreC = new BasisSet("infiles/turbomole/H_4-31G");
        basisCoreD = new BasisSet("infiles/turbomole/H_4-31G");
        basisCoreE = new BasisSet("infiles/turbomole/H_4-31G");

    }else{
        cerr << "dead!" << endl;
    }
    /********************************************************/

    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    basisCoreA->setCoreCharge(coreCharges(0));
    basisCoreA->setCoreMass(coreMass(0));
    basisCoreA->setCorePosition(A);

    basisCoreB->setCorePosition(B);
    basisCoreB->setCoreCharge(coreCharges(1));
    basisCoreB->setCoreMass(coreMass(1));


    if(m_case == 6 ||m_case == 7 ){
        maxAngularMomentum = basisCoreC->getAngularMomentum();
        basisCoreC->setCorePosition(C);
        basisCoreC->setCoreCharge(coreCharges(2));
        basisCoreC->setCoreMass(coreMass(2));

    }

    if(m_case==8 ||m_case == 9){
        basisCoreC->setCorePosition(C);
        basisCoreC->setCoreCharge(coreCharges(2));
        basisCoreC->setCoreMass(coreMass(2));

        basisCoreD->setCorePosition(D);
        basisCoreD->setCoreCharge(coreCharges(3));
        basisCoreD->setCoreMass(coreMass(3));

        basisCoreE->setCorePosition(E);
        basisCoreE->setCoreCharge(coreCharges(4));
        basisCoreE->setCoreMass(coreMass(4));
    }


    system = new System(nElectrons, maxAngularMomentum);
    system->addBasisSet(basisCoreA);
    system->addBasisSet(basisCoreB);

    if(m_case == 6 || m_case == 7 ){
        system->addBasisSet(basisCoreC);
    }

    if(m_case == 8 ||m_case == 9 ){
        system->addBasisSet(basisCoreC);
        system->addBasisSet(basisCoreD);
        system->addBasisSet(basisCoreE);
    }

    HFsolver *solver = new HFsolver(system, rank, nProcs);
    solver->runSolver();
    coeffs = solver->getC();

    clock_t end = clock();
    if(rank==0){
        cout << "-------------------------------"  << endl;
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << "s" << endl;
    }

    runQMC(coeffs, system);


    MPI_Finalize();

    return 0;

}


void runQMC(mat &coeffs, System *system)
{

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
    QMC2::VMCparams vmcP;

    gP.n_p = system->getNumOfElectrons();

    HFOrbitals hfOrbitals(system);
    QMC2::ExpandedBasis expHFBasis(&hfOrbitals, coeffs);
    sO.SP_basis = &expHFBasis;


    sO.system = new QMC2::Fermions(gP, sO.SP_basis);

    QMC2::MolecularCoulomb MCO(gP);
    QMC2::Coulomb COL(gP);

    sO.system->add_potential(&MCO);
    sO.system->add_potential(&COL);

    //    sO.jastrow = new Pade_Jastrow(gP, vP);


    QMC2::VMC vmc(gP, vmcP, sO, parParams);
    vmc.set_error_estimator(new QMC2::SimpleVar(parParams));

    vmc.run_method();


}
