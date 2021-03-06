#include "System.h"

#include <math.h>
#include <sstream>
#include <iomanip>

#include "../Potential/Potential.h"
#include "../ErrorEstimator/ErrorEstimator.h"


using namespace QMC2;

System::System(int n_p, int dim, Orbitals* orbital) {
    this->n_p = n_p;
    n2 = ceil(n_p / 2.0);
    this->dim = dim;

    this->orbital = orbital;
}

System::System() {

}

void System::add_potential(Potential* pot) {
    potentials.push_back(pot);
}

void System::forEachPotentialDo(std::function<void(Potential *pot)> func)
{
    for (Potential *pot : potentials)
    {
        func(pot);
    }
}

double System::get_potential_energy(const Walker* walker) {
    double potE = 0;
    double potE_i;

    for (Potential *pot : potentials)
    {
        potE_i = pot->get_pot_E(walker);
        pot->pot_sampler->push_value(potE_i);
        potE += potE_i;
    }

    return potE;
}

void System::update_potential_samples(double weight) {

    for (Potential *pot : potentials)
    {
        pot->pot_sampler->update_mean(weight);
    }
}

void System::push_potential_samples() {

    for (Potential *pot : potentials)
    {
        pot->pot_sampler->push_mean();
    }
}

void System::initializeSamplingErrorEstimators(int type, int n_c)
{
    for (Potential *pot : potentials)
    {
        pot->pot_sampler->initializeErrorEstimator(type, n_c);
    }
}

std::string System::dump_samples()
{
    using namespace std;

    stringstream s;

    s << setprecision(6) << fixed;
    
    for (Potential *pot : potentials) {
        
        s << pot->get_name() << " " << pot->pot_sampler->result() << endl;

    }

    return s.str();

}

void System::reset_potential_samples()
{
    for (Potential * potential : potentials)
    {
        potential->pot_sampler->reset();
    }
}

void System::finalize_potential_samples()
{
    for (Potential * potential : potentials)
    {
        potential->pot_sampler->getErrorEstimator()->finalize();
    }
}
