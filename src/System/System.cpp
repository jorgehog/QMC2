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

double System::get_potential_energy(const Walker* walker) {
    double potE = 0;
    double potE_i;

    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        potE_i = (*pot)->get_pot_E(walker);
        (*pot)->pot_sampler->queue_value(potE_i);
        potE += potE_i;
    }

    return potE;
}

void System::update_potential_samples(double weight) {
    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        (*pot)->pot_sampler->update_mean(weight);
    }
}

void System::push_potential_samples() {
    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        (*pot)->pot_sampler->push_mean();
    }
}

void System::setMeanErrorEstimatorNumberOfCycles(int n_c)
{
    for (Potential *pot : potentials)
    {
        pot->pot_sampler->getMeanErrorEstimator()->setNumberOfCycles(n_c);
    }
}

void System::setMeanOfMeansErrorEstimatorNumberOfCycles(int n_c)
{
    for (Potential *pot : potentials)
    {
        pot->pot_sampler->getMeanOfMeansErrorEstimator()->setNumberOfCycles(n_c);
    }
}

std::string System::dump_samples(bool mean_of_means) {
    using namespace std;

    stringstream s;

    s << setprecision(6) << fixed;
    
    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        
        s << (*pot)->get_name() << " ";

        if (mean_of_means) {

            s << (*pot)->pot_sampler->extract_mean_of_means();

        } else {

            s << (*pot)->pot_sampler->extract_mean();

        }

        s << endl;

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
        potential->pot_sampler->getMeanErrorEstimator()->finalize();
        potential->pot_sampler->getMeanOfMeansErrorEstimator()->finalize();
    }
}
