/* 
 * File:   System.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 16:49
 */

#include "../QMCheaders.h"

System::System(int n_p, int dim, Orbitals* orbital) {
    this->n_p = n_p;
    n2 = n_p / 2;
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
        (*pot)->pot_sampler.queue_value(potE_i);
        potE += potE_i;
    }

    return potE;
}

void System::update_potential_samples(double weight) {
    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        (*pot)->pot_sampler.update_mean(weight);
    }
}

void System::push_potential_samples() {
    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        (*pot)->pot_sampler.push_mean();
    }
}

std::string System::dump_samples(bool mean_of_means) {
    using namespace std;

    stringstream s;

    s << setprecision(6) << fixed;
    
    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        
        s << (*pot)->get_name() << " ";

        if (mean_of_means) {

            s << (*pot)->pot_sampler.extract_mean_of_means();

        } else {

            s << (*pot)->pot_sampler.extract_mean();

        }

        s << endl;

    }

    return s.str();

}