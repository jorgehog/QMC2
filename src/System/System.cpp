/* 
 * File:   System.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 16:49
 */

#include "../QMCheaders.h"

System::System(int n_p, int dim, Orbitals* orbital) {
    this->n_p = n_p;
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

    for (std::vector<Potential*>::iterator pot = potentials.begin(); pot != potentials.end(); ++pot) {
        potE += (*pot)->get_pot_E(walker);
    }

    return potE;
}
