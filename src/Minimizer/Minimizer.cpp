/* 
 * File:   Minimizer.cpp
 * Author: jorgehog
 * 
 * Created on 23. august 2012, 16:52
 */

#include "../QMCheaders.h"

Minimizer::Minimizer(VMC* vmc, const ParParams & pp, const arma::rowvec & alpha, const arma::rowvec & beta) {

    this->vmc = vmc;

    this->Nspatial_params = alpha.n_rows;
    this->Njastrow_params = beta.n_rows * vmc->get_jastrow_ptr()->active;
    this->Nparams = this->Njastrow_params + this->Nspatial_params;

    for (int i = 0; i < Nspatial_params; i++) {
        vmc->get_orbitals_ptr()->set_parameter(alpha(i), i);
    }

    for (int i = 0; i < Njastrow_params; i++) {
        vmc->get_jastrow_ptr()->set_parameter(beta(i), i);
    }
    
    if (pp.is_master){
        std_out = new STDOUT();
    } else {
        std_out = new NO_STDOUT();
    }
    
    is_master = pp.is_master;
    n_nodes = pp.n_nodes;

}

void Minimizer::output(std::string message, double number) {
    using namespace std;

    if (number != -1) {
        s << message << " " << number << endl;
    } else {
        s << message << endl;
    }


    if (Nspatial_params != 0) s << "\nAlpha:\n";
    for (int alpha = 0; alpha < Nspatial_params; alpha++) {
        s << "\t" << vmc->get_orbitals_ptr()->get_parameter(alpha) << endl;
    }

    if (Njastrow_params != 0) s << "\nBeta:\n";
    for (int beta = 0; beta < Njastrow_params; beta++) {
        s << "\t" << vmc->get_jastrow_ptr()->get_parameter(beta) << endl;
    }

    s << endl;
    std_out->cout(s);
    
}

void Minimizer::add_output(OutputHandler* output_handler) {
    output_handler->set_min_ptr(this);
    this->output_handler.push_back(output_handler);
}

void Minimizer::dump_output() {

    for (std::vector<OutputHandler*>::iterator output_obj = output_handler.begin(); output_obj != output_handler.end(); ++output_obj) {
        (*output_obj)->dump();
    }

}

void Minimizer::finalize_output() {

    for (std::vector<OutputHandler*>::iterator output_obj = output_handler.begin(); output_obj != output_handler.end(); ++output_obj) {
        (*output_obj)->finalize();
    }

}

void Minimizer::error_output() {
    if (error_estimators.at(0)->do_output) {
        for (int i = 0; i < Nparams; i++) {
            double error = error_estimators.at(i)->estimate_error();
            std::cout << "Error" << i << ": " << error << std::endl;
        }
    }
}