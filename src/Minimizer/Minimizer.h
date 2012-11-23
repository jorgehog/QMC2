/* 
 * File:   MINIMIZER.h
 * Author: jorgehog
 *
 * Created on 23. August 2012, 16:52
 */

#ifndef MINIMIZER_H
#define	MINIMIZER_H

class Minimizer {
protected:
    
    int n_nodes;
    bool is_master;

    VMC* vmc;

    STDOUT* std_out;
    std::stringstream s;
    
    int Nspatial_params;
    int Njastrow_params;
    int Nparams;

    std::vector<OutputHandler*> output_handler;
    std::vector<ErrorEstimator*> error_estimators;

    void dump_output();
    void finalize_output();
    void error_output();
    
    virtual void node_comm() = 0;


    //virtual void get_variational_gradients(); is this a general req?

public:

    Minimizer(VMC* vmc, const ParParams &, const arma::rowvec & alpha, const arma::rowvec & beta);

    void add_output(OutputHandler* output_handler);

    Orbitals* get_orbitals() {
        return vmc->get_orbitals_ptr();
    }

    Jastrow* get_jastrow() {
        return vmc->get_jastrow_ptr();
    }

    virtual VMC* minimize() = 0;

    void output(std::string message, double number = -1);

    void add_error_estimator(ErrorEstimator* error_estimator) {
        this->error_estimators.push_back(error_estimator);
    }
};


#endif	/* MINIMIZER_H */