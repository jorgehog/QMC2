/* 
 * File:   ErrorEstimator.h
 * Author: jorgmeister
 *
 * Created on October 29, 2012, 3:58 PM
 */

#ifndef ERRORESTIMATOR_H
#define	ERRORESTIMATOR_H

class ErrorEstimator {
public:
    
    ErrorEstimator();
    
    ErrorEstimator(int n_c, std::string filename,
            std::string path,
            bool parallel,
            int my_rank, int num_procs);

    virtual void update_data(double val, int i) {
        data(i) = val;
    }

    virtual double estimate_error() = 0;

    void set_n_c(int n_c) {
        this->n_c = n_c;
    }


protected:
    int n_c;

    bool parallel;
    int my_rank;
    int num_procs;

    std::string filename;
    std::string path;
    std::ofstream file;

    arma::rowvec data;
    
    void finalize(){
        file.close();
    }

};

#endif	/* ERRORESTIMATOR_H */

