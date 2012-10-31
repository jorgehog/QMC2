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

    bool do_output;
    bool to_file;

    ErrorEstimator();

    ErrorEstimator(int n_c, std::string filename,
            std::string path,
            bool parallel,
            int my_rank, int num_procs, std::string* infile = NULL);

    virtual void update_data(double val) {
        data(i) = val;
        i++;
    }

    virtual double estimate_error() = 0;

protected:
    int n_c;
    int i;

    bool parallel;
    int my_rank;
    int num_procs;

    std::string filename;
    std::string path;
    std::ofstream file;

    arma::rowvec data;

    void finalize();

};

#endif	/* ERRORESTIMATOR_H */

