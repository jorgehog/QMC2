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

    bool data_to_file;
    bool output_to_file;

    ErrorEstimator();

    ErrorEstimator(int n_c, std::string filename,
            std::string path,
            bool parallel,
            int node, int n_nodes, bool rerun = false);

    void finalize();

    void node_comm();
    
    void init_file();

    virtual double estimate_error() = 0;

    virtual void update_data(double val) {
        data(i) = val;
        i++;
    }
    
    void clear(){
        data.clear();
    }

protected:
    int n_c;
    int i;

    bool parallel;
    int node;
    int n_nodes;

    bool rerun;

    std::string filename;
    std::string path;
    std::ofstream file;

    arma::rowvec data;

};

#endif	/* ERRORESTIMATOR_H */

