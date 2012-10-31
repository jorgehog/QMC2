/* 
 * File:   Blocking.h
 * Author: jorgmeister
 *
 * Created on October 29, 2012, 3:58 PM
 */

#ifndef BLOCKING_H
#define	BLOCKING_H

class Blocking : public ErrorEstimator {
public:
    Blocking(int n_c, std::string filename = "blocking_out",
            std::string path = "./",
            bool parallel = false,
            int my_rank = 0,
            int num_procs = 1 );
    
    double estimate_error();
    
protected:
    
    arma::rowvec local_block;
    
    int min_block_size;
    int max_block_size;
    int n_block_samples;
    
    double block_data(int block_size);

};

#endif	/* BLOCKING_H */

