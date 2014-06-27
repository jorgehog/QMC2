#pragma once


#include "../ErrorEstimator.h"


namespace QMC2
{

struct ParParams;

class Blocking : public ErrorEstimator {
public:

    Blocking(int n_c, ParParams & pp,
            std::string filename = "blocking_out",
            std::string path = "",
            int n_b = 100,
            int maxb = 10000,
            int minb = 10,
            bool rerun = false);

    double estimate_error();
    
    /*!
     * Calculates the variance as in SimpleVar
     */
    void get_initial_error();
    
    //! Calculates the block sizes.
    /*!
     * Due to integer division, alot of sizes
     * becomes equal. Only unique block sizes are returned.
     * @param block_sizes Vector containing the block sizes
     * @param n The number of unqiue block sizes
     */
    void get_unique_blocks(arma::Row<int> & block_sizes, int & n);

protected:

    arma::rowvec local_block;

    int min_block_size; //!< The minimum amount of samples in one block
    int max_block_size; //!< The maximum amount of samples in one block
    int n_block_samples; //!< The total amount of different block sizes.

    //! Calculates the variance and mean of the dataset with the specified block size.
    /*!
     * @param block_size The number of samples in each block.
     * @param var Reference to the variance of the block's means.
     * @param mean Reference to the mean of the block's means. Needed to combine the variances from different processes.
     */
    void block_data(int block_size, double &var, double &mean);

};

}
