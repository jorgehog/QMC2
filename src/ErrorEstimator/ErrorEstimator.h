/* 
 * File:   ErrorEstimator.h
 * Author: jorgmeister
 *
 * Created on October 29, 2012, 3:58 PM
 */

#ifndef ERRORESTIMATOR_H
#define	ERRORESTIMATOR_H

#include <string>
#include <fstream>
#include <armadillo>

/*! \brief Class handling error estimations of the QMC methods.
 * The QMC class holds an object of this type, calling the update_data
 * function in order to update the sampling pool.
 * finalize() then either dumps the samples to file for later processing,
 * or calulates an estimate.
 */
class ErrorEstimator {
public:

    bool data_to_file; //!< If true, the data vector are stored to file.
    bool output_to_file; //!< If init_file() method is called, this flag is true.
    
    ErrorEstimator();

    //! Constructor.
    /*!
     * @param n_c The expected number of samples to be stored.
     * @param filename The name of the file. Only necessary if init_file() is called.
     * @param path The path where the data and/or the file is stored (or read).
     */
    ErrorEstimator(int n_c, std::string filename,
            std::string path,
            bool parallel,
            int node, int n_nodes, bool rerun = false);

    //! Calculates the combined mean of n_nodes means.
    /*!
     * Only useful for parallel calls.
     * @param mean The local mean on an individual node
     * @param n The number of samples used to calculate the local mean.
     * @param n_tot The total number of samples on all nodes. Calculated if not supplied.
     */
    static double combine_mean(double mean, int n, int n_tot = 0);
    
    //! Calculates the combined variance of n_nodes variances.
    /*!
     * Only useful for parallel calls.
     * @param var The local variance on an individual node
     * @param mean The local mean on an individual node
     * @param n_tot The total number of samples used on all nodes.
     */
    double combine_variance(double var, double mean = 0, int n = 0);
    
    /*!
     * if [output_to_file]: Closes opened files
     * if [data_to_file]: Stores accumulated data.
     * if data vector was used, it's memory is freed.
     */
    void finalize();

    /*!
     * Gathers the data vectors from all processes into a single one on the master node.
     */
    void node_comm_gather_data();
    
    /*!
     * Exact reverse of node_comm_gather_data()
     */
    void node_comm_scatter_data();
    
    //! Opens a file with filename at path supplied in constructor.
    /*!
     * Subclass implementations can call this function. Superclass does not.
     */
    void init_file();

    //!Estimates the error based on the subclass implementation.
    virtual double estimate_error() = 0;
    
//    virtual void normalize();

    //! Adds values to the data vector.
    /*!
     * Can be overridden if storage is not wanted.
     * @param val A local sample of the quantity of which the error is calculated
     */
    virtual void update_data(double val); 


protected:
    int n_c; //!< Size of the data vector
    int i; //!< Count variable for the data vector.

    bool parallel;
    bool is_master;
    int node;
    int n_nodes;

    bool rerun; //!< If false, data is assumed to already exist.

    std::string filename;
    std::string path;
    std::ofstream file;

    arma::rowvec data; //!< The vector containing the samples used in error calculation.

};

#endif	/* ERRORESTIMATOR_H */

