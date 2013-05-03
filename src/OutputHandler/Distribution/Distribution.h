/* 
 * File:   Distribution.h
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#ifndef DISTRIBUTION_H
#define	DISTRIBUTION_H

/*! \brief Class for calculating distribution functions such as the one-body density.
 *  Does not collect data itself, but works merely as a control organ for the QMC class,
 * calling it's methods for storing position data.
 */
class Distribution : public OutputHandler {
public:

    //! Constructor.
    /*!
     * @param path The path where output is stored (or read).
     * @param name Name of the file.
     */
    Distribution(ParParams &, std::string path, std::string name);

    /*!
     * Signals QMC solver to store position data.
     */
    void dump();

    //! Method for calculating the distribution.
    /*!
     * Overrides the superclass implementation.
     */
    void finalize();

    //! Method for re-calculating the distribution given a stores set of position data. 
    /*!
     * Scatters the data across nodes.
     * @param n_p Number of particles in the set.
     * @param N Number of mesh points used in the histogram.
     * @param bin_edge The Cartesian position of the end points of the histogram. 
     */
    void rerun(int n_p, int N, double bin_edge_xy, double bin_edge_xz = 0, double bin_edge_yz = 0);

private:

    int dim;

    std::string name;


    //! Method for generating the one-body density and projected one-axis distribution.
    /*! 
     * Calculations done by making histograms of stored position data.
     * If rerun is true, the position data is not written to file.
     * @param n_p Number of particles in the set.
     * @param N Number of mesh points used in the histogram.
     * @param bin_edge The Cartesian position of the end points of the histogram. 
     */
    void generate_distribution2D(arma::mat & dist,
            int n_p,
            double bin_edge = 0,
            int N = 200,
            bool rerun = false);

    void generate_distribution3D(arma::mat & dist,
            int n_p,
            double bin_edge = 0,
            int N = 200,
            bool rerun = false);

    /*!
     * In order to calculate the radial distribution, the dimension is needed.
     * \see OutputHandler::post_pointer_init()
     */
    void post_pointer_init() {
        this->dim = qmc->dim;
    }

};


#endif	/* DISTRIBUTION_H */
