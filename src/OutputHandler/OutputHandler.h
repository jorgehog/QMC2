/* 
 * File:   OutputHandler.h
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#ifndef OUTPUTHANDLER_H
#define	OUTPUTHANDLER_H

#include <string>
#include <sstream>
#include <fstream>

class QMC;
class DMC;
class VMC;
class Minimizer;
class ASGD;

/*! \brief Class for handling output-methods.
 * Designed to avoid rewriting code, as well as avoid if-tests if output is
 * not desired.
 * \see QMC::output_handler, Minimizer::output_handler
 */
class OutputHandler {
protected:
    
    bool is_vmc;  //!< Switch used to typecast the QMC object to a VMC object.
    bool is_dmc;  //!< Switch used to typecast the QMC object to a DMC object.
    bool is_ASGD; //!< Switch used to typecast the Min object to an ASGD object.

    bool parallel;
    int node;
    int n_nodes;

    bool use_file; //!< If init_file() is called, this flag is set true. Assures correct finalization.

    std::stringstream s;

    std::string filename;
    std::string path;

    std::ofstream file;

    QMC* qmc;
    DMC* dmc;
    VMC* vmc;
    Minimizer* min;
    ASGD* asgd;

    /*!
     * Opens a file with filename at path supplied in constructor.
     * Subclass implementations can call this function. Superclass does not.
     */
    void init_file();

    //! Method for initialization requires once the correct QMC/Min pointer type is set.     
    /*!
     * Defaults to nothing.
     */
    virtual void post_pointer_init() {
    };

public:

    OutputHandler();
    
    //! Constructor.
    /*!
     * @param filename The name of the output file.
     * @param path The path of the output.
     */
    OutputHandler(std::string filename,
            std::string path,
            bool parallel,
            int node,
            int n_nodes
            );

    //!Methods for updating the output. 
    /*!
     * Typically retrieves information
     * through the solver pointers (given correct accessibility levels/friend)
     */
    virtual void dump() = 0;
    
    //!Finalizes the output. 
    /*!
     * Closes file if use_file flag is true.
     * Can be overridden if more complex tasks needs to be done, such as calculating
     * histograms etc.
     * \see Distribution::finalize()
     */
    virtual void finalize();

    /*!
     * Sets the QMC pointer and typecasts it according to the solver flags.
     */
    void set_qmc_ptr(QMC* qmc);
    
    /*!
     * Sets the Min pointer and typecasts it according to the minimizer flags.
     */
    void set_min_ptr(Minimizer* min);

};


#endif	/* OUTPUTHANDLER_H */