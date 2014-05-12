#pragma once


#include <string>
#include <sstream>
#include <fstream>


namespace QMC2
{


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
    
    bool parallel;
    int node;
    int n_nodes;

    bool use_file; //!< If init_file() is called, this flag is set true. Assures correct finalization.

    std::stringstream s;

    std::string filename;
    std::string path;

    std::ofstream file;

    /*!
     * Opens a file with filename at path supplied in constructor.
     * Subclass implementations can call this function. Superclass does not.
     */
    void init_file();

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

    virtual void reset() {
        if (use_file) init_file();
    }

};

}
