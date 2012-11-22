/* 
 * File:   OutputHandler.cpp
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#include "../QMCheaders.h"

OutputHandler::OutputHandler() {

}

OutputHandler::OutputHandler(std::string filename,
        std::string path,
        bool parallel,
        int my_rank,
        int n_nodes) {

    this->is_dmc = false;
    this->is_vmc = false;

    this->my_rank = my_rank;
    this->n_nodes = n_nodes;
    this->parallel = parallel;

    this->filename = filename;
    this->path = path;

    if (parallel) {
        filename = filename + boost::lexical_cast<std::string > (my_rank);
    }

    this->file.open(((path + filename) + ".dat").c_str());

}

void OutputHandler::set_qmc_ptr(QMC* qmc) {
    if (this->is_dmc) {
        dmc = (DMC*) qmc;
    } else if (this->is_vmc) {
        vmc = (VMC*) qmc;
    } else {
        this->qmc = qmc;
    }
}

void OutputHandler::set_min_ptr(Minimizer* min) {
    if (this->is_ASGD) {
        this->asgd = (ASGD*) min;
    } else {
        this->min = min;
    }
}

void OutputHandler::finalize() {
    file.close();

    if (parallel & (my_rank == 0)) {
        std::string compressorPath = "~/MASTER/QMC2/tools/compressData.sh";

        //Bash script "~/compressData ~/test/ blocking 4"
        std::system(((((((((std::string)"bash " +
                compressorPath) +
                " ") +
                path) +
                " ") +
                filename) +
                " ") +
                boost::lexical_cast<std::string > (n_nodes)).c_str()
                );

    }
}