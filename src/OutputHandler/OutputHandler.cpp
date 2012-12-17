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
        int node,
        int n_nodes) {

    this->is_dmc = false;
    this->is_vmc = false;

    this->node = node;
    this->n_nodes = n_nodes;
    this->parallel = parallel;

    this->filename = filename;
    this->path = path;

    if (parallel) {
        filename = filename + boost::lexical_cast<std::string > (node);
    }

    //default:
    use_file = false;

}

void OutputHandler::init_file(){
    use_file = true;
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
        post_pointer_init();
    } else {
        this->min = min;
    }
}

void OutputHandler::finalize() {
    if (!use_file){
        return;
    }
    file.close();

    if (parallel & (node == 0)) {
        std::string compressorPath = "~/MASTER/QMC2/tools/compressData.sh";

        //Bash script "~/compressData ~/test/ blocking 4"
        int success = std::system(((((((((std::string)"bash " +
                                  compressorPath) +
                                  " ") +
                                  path) +
                                   " ") +
                                   filename) +
                                  " ") +
                                  boost::lexical_cast<std::string > (n_nodes)).c_str()
                                 );
        if (success != 0){
            std::cout << "compressData failed. " << std::endl;
            exit(1);
        }

    }
}