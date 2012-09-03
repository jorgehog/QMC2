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
        int num_procs) {

    this->my_rank = my_rank;
    this->num_procs = num_procs;
    this->parallel = parallel;

    this->filename = filename;
    this->path = path;

    if (parallel) {
        filename = filename + boost::lexical_cast<std::string > (my_rank);
    }

    this->file.open(((path + filename) + ".dat").c_str());

}

void OutputHandler::finalize() {
    file.close();

    if (parallel) {
        if (my_rank == 0) {
            std::string compressorPath = "~/MASTER/QMC2/tools/compressData.sh";

            //Bash script "~/compressData ~/test/ blocking 4"
            std::system(((((((((std::string)"bash " +
                    compressorPath) +
                    " ") +
                    path) +
                    " ") +
                    filename) +
                    " ") +
                    boost::lexical_cast<std::string > (num_procs)).c_str()
                    );

        }
    }
}