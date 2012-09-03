/* 
 * File:   BlockingData.cpp
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#include "../../QMCheaders.h"

BlockingData::BlockingData(std::string filename,
        std::string path,
        bool parallel,
        int my_rank,
        int num_procs)
: OutputHandler(filename, path, parallel, my_rank, num_procs) {

}

void BlockingData::dump() {
    file << qmc->local_E << endl;
}