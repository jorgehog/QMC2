#include "OutputHandler.h"
#include "../defines.h"

using namespace QMC2;

OutputHandler::OutputHandler() {

}

OutputHandler::OutputHandler(std::string filename,
                             std::string path,
                             bool parallel,
                             int node,
                             int n_nodes) {

    this->node = node;
    this->n_nodes = n_nodes;
    this->parallel = parallel;

    this->filename = filename;
    this->path = path;

    if (parallel) {
        filename = filename + TOSTR(node);
    }

    //default:
    use_file = false;

}

void OutputHandler::init_file() {

    using std::stringstream;

    use_file = true;

    stringstream s;

    s << path << filename;

    if (!suffix.empty()) {
        s << suffix;
    }

    s << ".dat";

    this->file.open(s.str().c_str());
}

void OutputHandler::finalize() {
    if (!use_file) {
        return;
    }
    file.close();
}

std::string OutputHandler::suffix;
