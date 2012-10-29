/* 
 * File:   Blocking.cpp
 * Author: jorgmeister
 * 
 * Created on October 29, 2012, 3:58 PM
 */

#include "../../QMCheaders.h"

Blocking::Blocking(int n_c,
        std::string filename,
        std::string path,
        bool parallel,
        int my_rank,
        int num_procs)
: ErrorEstimator(n_c, filename, path, parallel, my_rank, num_procs) {

}

//double mean(double *vals, int n_vals) {
//
//    double m = 0;
//    for (int i = 0; i < n_vals; i++) {
//        m += vals[i];
//    }
//
//    return m / double(n_vals);
//}
//
//// calculate mean and variance of vals, results stored in res
//
//void meanvar(double *vals, int n_vals, double *res) {
//    double m2 = 0, m = 0, val;
//    for (int i = 0; i < n_vals; i++) {
//        val = vals[i];
//        m += val;
//        m2 += val*val;
//    }
//
//    m /= double(n_vals);
//    m2 /= double(n_vals);
//
//    res[0] = m;
//    res[1] = m2 - (m * m);
//
//}
//
//// find mean and variance of blocks of size block_size.
//// mean and variance are stored in res
//
//void blocking(double *vals, int n_vals, int block_size, double *res) {
//
//    // note: integer division will waste some values
//    int n_blocks = n_vals / block_size;
//
//    /*
//    cerr << "n_vals=" << n_vals << ", block_size=" << block_size << endl;
//    if(n_vals%block_size > 0)
//      cerr << "lost " << n_vals%block_size << " values due to integer division" 
//           << endl;
//     */
//
//    double* block_vals = new double[n_blocks];
//
//    for (int i = 0; i < n_blocks; i++) {
//        block_vals[i] = mean(vals + i*block_size, block_size);
//    }
//
//    meanvar(block_vals, n_blocks, res);
//
//    delete block_vals;
//}
//
//int main(int nargs, char* args[]) {
//
//    int n_procs, min_block_size, max_block_size, n_block_samples;
//    //   Read from screen a possible new vaue of n
//    if (nargs > 4) {
//        n_procs = atoi(args[1]);
//        min_block_size = atoi(args[3]);
//        max_block_size = atoi(args[4]);
//        n_block_samples = atoi(args[5]);
//    } else {
//        cerr << "usage: ./mcint_blocking.x <n_procs> <min_bloc_size> "
//                << "<max_block_size> <n_block_samples>" << endl;
//        exit(1);
//    }
//
//
//    int local_n, n;
//
//    n = pow(10., atof(argv[2]);
//            local_n = n / 4;
//
//            // get all mc results from files
//            double* mc_results = new double[n];
//
//
//            string line;
//            ifstream myfile;
//            myfile.open("energies.dat", ios::in);
//
//            double lol = 0;
//
//
//    for (int i = 0; i < n; i++) {
//        getline(myfile, line);
//                mc_results[i] = atof(line.c_str());
//    }
//    myfile.close();
//
//
//
//
//
//            // and summarize
//            double mean, sigma;
//            double res[2];
//
//            meanvar(mc_results, n, res);
//
//            mean = res[0]; sigma = res[1];
//
//            cout << "Value of integral = " << mean
//            << " Analytic result =  " << 3 << endl;
//            cout << "Value of variance = " << sigma
//            << endl; ///n-(total_sum*total_sum/(n*n)) << endl;
//            cout << "Standard deviation = "
//            << sqrt(sigma / (n - 1.0))
//            << endl;
//
//
//            // Open file for writing, writing results in formated output for plotting:
//            ofstream outfile;
//            outfile.open("blockres.dat", ios::out);
//
//            outfile << setprecision(10);
//
//            double* block_results = new double[n_block_samples];
//
//            int block_size, block_step_length;
//
//            block_step_length = (max_block_size - min_block_size) / n_block_samples;
//
//            // loop over block sizes
//    for (int i = 0; i < n_block_samples; i++) {
//        block_size = min_block_size + i*block_step_length;
//                blocking(mc_results, n, block_size, res);
//
//                mean = res[0];
//                sigma = res[1];
//
//                // formated output
//                outfile << block_size << "\t" << sqrt(sigma / ((n / block_size) - 1.0))
//                << endl;
//
//    }
//
//    outfile.close();
//
//    return 0;
//}
