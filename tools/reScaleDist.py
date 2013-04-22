# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:53:26 2013

@author: jorgehog
"""

import sys, re, subprocess, os, shutil
from os.path import join as pjoin
from pyLibQMC import parseCML, paths, misc


def recalcVMCdist(VMC_raw_path, VMC_raw_name, n_p, bin_edge, N, n_cores, MPI_flag):
    
    exe = pjoin(paths.programPath, misc.QMC2programName)    
    
    path = VMC_raw_path.split("walker_positions")[0]
    name = re.findall("dist_rawdata_(.+?)\.arma", VMC_raw_name)[0]  
    
    MPI = ""
    if MPI_flag:
        MPI = "mpiexec -n %d" % n_cores
    
    args = MPI.split() + [exe, "redist", n_p, path, name, N, bin_edge]

    subprocess.call(args)

def main():
        
    stdoutToFile, MPI_flag, openGUI, n_cores = parseCML(sys.argv)
    print "MPI Nodes: ", n_cores

    if len(sys.argv) < 3:
        sys.exit("Sufficinet cml args not supplied. (VMC rawdata, DMC dist file)")
    if not "vmc" in sys.argv[1] or not "dmc" in sys.argv[2]:
        sys.exit("Incorrect input type. (VMC rawdata, DMC dist file)")
    
    
    DMC_dist = sys.argv[2]
    VMC_raw_path, VMC_raw_name = os.path.split(sys.argv[1])
    DMC_dist_path, DMC_dist_name = os.path.split(DMC_dist)
    
    if len(sys.argv) > 3:
        N = sys.argv[3]
    else:
        N = "200"
    
    VMC_name = re.findall("dist_rawdata_(.+?)vmc\.arma", VMC_raw_name)[0]
    DMC_name = re.findall("dist_out_(.+?)dmc_", DMC_dist_name)[0]
    
    if VMC_name != DMC_name:
        raise Exception("Distributions does not match. (%s != %s)" % (VMC_name, DMC_name))
    
    n_p = re.findall("(\d+)c\d+", VMC_name)[0]
    bin_edge = re.findall("dist_out.+?dmc_edge(\d+\.?\d+)\.arma", DMC_dist_name)[0]
    
#    print n_p
#    print bin_edge
#    print N    
#    print VMC_raw_path
#    print VMC_raw_name
#    print "-----"
   
    
    recalcVMCdist(VMC_raw_path, VMC_raw_name, n_p, bin_edge, N, n_cores, MPI_flag)
    
    new_dir = pjoin("OneBodyDensities", DMC_name)
    new_path = pjoin(paths.scratchPath, new_dir)

    new_VMC_dist = "dist_out_%svmc_edge%s.arma" % (VMC_name, bin_edge)
    new_VMC_dist_rad = "radial_out_%svmc_edge%s.arma" % (VMC_name, bin_edge)
    
    DMC_dist_rad = "radial_out_%sdmc_edge%s.arma" % (DMC_name, bin_edge)
    
#    print new_path
#    print new_VMC_dist
#    print new_VMC_dist_rad
#    print DMC_dist    
#    print DMC_dist_rad
    
    if not os.path.exists(new_path):
        os.mkdir(new_path)
        
    shutil.copy(pjoin(VMC_raw_path, new_VMC_dist), new_path)
    shutil.copy(pjoin(VMC_raw_path, new_VMC_dist_rad), new_path)
    shutil.copy(pjoin(DMC_dist_path, DMC_dist_rad), new_path)
    shutil.copy(DMC_dist, new_path)    
    

    
    
if __name__ == "__main__":
    main()