# -*- coding: utf-8 -*-

import sys, re, subprocess
from os.path import join as pjoin
from pyLibQMC import parseCML, paths, misc

sys.path.append(pjoin(paths.toolsPath, "DCViz", "src"))
    
try:
    import DCVizWrapper as viz
    canViz = True
except:
    canViz = False


def initRun(n_p, path, name, N, bin_edge, n_cores, mpiFlag):
    exe = pjoin(paths.programPath, misc.QMC2programName)
    
    mpi = ""
    if mpiFlag:
        mpi = "mpiexec -n %d" % n_cores
    
    args = mpi.split() + [exe, "redist", n_p, path, name, N, bin_edge]

    subprocess.call(args)
    
    outPath = "%swalker_positions/__type___out_%s_edge%s.arma" % (path, name, bin_edge)
    dist_path = outPath.replace("__type__", "dist")
    radial_path = outPath.replace("__type__", "radial")
    
    if canViz:
        viz.main(dist_path, False)
        viz.main(radial_path, False)


def main():
    
    stdoutToFile, mpiFlag, openGUI, n_cores = parseCML(sys.argv)
    print "MPI Nodes: ", n_cores
    
    #path N binEdge
    if len(sys.argv) < 3:
        print "Error in command line arguments. [rawdata, N, bin egde]"
        sys.exit(1)    
        
    rawfile = sys.argv[1]

    if 'rawdata' not in rawfile:
        print "Selected file is not a distribution data file."
        sys.exit(1)
    
    N = sys.argv[2]
    bin_edge = sys.argv[3]
    
    name = re.findall("dist_rawdata_(.+).arma", rawfile)[0]
    n_p = re.findall("(\d+)c\d", name)[0]
    path = rawfile.split("walker_positions")[0]

    initRun(n_p, path, name, N, bin_edge, n_cores, mpiFlag)

if __name__ == "__main__":
    main()

