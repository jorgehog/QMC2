import tarfile
import os, sys, re, inspect
from os.path import join as pjoin
from subprocess import Popen, PIPE, call

from pyLibQMC import paths



def update():
    p = Popen("ifconfig", stdout=PIPE)
    
    ip_eth0, ip_lo = re.findall("inet addr\:(.+?)\s", p.stdout.read())
    
    thisFilePath = os.path.abspath(inspect.getfile(inspect.currentframe()))
    
    with open(thisFilePath, 'r') as thisFile:
        new = re.sub("(ssh_ip\s*=\s*).*?\n", "\g<1>\"%s\"\n" % ip_eth0, thisFile.read())
    thisFile.close()

    with open(thisFilePath, 'w') as thisFile:
        thisFile.write(new)
    thisFile.close()
    
    print "ip successfully updated."


def excludeFunc(filename):
    
    filesizeLimit = 10 #MB
    filesizeLimit *= 1024*1024 #convert to MB
    
    excludes = []
    
    size = os.path.getsize(filename)
    file_ = os.path.split(filename)[1]                
                
    notExcluded = True
    for key in excludes:
        notExcluded = notExcluded&(not key in file_)
 
    return not(0 < size < filesizeLimit and notExcluded)
                   

def packFiles(runPath):

    tarPath = pjoin(paths.scratchPath, "tarfiles")

    if not os.path.exists(tarPath):
        print "initialized dir: ", tarPath
        os.mkdir(tarPath)
        
    
    if runPath[-1] == "/":
        runPath = runPath[:-1]
    arcname = os.path.split(runPath)[1] 
    
    
    thisTarpath = pjoin(tarPath, arcname  + ".tar.gz")
    print "Path: ", thisTarpath
    
    
    tar = tarfile.open(thisTarpath, 'w:gz')
    
    for root, dirs, files in os.walk(runPath):
        
        print "Entring directory: ", root
        for file_ in files:
            abspath = pjoin(root, file_)
            print "packing ", file_
            tar.add(abspath, exclude=excludeFunc)                
    
    tar.close()
    return thisTarpath       

def sendTar(path):    
    
    ssh_ip = "193.157.210.112"
 
    base, name = os.path.split(path)
    recv = "jorgmeister@%s:/home/jorgmeister/recv/%s" % (ssh_ip, name)

    arg = " ".join(["scp", path, recv])
    
    call(arg, shell=True)


def main():

    if len(sys.argv) < 2:
        print "insufficient commandline arguments."
        print "python packFiles.py update | [runpath]" 
    
    if sys.argv[1] == "update":
        update()
    else:
        tarPath = packFiles(sys.argv[1])
        sendTar(tarPath)    

if __name__ == "__main__":
    main()
    

#echo scp $rundir $recvID@$ip:$recvDir