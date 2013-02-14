import re, os
from pyLibQMC import paths

mainCppFile = open(os.path.join(paths.CODE, "src", "QMCmain.cpp"), 'r')
raw = mainCppFile.read()
mainCppFile.close()

generateSTDOUTtests = True
generateCMLMaps = True


if generateSTDOUTtests:

    p = "\s*if \(def\.compare\(argv\[(\d+)\]\) != 0\) (.+) = "
    data = re.findall(p, raw)
    	
    print "    if (initOut) {"
    print "        if (parParams.is_master) {"
    
    for index, var in data:
        print '            std::cout << "%2s" << " %s" << " = " << %s << std::endl;' % (index, var.ljust(32), var.ljust(32))
	
    print "        }"
    print
    print "#ifdef MPI_ON"
    print "        MPI_Finalize();"
    print "#endif"
    print    
    print "        exit(0);\n    }\n\n\n"

if generateCMLMaps:

	p = "^(\s*if \(def\.compare\(argv\[\d+\]\) != 0\) .+ = .*argv\[\d+\].*\;\s*)$"
	data = re.findall(p, raw, re.MULTILINE)

	pyMap = {"ge" : [], "vm" : [], "dm" : [], "ou" : [], "mi" : [], "va" : []}

	cppMap = ""
	for i in range(len(data)):
		if "vmcParams.dt" in data[i]:
			vmcDt = i + 1
		elif "dmcParams.dist_in" in data[i]:
			dmcIn = i + 1
		elif "outputParams.dist_out" in data[i]:
			distOut = i + 1

		cppMap += re.sub("argv\[\d+\]", "argv[%d]" % (i+1), data[i]) + "\n"


		pyData = re.findall("\s*if \(def.compare\(argv\[\d+\]\) != 0\) (.+)\.(.+) = ", data[i])[0]
	
		pyMap[pyData[0][:2]].append((pyData[1], i))

		

	cppMap += "\tint vmc_dt_loc = %d;\n\tint dist_in_loc = %d;\n\tint dist_out_loc = %d;\n" % (vmcDt, dmcIn, distOut)
	
	l = None
	pyMapS = ""
	for key in ['ou', 'ge', 'vm', 'dm', 'mi', 'va']:
		if key == 'va':
			pyKey = 'vp'
		else:
			pyKey = key[0]
		
		start = "cmlMAP%s = {" % pyKey

		l = len(start)
		pyMapS += start
		
		pyMapS += '%s : %d,\n' % (('"%s"' % pyMap[key][0][0]).ljust(15), pyMap[key][0][1])
		for pair in pyMap[key][1:]:
			name, index = pair
			pyMapS += ''.ljust(l) + '%s : %d,\n' % (('"%s"' % name).ljust(15), index)
		pyMapS = pyMapS[:-2] + "}\n\n"

	print "[CPP]\n\n"
	print cppMap
	print "\n\n[PYTHON]\n\n"
	print pyMapS




		





