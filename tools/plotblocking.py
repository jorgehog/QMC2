from scitools.std import plot, figure
import re, os

allfiles = os.listdir(os.getcwd())
dirs = []
for f in allfiles:
	if os.path.isdir(f):
		dirs.append(f)

i = 1
for dirname in dirs:
	os.chdir(dirname)

	files = " ".join(os.listdir(os.getcwd()))
	files = re.findall("blocking_\w+?_out\d?.dat", files)

	print files
	
	for bfile in files:
		f = open(bfile, 'r')

		blocks = []
		error = []

		for line in f:
			blocks.append(int(line.split()[0]))
			error.append(float(line.split()[1]))
		f.close()
		
		if blocks and error:
			figure(i)
			title = "blocking " + bfile.split("_")[1] + str(re.findall(".*(\d+).dat",bfile.split("_")[2])) + " " + dirname
			plot(blocks, error, "k*", title=title.replace("['0']", " alpha").replace("['1']", " beta").replace("[]", ""), xlabel="block size", ylabel="error")
		else:
			print "no entries in file: ", dirname + bfile, " -- skipping." 
	
		i+=1

	os.chdir("../")

raw_input()
