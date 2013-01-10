### DCViz (with lazers) ###

... is a library developed by the wizard jorgen, written in horrible Python.
It is designed to work, not be pretty. And speaking of which, this is how it works:

-Plots column data in matplotlib either
--real time: Given a delay interval [0,10], data is replotted in the same windows
             after the given time.
--not real time: Data is plotted!

### Example
    An example is located in the Example folder! Simply run 'Example.py' to generate 3 testcases and load the GUI.
    Press 'ctrl+s' and open the folder that pops up (the example folder). Choose a case from the dropdown and press play.

### Usage

First of all, check the example

-I have a file containing columndata called 'testcase2.dat', and I want to plot it!

Howto: create a subclass 'DCvizPlotter', with nametag set as 'testcase\d\.dat' (regex) and it will be recognized, so:

       def myNewClass(dcv_sup):
           nametag = 'testcase\d\.dat'

       set a 'figMap', a python dictionary, on how you want figures and subfigures set up, e.g.
       
       def myNewClass(dcv_sup):
           nametag = 'OUT\d\.dat'
           figMap = {'fig1': ['subfig1', 'subfig2'], 'fig2': ['subfig3'] ....}
           
       
       you will then have access to 'self.fig1', self.subfig1 etc. when you implement the virtual 
       plot function:
       
       def myNewClass(dcv_sup):
           nametag = 'testcase\d\.dat'
           figureMap = {'fig1': ['subfig1', 'subfig2'], 'fig2': ['subfig3'] ....}
           
           
           def plot(self):
               
               fig1, subfig1 = self.fig1, self.subfig1
               
               #Do all the plotting in the world to make your plots look nice.
               
        This function should work without a GUI, so simply test it by running
        
        c = myNewClass(somePathToTheFile, dynamic = True/False) #OUT2.dat file
        c.mainloop() #If dynamic is set true, you interrupt the loop with ctrl+c
        
        
        
        
Ok I have done so, now I want to use the GUI!


    Allright, you can start the DCVizGUI.main(masterDir) function with a chosen masterDir
    (the place all your open dialoges windows will 'look first'. Typically you run folder)
    
    hotkeys (also aviable from file/options menu):
    -ctrl+s : opens a dialogue window to set the path (and load files in this path)
              File detection is based on the nametag. It will be matched with a subclass of dcv_sup
              (see prev section)
   
    -ctrl+o : opens a single file as above
   
    -ctrl+r : reloads the config file
    
    -ctrl+c : display the config file
    
     The config (located in the GUI folder):
     
     you can set the delay time (also setable from the slider in the GUI)
     you can silence warnings, and also all terminal output if you want to
     
     
Ok I loaded some files, what now?

    The dropdown menu should be selectable; choose your file and press play!
    The rest should be intuitive!
    

Cheers!
    
        
        
       

