# -*- coding: utf-8 -*-

from Tkinter import *
import tkMessageBox
import tkFileDialog

import sys, os, re, threading, time

from classfile import *

class spawned_thread(threading.Thread):
    def __init__(self, mode):
        threading.Thread.__init__(self)
        
        self.mode = mode
        self._stop = threading.Event()
        
    def run(self):
        self.mode.mainloop(self._stop)
            
    def stop(self):
        self._stop.set()

    def stopped(self):
        return self._stop.isSet()
        



class QMCGUI:

    def __init__(self, root):
        
        self.root = root
        self.path = os.getcwd()
        self.load_images()
        
        self.active_mode = None
        self.job = None
        
        #Initial terminal output flag
        self.hide_source = False
        self.terminal_silence = False
        self.warning_silence = False

        self.load_ext()
        
        #~ Menu :::::::::::::::::::::::::::::::::::::::::::
        menu = Menu(root)
        root.config(menu=menu)

        filemenu = Menu(root)
        optmenu = Menu(root)
        menu.add_cascade(label="File", menu=filemenu)
        menu.add_cascade(label="Options", menu=optmenu)

        #File menu
        filemenu.add_command(label="Open...", command=self.openfile)
        filemenu.add_command(label="Open path...", command=self.setpath)
        filemenu.add_separator()
        filemenu.add_command(label="Clear data", command=self.clean_modeselector)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=root.destroy)
        
        #Options menu
        optmenu.add_command(label="Display config", command=lambda : self.showInfo("Config", self.config))
        optmenu.add_command(label="Reload config", command=self.load_ext)
        #::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #~ StartButton ::::::::::::::::::::::::::::::::::::
        self.start_button = Button(root, image=self.img["play"], command=self.start, relief=FLAT, height=50, width=50)
        #::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #~ DynamicCheckBox ::::::::::::::::::::::::::::::::
        self.dynamic = BooleanVar()
        self.dynamicCheckBox = Checkbutton(root, text="Dynamic", variable=self.dynamic, command=self.change_dynamic_state)
        #::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #~ Mode selection dropdown menu :::::::::::::::::::
        self.mode = StringVar()
        self.mode.set("No data...")
        self.modeMap = {"" : "No data..."}
        self.modeSelector = OptionMenu(root, self.mode, *self.modeMap.values())
        self.modeSelector.config()
        self.modeMap = {}
        #::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #~ Refresh timer ::::::::::::::::::::::::::::::::::
        self.timer_text = Label(root, text="Refresh dt [s]:")
        self.refresh_dt = StringVar(value=str(self.refresh_dt_config))
        self.timer_field = Entry(root, textvariable=self.refresh_dt, width=10)
        self.timer_field.bind("<Return>", self.update_dt)
        #::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #~ grid Setup :::::::::::::::::::::::::::::::::::::
        self.start_button.grid(row=0, column=0, sticky=E + W + N + S)
        self.dynamicCheckBox.grid(row=1,column=0, sticky=E)
        self.modeSelector.grid(row=0, column=1, columnspan=3, sticky = E + W)
        self.timer_text.grid(row=2, column=0, sticky=W)
        self.timer_field.grid(row=2, column=2, sticky = W)
        #::::::::::::::::::::::::::::::::::::::::::::::::::
        
        #~ Initial states :::::::::::::::::::::::::::::::::
        self.start_button['state'] = DISABLED
        self.modeSelector['state'] = DISABLED
        self.timer_text['state'] = DISABLED
        self.timer_field['state'] = DISABLED
        #::::::::::::::::::::::::::::::::::::::::::::::::::
        
    def start(self):
            
        try:
            self.active_mode = self.modeMap[self.mode.get()]
        except:
            self.raiseWarning("Select a data set.")
            return
        
        self.active_mode.dynamic = self.dynamic.get()
        if self.active_mode.dynamic:
            self.start_button.configure(image=self.img["stop"], command=self.stop, relief=FLAT)
            
        self.job = spawned_thread(self.active_mode)
        self.job.setDaemon(True)
            
        if not self.job.is_alive():
            self.terminal_tracker("Job", "Starting.")
            self.job.start()
        else:
            print "BUG THREAD: This should never happen..."
    
        
    def stop(self):
        
        if self.job.is_alive():

            self.terminal_tracker("Job", "Stoping... ", hold="on")
            self.job.stop()
            self.job.join()
            
            self.terminal_tracker("Job", "Stopped.")

        else:
            print "BUG THREAD: This should never happen..."
            
        self.start_button.configure(image=self.img["play"], command=self.start, state=ACTIVE, relief=FLAT)


    def openfile(self):
        
        filename = tkFileDialog.askopenfilename(title="Choose file to display")

        if not filename:
            return
        
        self.detect_modetype(filename)
        self.update_modeselector()
        
    def setpath(self):

        self.path = tkFileDialog.askdirectory(title="Choose path ...")
    
        if not self.path:
            return
        
        for content in os.listdir(self.path):
            filename = self.path + "/" + content
            self.detect_modetype(filename)

        self.update_modeselector()
        
    def detect_modetype(self, filename):

        s = 20

        for mode in self.unique_modes:
            if re.findall(mode.nametag, filename):
                
                self.terminal_tracker("Detector", "matched [%s] with [%s]" %  \
                                          (filename.split("/")[-1].center(s), \
                                               self.unique_modes_names[self.unique_modes.index(mode)].center(s)))
                
                args = [0]
                mode = mode(*args)
                mode.filepath = filename
                
                if self.check_consistency_fail(mode):
                    self.raiseWarning("Similar dataset previously selected: " + str(mode))
                    break
                
                self.modeMap[str(mode)] = mode
                break

    def check_consistency_fail(self, mode):
        return str(mode) in self.modeMap.keys()

    def load_ext(self):

        config = open("config.txt", 'r')
        self.config = config.read()
        config.close()

        self.autodetect_modes()

        self.refresh_dt_raw = re.findall(r"dynamic refresh interval \[seconds\]\s*=\s*(\d+\.?\d*)", self.config)[0]
        self.refresh_dt_config = float(self.refresh_dt_raw)
        self.terminal_silence = not bool(int(re.findall("terminal tracker\s*=\s*([01])", self.config)[0]))
        self.warning_silence = not bool(int(re.findall("warnings on\s*=\s*([01])", self.config)[0]))

        #Static overriding
        for mode in self.unique_modes:
            mode.refresh_dt = self.refresh_dt_config

    def update_dt(self, event):

        new_dt = self.refresh_dt.get()
        if float(new_dt) <= 0.5:
            self.raiseWarning("Illeagal refresh interval. Choose one > 0.5")
            return

        config = open("config.txt", 'w')
        for line in self.config.split("\n"):
            if line.startswith("dynamic refresh interval [seconds]"):
                config.write(self.config.replace(line, line.replace(self.refresh_dt_raw, new_dt)))
                config.close()
                self.load_ext()
                return

        

    def autodetect_modes(self):
        classfile = open('classfile.py', 'r')
        raw = classfile.read()
        classfile.close()

        self.unique_modes_names = re.findall('^class (\w+)\(dummy\):', raw, re.MULTILINE)
        self.unique_modes = [eval(subclass) for subclass in self.unique_modes_names]
        
        if not self.unique_modes_names:
            self.raiseWarning("No subclass implementations found.")

        for mode in self.unique_modes:
            instance = mode(0)
            try:
                nametag = instance.nametag
            except:
                self.raiseWarning("Subclass %s has no attribute 'nametag' (output filename identifier)." % \
                                  self.unique_modes_names[self.unique_modes.index(mode)])
        
    
    def update_modeselector(self, silent=False):

        self.modeSelector['menu'].delete(0,END)

        if not self.modeMap:

            if not silent:
                self.raiseWarning("Unable to fetch seleceted data")
            self.mode.set("No data...")
            self.start_button['state'] = DISABLED
            self.modeSelector['state'] = DISABLED

            return


        for key in self.modeMap.keys():
            self.modeSelector['menu'].add_command(label=key, \
                                                  command=lambda temp = self.modeMap[key]: \
                                                  self.modeSelector.setvar(self.modeSelector.cget("textvariable"), \
                                                                           value = temp))
        
        self.start_button['state'] = ACTIVE
        self.modeSelector['state'] = ACTIVE
        self.mode.set("Select data...")

    def clean_modeselector(self):
        self.modeMap = {}
        self.update_modeselector(silent=True)
        self.terminal_tracker("GUI", "Cleaned loaded data.")

    def change_dynamic_state(self):

        if self.dynamic.get():
            self.timer_text['state'] = ACTIVE
            self.timer_field['state'] = NORMAL
        else:
            self.timer_text['state'] = DISABLED
            self.timer_field['state'] = DISABLED

    def raiseWarning(self, s):

        self.terminal_tracker("Warning", s)
        
        if self.warning_silence:
            return
        
        tkMessageBox.showwarning(
                    "QMC GUI",
                    s
        )

    def showInfo(self, h, s):
        tkMessageBox.showinfo(h, s)
        
        
    def terminal_tracker(self, source, message, hold='off'):
        if self.terminal_silence:
            return

        s = 10

        if hold=='on':
            self.hide_source = True
            print "[%s]:  %s" % (source.center(s), message),
        elif hold=='off':
            if not self.hide_source:
                print "[%s]: " % source.center(s),
            print message
            self.hide_source = False
        else:
            print "Error: hold=%s is not a leagal value." % hold
        
    def load_images(self):
        
        play = PhotoImage(file='Images/play.gif')
        stop = PhotoImage(file='Images/stop.gif')
        
        self.img = {}
        self.img['play'] = play
        self.img['stop'] = stop

def main():
    
    root = Tk()
    root.title(string="QMC GUI")
    root.wm_iconbitmap('@Images/QMC.xbm')
    
    master = QMCGUI(root)

    root.mainloop()




if __name__ == "__main__":
    main()
