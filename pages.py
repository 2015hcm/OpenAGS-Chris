# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 12:47:15 2020

@author: chris
"""
from tkinter import Frame, DoubleVar, StringVar, IntVar, Checkbutton, Spinbox, messagebox, filedialog
from tkinter.ttk import Button, Entry, Label, OptionMenu
import os, math
from util import var_mul, binary_search_find_nearest, binary_search_buried
from itertools import chain, groupby
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from scipy.optimize import curve_fit
#Pages
#TODO: variable names camelCase, function names underscore separated
#TODO: possibly resolve star import

class ReviewFitPage(Frame): 
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.controller = parent
        title = Label(self,text="Fit Review",font=("Arial",14,"bold"))
        title.grid(row=0,column=0,columnspan=6)
        self.f = Figure(figsize=(5,5), dpi=100)
        self.a = self.f.add_subplot(111)             
        Label(self, text="Left Bound:").grid(row=1,column=6)
        Label(self, text="Right Bound:").grid(row=1,column=7)
        peakSelectorButton =  Button(self, text="Manually Choose Element", command=self.send_to_mes)
        peakSelectorButton.grid(row=5,column=6)
        integrateButton = Button(self, text="Manually Integrate and Choose Element", command = self.send_to_mi)
        integrateButton.grid(row=5,column=7)
        acceptButton = Button(self, text="Fit looks good, let the program continue", command=self.submit)
        acceptButton.grid(row=6,column=6)
        ignoreButton = Button(self, text="Ignore this data", command = lambda:self.controller.peakCounter.set(self.controller.peakCounter.get()+1))
        ignoreButton.grid(row=6,column=7)
        self.KnownPeaksFrame = Frame(self)
        Label(self.KnownPeaksFrame, text="Peak Information (Center, Amplitude, Width)").grid(row=0, column=0, columnspan=2)
        
        self.newPeakEntry = Entry(self.KnownPeaksFrame)
        self.newPeakEntry.grid(row=1,column=0)
        self.newPeakEntry.insert(0,"Enter Peak Energy")
        self.newPeakEntry.bind("<FocusIn>",self.clear_new_entry)
        
        Button(self.KnownPeaksFrame, text="Add", command=lambda:self.add_peak(self.newPeakEntry.get())).grid(row=1,column=1)
        
        self.KnownPeaksFrame.grid(row=3,column=6, columnspan=2)
        Button(self, text="Reanalyze", command=self.reanalyze).grid(row=4,column=6,columnspan=2)
        self.energies = []
        self.cps = []
        self.slope = 0
        self.intercept = 0
        self.ctrs = []
        self.annots = []
        self.fill = None
        self.peakGUIList = []
        self.removeBtns = []
    def populate_values(self, startInd, endInd, peaks, variances):
        self.a.cla()
        self.startInd = startInd
        self.energies = self.controller.file1Energies[startInd-20:endInd+21]
        self.cps = self.controller.file1CPS[startInd-20:endInd+21]
        self.variances = list(variances)
        if len(peaks) > 0:
            self.ctrs = np.array(peaks[2::3])
            self.amps = np.array(peaks[3::3])
            self.wids = np.array(peaks[4::3])
            self.slope = peaks[0]
            self.intercept = peaks[1]
            poss = chain.from_iterable(self.get_possibilites(True))
            poss = list(k for k,_ in groupby(poss))
            self.peaks = list(peaks[2:])
            fitX = np.arange(self.energies[20],self.energies[-20],.01)
            self.fitX = fitX
            backgroundY = fitX * peaks[0] + np.array([peaks[1]]*len(fitX))
            self.backgroundY = backgroundY
            fitY = self.controller.multiple_gaussian_and_secant(fitX,*peaks)
            self.a.plot(fitX, fitY, "r-", label = "Fit")
            fullAmps = self.amps + (peaks[0] * self.ctrs + peaks[1]) 
            self.peakPoints = self.a.plot(self.ctrs, fullAmps,"go",label="Found Peaks",picker=6)[0]
            self.peakPoints.set_pickradius(10)
            self.backgroundLine = self.a.plot(fitX,backgroundY,"k-", label = "Background")
            self.f.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)
        self.a.plot(self.energies, self.cps,"bo", label="Observed Data")
        self.a.set_xlim(self.energies[20],self.energies[-20])
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=1,column=0,columnspan=6, rowspan=6)
        self.a.legend(loc='upper right', prop={'size': 8})
        toolbarFrame = Frame(self)
        toolbarFrame.grid(row=7,column=0,columnspan=6, sticky="w")
        toolbar = NavigationToolbar2Tk(self.canvas, toolbarFrame)
        self.spinboxLeftVal = DoubleVar()
        self.spinboxLeft = Spinbox(self, values=self.energies[:-20], command = self.update_right, textvariable = self.spinboxLeftVal)
        self.spinboxLeft.bind("<Return>",lambda x:self.update_right())
        self.spinboxLeft.bind("<FocusOut>",lambda x:self.update_right())
        self.spinboxLeftVal.set(self.energies[20])
        self.spinboxLeft.grid(row=2,column=6)
        self.spinboxRightVal = DoubleVar()
        self.spinboxRight = Spinbox(self, values = self.energies[21:], command = self.update_left, textvariable = self.spinboxRightVal)
        self.spinboxRight.bind("<Return>",lambda x:self.update_left())
        self.spinboxRight.bind("<FocusOut>",lambda x:self.update_left())
        self.spinboxRightVal.set(self.energies[-20])
        self.spinboxRight.grid(row=2,column=7)
        for i in self.peakGUIList:
            i.destroy()
        for i in self.removeBtns:
            i.destroy()
        self.peakGUIList = []
        self.removeBtns = []
        for i in range(len(self.ctrs)):
            self.peakGUIList.append(Label(self.KnownPeaksFrame, text="{:5f}, {:5f}, {:5f}".format(self.ctrs[i],self.amps[i],self.wids[i])))
            self.removeBtns.append(Button(self.KnownPeaksFrame, text="Remove", command=lambda temp=i:self.remove_peak(temp)))
            self.peakGUIList[i].grid(row=i+2,column=0)
            self.removeBtns[i].grid(row=i+2,column=1)
    def remove_peak(self, ind):
        self.peakGUIList[ind].destroy()
        del self.peakGUIList[ind]
        self.removeBtns[-1].destroy()
        del self.removeBtns[-1]
    def add_peak(self, loc):
        try:
            loc = float(loc)
        except ValueError:
            return None
        self.reset_new_entry()
        i=0
        while self.energies[i] < loc:
            i += 1
        self.peakGUIList.append(Label(self.KnownPeaksFrame, text="{:5f}, {:5f}, 1.00".format(self.energies[i], self.cps[i])))
        self.removeBtns.append(Button(self.KnownPeaksFrame, text="Remove", command=lambda temp=len(self.peakGUIList)-1:self.remove_peak(temp)))
        self.peakGUIList[-1].grid(row=len(self.peakGUIList)+1, column=0)
        self.removeBtns[-1].grid(row=len(self.peakGUIList)+1, column=1)
    def on_plot_hover(self, e):
        cont, ind = self.peakPoints.contains(e)
        if cont and self.fill == None:
            x = self.peakPoints.get_xdata()[ind["ind"][0]]
            y = self.peakPoints.get_ydata()[ind["ind"][0]] - (self.slope * x + self.intercept)
            self.fill = self.a.fill_between(self.fitX, self.backgroundY, self.backgroundY + self.controller.multiple_gaussian(self.fitX, x, y, self.wids[np.where(self.ctrs == x)[0][0]]), facecolor='red', alpha=0.5)
            self.canvas.draw()
        else:
            if self.fill != None:
                self.fill.remove()
                self.fill = None
                self.canvas.draw()
    #TODO: add manual entry flag, and call it with a keybind from the spinboxes
    def update_left(self):
        temp = self.spinboxLeftVal.get()
        if self.spinboxRightVal.get() > self.energies[-1]:
            messagebox.showinfo("Out of Bounds","Please enter an X value within the currently plotted points")
            self.spinboxRightVal.set(self.energies[-20])
            return None
        try:
            newRange = self.energies[:np.where(self.energies==self.spinboxRightVal.get())[0][0]]
        except IndexError:
            newRange = self.energies[:binary_search_find_nearest(self.energies, self.spinboxRightVal.get())]
        self.spinboxLeft.configure(values=newRange)
        self.spinboxLeftVal.set(temp)
    def update_right(self):
        temp = self.spinboxRightVal.get()
        if self.spinboxLeftVal.get() < self.energies[0]:
            messagebox.showinfo("Out of Bounds","Please enter an X value within the currently plotted points")
            self.spinboxLeftVal.set(self.energies[20])
            return None
        try:
            newRange = self.energies[np.where(self.energies==self.spinboxLeftVal.get())[0][0]+1:]
        except IndexError:
            newRange = self.energies[binary_search_find_nearest(self.energies, self.spinboxLeftVal.get()):]
        self.spinboxRight.configure(values=newRange)
        self.spinboxRightVal.set(temp)
    def reanalyze(self):
        try:
            lowIndex = np.where(self.energies==self.spinboxLeftVal.get())[0][0]
            highIndex = np.where(self.energies==self.spinboxRightVal.get())[0][0] + 1
        except IndexError:
            lowIndex = binary_search_find_nearest(self.energies, self.spinboxLeftVal.get())
            highIndex = binary_search_find_nearest(self.energies, self.spinboxRightVal.get())
        guesses = [[float(j) for j in i.cget("text").split(", ")] for i in self.peakGUIList]
        guesses = list(chain.from_iterable(guesses))
        guesses = [self.slope, self.intercept] + guesses
        try:
            popt, pcov = curve_fit(self.controller.multiple_gaussian_and_secant, xdata = self.energies[lowIndex:highIndex], ydata = self.cps[lowIndex:highIndex], p0 = np.array(guesses))
            self.populate_values(self.startInd - 20 + lowIndex, self.startInd - 20 + highIndex, popt, np.diag(pcov)[2:])
        except RuntimeError:
            messagebox.showwarning("No fit found", "No good peak fit could be determined. Try entering fewer centroids.")
        
    def clear_new_entry(self, _):
        self.newPeakEntry.delete(0,"end")
    def reset_new_entry(self):
        self.newPeakEntry.delete(0,"end")
        self.newPeakEntry.insert(0,"Enter Peak Energy")
    def get_possibilites(self, MP):
        if MP:#Multiple peaks
            return [self.controller.get_possibilites_list(c) for c in self.ctrs]
        else:
            return self.controller.get_possibilites_list((self.energies[20]+self.energies[-20])/2)
    def send_to_mi(self):
        #TODO: get left and right selectors
        try:
            left = float(self.spinboxLeftVal.get())
            right = float(self.spinboxRightVal.get())
        except ValueError:
            messagebox.showinfo("Bad Input","Please enter in-range floating point numbers for the bounds.")
        if left < self.energies[0] or right > self.energies[-1]:
            messagebox.showinfo("Bad Input","Please enter in-range floating point numbers for the bounds.")
            return None
        r = binary_search_find_nearest(self.energies, right)
        l = binary_search_find_nearest(self.energies, left)
        self.controller.show_frame(ManualIntegrationPage)
        self.controller.frames[ManualIntegrationPage].populate_values(self.energies[l:r], self.cps[l:r])
        self.controller.frames[ManualIntegrationPage].add_peak_selector(self.get_possibilites(False))
        
    def send_to_mes(self):
        cpsNoBgd = [j - (self.energies[i] * self.slope + self.intercept) for i,j in enumerate(self.cps)]
        self.controller.show_frame(ManualElementSelect)
        self.controller.frames[ManualElementSelect].populate_values(self.energies[20:-20], cpsNoBgd[20:-20])
        self.controller.frames[ManualElementSelect].add_peak_selectors(self.peaks,self.get_possibilites(True),self.variances)
    def submit(self):
        self.controller.add_auto_peaks(self.peaks, list(self.variances))
        self.controller.increment_peak_counter()

class ManualIntegrationPage(Frame): 
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.controller = parent
        title = Label(self,text="Manual Integration",font=("Arial",14,"bold"))
        title.grid(row=0,column=0,columnspan=6)
        self.f = Figure(figsize=(5,5), dpi=100)
        self.a = self.f.add_subplot(111)       
        Label(self, text="Match to Element", font=("Arial",14)).grid(row=1,column=5)
        Label(self, text="Lower Bound: ").grid(row=7,column=0)
        self.lowerBound = Entry(self)
        self.lowerBound.grid(row=7,column=1)
        Label(self, text="Upper Bound: ").grid(row=7,column=2)
        self.upperBound = Entry(self)
        self.upperBound.grid(row=7,column=3)
        integrateButton = Button(self, text="Integrate", command = self.submit)
        integrateButton.grid(row=7,column=4)
        Button(self, text="Back",command=self.back).grid(row=8,column=0)
        self.curMenu = None
    def populate_values(self, energies, cps):
        self.energies = energies
        self.cps = cps
        self.a.cla()
        self.a.plot(energies, cps,"b-")
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=1,column=0, columnspan=5, rowspan=5)
        
    def add_peak_selector(self, poss):
        self.poss = poss
        if self.curMenu != None:
            self.curMenu.destroy()
        self.selected = StringVar(self.controller)
        self.selected.set("Select Element")
        opts = ["Select Element"] + [p[0]+": "+str(p[1]) for p in poss] #Stupid tkinter workaround
        self.curMenu = OptionMenu(self,self.selected,*opts)
        self.curMenu.grid(row=2, column = 5)
    def back(self):
        self.controller.show_frame(ReviewFitPage)
    def submit(self):
        try: 
            lowBound = float(self.lowerBound.get())
            highBound = float(self.upperBound.get())
        except ValueError:
            messagebox.showinfo("Error","Use numbers for the bounds please!")
            return None
        if lowBound > highBound:
            messagebox.showinfo("Error","Lower bound cannot be greater than upper bound!")
            return None
        if lowBound < self.energies[0] or highBound > self.energies[-1]:
            messagebox.showinfo("Error","Bounds must be within data bound")
            return None
        if self.selected.get() == "Select Peak":
            messagebox.showinfo("Error","Select a peak match!")
            return None
        #do integration
        #TODO: Infinite loop here i think
        i=0
        while self.energies[i]<lowBound:
            i += 1
        j=i
        while self.energies[j]<highBound:
            j += 1
        #trapezoidal sum
        startEnergy = self.energies[i]
        endEnergy = self.energies[j]
        n = j-i+1
        area = (self.energies[j] - self.energies[i]) / (2 * (j-i+1)) * (sum(self.cps[i:j+1])+sum(self.cps[i+1:j]))
        peak = []
        element, loc = self.selected.get().split(": ")
        for i in self.poss:
            if i[0] == element and str(i[1]) == loc:
                peak = i
        self.controller.add_mi_peak(peak, area, startEnergy, endEnergy, n)
        self.controller.increment_peak_counter()
        self.curMenu.destroy()
class ManualElementSelect(Frame): 
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.controller = parent
        title = Label(self,text="Manual Peak Matching",font=("Arial",14,"bold"))
        title.grid(row=0,column=0,columnspan=6)
        self.f = Figure(figsize=(5,5), dpi=100)
        self.a = self.f.add_subplot(111)        
        Label(self, text="Peak Location: ").grid(row=1,column=4)
        Label(self, text="Match: ").grid(row=1,column=5)
        Button(self, text="Submit", command=self.submit).grid(row=20, column=5)
        self.newObjects = []
        Button(self, text="Back",command=self.back).grid(row=25,column=0)
    def populate_values(self, energies, cps):
        self.a.cla()
        self.a.plot(energies, cps,"b-")
        self.canvas = FigureCanvasTkAgg(self.f, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=1,column=0, columnspan=4, rowspan=20)
        
    def add_peak_selectors(self, peaks, poss, variances):
        for o in self.newObjects:
            o.destroy()
        self.newObjects = []
        self.selected = []
        self.peaks = peaks
        self.variances = variances
        centers = peaks[::3]
        for i in range(len(centers)):
            self.selected.append(StringVar(self.controller))
        self.poss = poss
        for i in range(len(centers)):
            self.selected[i].set("Ignore")
            opts = ["","Ignore"] + [j[0]+": "+str(j[1]) for j in poss[i]] #Stupid tkinter workaround
            #TODO: Dicts of these?
            copy = i
            temp = Label(self,text=str(centers[i]))
            temp.grid(row=i+2, column = 4)
            self.newObjects.append(temp)
            temp=OptionMenu(self,self.selected[copy],*opts)
            temp.grid(row=i+2, column = 5)
            self.newObjects.append(temp)
    def back(self):
        self.controller.show_frame(ReviewFitPage)
    def submit(self): #TODO: Remove all peak locations and option menus
        sel = False
        inp = []
        usedPoss = []
        usedPeaks = []
        usedVars = []
        for i in range(len(self.selected)):
            if self.selected[i].get() != "Ignore":
                inp.append(self.selected[i].get())
                usedPeaks += self.peaks[3*i:3*i+3]
                usedVars += self.variances[3*i:3*i+3]
                sel = True
        if not sel:
            massagebox.showinfo("Error","Please match at least one peak!")
            return None
        for p in inp:
            temp = p.split(": ")
            for q in self.poss:
                for r in q:
                    if r[0] == temp[0] and str(r[1]) == temp[1]:
                        usedPoss.append(r)
        toReturn = []
        for i in range(0, len(usedPeaks), 3):
            ctr = usedPeaks[i]
            amp = usedPeaks[i+1]
            ampVar = usedVars[i+1]
            wid = usedPeaks[i+2]
            widVar = usedVars[i+2]
            mass = amp * abs(wid) * math.sqrt(2*math.pi) / float(usedPoss[i//3][2])
            massDev = math.sqrt(var_mul(amp, ampVar, wid, widVar)) * math.sqrt(2*math.pi) / float(usedPoss[i//3][2])
            toReturn.append([usedPoss[i//3][0], [ctr, amp, wid, amp*wid*math.sqrt(2*math.pi)], float(usedPoss[i//3][1]), float(usedPoss[i//3][2]), mass, massDev])
        self.controller.add_ms_peaks(toReturn)
        self.controller.increment_peak_counter()
        
class ResultsViewer(Frame):
    def __init__(self, parent):
        Frame.__init__(self,parent)
        self.parent=parent
        self.header = Label(self, text="Results by File", font=("Arial",14,"bold"))
        self.header.grid(row=0,column=0,columnspan=3)
        self.currentRow = 1
    def display_file_results(self,fname,results):
        self.parent.show_frame(SingleFileViewer)
        self.parent.frames[SingleFileViewer].display_file_results(fname, results)
    def export_all_data(self,filesDict):
        fname = filedialog.asksaveasfilename(initialdir = ".",title = "Select export file",filetypes = [("CSV files","*.csv")])
        if fname[-4:] != ".csv":
            fname += ".csv"
        exportFile = open(fname,'w')
        for f in filesDict.keys():
            exportFile.write(f+"\n")
            exportFile.write("Peaks Used for Mass Calculations:\n")
            exportFile.write("Element, Theorhetical Peak Energy, Observed Peak Energy, Peak Height, Peak Width, Peak Area, Peak Sensitivity, Predicted Mass (mg), Std. Dev. (mg)\n")
            for element in filesDict[f][0].keys():
                exportFile.write(element+",")
                for peak in filesDict[f][0][element]:
                    outList = [peak[1]] + list(peak[0][:4]) + list(peak[2:])
                    exportFile.write(','.join([str(i) for i in outList])+"\n,")
                exportFile.seek(exportFile.tell() - 1, os.SEEK_SET)  #Don't indent the next element
            exportFile.write("Disregarded Peaks:\n")
            exportFile.write("Observed Peak Energy, Peak Height, Peak Width, Peak Area, Closest Element, Element's Peak Center, Element Sensitivity, Predicted Mass (mg), Std. Dev. (mg)\n")
            for element in filesDict[f][2].keys():
                for peak in filesDict[f][2][element]:
                    outList = list(peak[0][:4]) + [element] + list(peak[1:])
                    exportFile.write(','.join([str(i) for i in outList])+"\n")
                exportFile.seek(exportFile.tell() - 1, os.SEEK_SET)  #Don't indent the next element
            exportFile.write("Masses Calculated:"+"\n")
            exportFile.write("Element, Predicted Mass (mg), Std. Dev. (mg), Number of Atoms (mol), Std. Dev. (mol)\n")
            for element in filesDict[f][1].keys():
                exportFile.write(element + ',' + ','.join([str(i) for i in filesDict[f][1][element]])+"\n")
            exportFile.write("\n\n")
        exportFile.close()
    def show_files(self,filesDict):
        for f in filesDict.keys():
            fileLabel = Label(self, text=f)
            fileLabel.grid(row=self.currentRow,column=0)
            showButton = Button(self, text = "Show", command=lambda tmp=f:self.display_file_results(tmp,filesDict[tmp]))
            showButton.grid(row=self.currentRow,column=1)
            self.currentRow += 1
        exportAllButton = Button(self, text="Export All Data", command=lambda:self.export_all_data(filesDict))
        exportAllButton.grid(row=self.currentRow, column=0,columnspan=2)
class SingleFileViewer(Frame):
    def __init__(self, parent):
        Frame.__init__(self,parent)
        self.parent=parent
        self.currentRow = 0
    
    def write_row(self, itemsList):
        for i in range(len(itemsList)):
            Label(self, text=str(itemsList[i]), borderwidth=4, relief="solid").grid(row=self.currentRow, column=i, sticky="nsew")
        self.currentRow += 1
        
    def display_file_results(self,fname, results):
        self.currentRow = 0
        Button(self, text="Back",command=self.go_back).grid(row=1000,column=0)
        Label(self, text=fname,font=("Arial",14)).grid(row=0,column=0,columnspan=8)
        Label(self, text="Peaks Used for Mass Calculations:").grid(row=1,column=0, columnspan=2, pady=(20,0), sticky="w")
        self.currentRow = 2
        self.write_row(['Element', 'Theorhetical Peak Energy', 'Observed Peak Energy', 'Peak Height', 'Peak Width', 'Peak Area','Peak Sensitivity', 'Predicted Mass (mg)', 'Std. Dev. (mg)'])
        for element in results[0].keys():
            for peak in results[0][element]:
                outList = [element, peak[1]] + list(peak[0][:4]) + list(peak[2:])
                self.write_row(outList)
        Label(self, text="Disregarded Peaks:").grid(row=self.currentRow,column=0, columnspan=2, pady=(20,0), sticky="w") 
        self.currentRow += 1
        self.write_row(['Observed Peak Energy', 'Peak Height', 'Peak Width', 'Peak Area','Closest Element', "Element's Peak Center", 'Element Sensitivity', 'Predicted Mass (mg)', 'Std. Dev. (mg)'])
        for element in results[2].keys():
            for peak in results[2][element]:
                outList = list(peak[0][:4]) + [element] + list(peak[1:])
                self.write_row(outList)
        Label(self, text="Masses Calculated:").grid(row=self.currentRow,column=0, columnspan=2, pady=(20,0), sticky="w") 
        self.currentRow += 1
        self.write_row(['Element', 'Predicted Mass (mg)', 'Std. Dev. (mg)', 'Number of Atoms (mol)', 'Std. Dev. (mol)'])
        for element in results[1].keys():
            self.write_row([element]+results[1][element])
            
    def go_back(self):
        for c in self.winfo_children():
            c.destroy()
        self.parent.show_frame(ResultsViewer)
class BaseFileSelectorFrame(Frame):
    def __init__(self, parent, files):
        Frame.__init__(self, parent)
        self.files = files
        self.controller = parent
        Label(self, text="Select Base Files").grid(row=0,column=0)
        i=1
        self.checkBoxList = []
        for f in files:
            tmpVar = IntVar()
            Checkbutton(text=f,variable=tmpVar).grid(row=i,column=0)
            self.checkBoxList.append(tmpVar)
        Button(text="Next", command=self.send_to_sample_select)
    def send_to_sample_select(self):
        outFiles = []
        for i in range(len(self.files)):
            if self.checkBoxList[i].get():
                outFiles.append(self.files[i])
        self.controller.show_frame(SampleSelectFrame)
        self.controller.frames[SampleSelectFrame].populate_list(outFiles)
class SampleSelectFrame(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.filesList = []
        self.checkBoxList = []
        Label(self, text="Select Sample Files").grid(row=0,column=0)
        Button(self, text="Analyze",command=self.do_analysis).grid(row=1000,column=0)
        
    def populate_list(self, files):
        for currentCheckbox in self.filesList:
            currentCheckbox.grid_forget()
        self.filesList = []
        self.checkBoxList = []
        for newFile in files:
            tmpVar = IntVar()
            tmpBtn = Checkbutton(text=newFile, variable=tmpVar)
            self.fileslist.append(tmpBtn)
            self.checkBoxList.append(tmpVar)
    def do_analysis(self):
        outFiles = []
        for i in range(len(self.filesList)):
            if self.checkBoxList[i].get():
                outFiles.append(self.filesList[i].cget("text"))
        