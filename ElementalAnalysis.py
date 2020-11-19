# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:41:44 2020

@author: chris
"""
import tkinter as tk
from pages import ManualIntegrationPage, ReviewFitPage, ManualElementSelect, ResultsViewer, SingleFileViewer
import numpy as np
import copy, math
from util import peakdet, binary_search_find_nearest, var_mul
from scipy.optimize import curve_fit
from scipy.signal import find_peaks_cwt
class ElementalAnalysisFrame(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self,parent)
        self.parent = parent
        self.peakCounter = tk.IntVar(0)
        self.manualAnalysisResults = {"D":{},"ND":[],"V":[]}
        """
        manualAnalysisResults={"D":{finalPredictions},"ND":[peak1,peak2,peak3], "V":[var1,var2,var3]}
        D: Done, no additional analysis required. add to finalPredictions immediately
        ND: Not Done, send to old matcher function
        V: variances for peaks in ND list
        """
        self.manualChoices = []
        self.frames = {}
        for F in (ResultsViewer, ManualIntegrationPage, ReviewFitPage, ManualElementSelect, SingleFileViewer):
            frame = F(self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(ReviewFitPage)
        self.filesDict = dict()
    def clean_close(self):
        self.peakCounter.set(len(self.fitRegions))
        p = self.parent.master
        self.parent.destroy()
        p.protocol("WM_DELETE_WINDOW", p.destroy)
    def clean_close_all(self):
        self.peakCounter.set(len(self.fitRegions))
        p = self.parent.master
        self.parent.destroy()
        p.destroy()
    def add_all_data(self, filesList, fileInfo, ROIs, targets, reference):
        self.filesList = filesList
        self.fileInfo = fileInfo
        self.ROIs = ROIs
        self.targets = targets
        self.all_peaks_sens = reference
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()
    def increment_peak_counter(self):
        self.peakCounter.set(self.peakCounter.get()+1)
    def add_mi_peak(self,peak, area, lb, ub, n):
        """This function is called when the user chooses the option to manually integrate under a peak, when the fit does not seem good enough. 
        The function performs a trapezoidal area sum and calculates the standard deviation of that sum, then adds the result to the "Done" section of the manualAnalysisResults dictionary."""
        curReg = self.fitRegions[self.peakCounter.get()]
        curReg[0] = self.file1Energies[curReg[0]]
        curReg[1] = self.file1Energies[curReg[1]]
        self.manualChoices.append((curReg, ["MI",lb,ub,peak]))
        i=0
        data = []
        while self.file1Energies[i] < lb - 5:
            i += 1
        while self.file1Energies[i] < lb:
            data.append(self.file1CPS[i])
            i += 1
        while self.file1Energies[i] < ub:
            i += 1
        while self.file1Energies[i] < ub + 5:
            data.append(self.file1CPS[i])
            i += 1
        mean = sum(data) / len(data)
        stdev = (sum([((x - mean) ** 2) for x in data]) / len(data)) ** 0.5
        areaStdev = (ub-lb) / n * (2*n-2) * stdev
        mass = area/float(peak[2])
        massStdev = areaStdev/float(peak[2])
        if peak[0] not in self.manualAnalysisResults["D"].keys():
            self.manualAnalysisResults["D"][peak[0]] = [(["MI","MI","MI", area], float(peak[1]), float(peak[2]), mass, massStdev)]
        else:
            self.manualAnalysisResults["D"][peak[0]].append((["MI","MI","MI", area], float(peak[1]), float(peak[2]), mass, massStdev))
    def add_ms_peaks(self, peaks):
        """This function is for peaks that were manually matched (alternately: manually selected, or ms) to a specific element.
           It will add the peaks to the "Done" section of the manualAnalysisResults dictionary."""
        curReg = self.fitRegions[self.peakCounter.get()]
        curReg[0] = self.file1Energies[curReg[0]]
        curReg[1] = self.file1Energies[curReg[1]]
        self.manualChoices.append((curReg, ["MS", peaks])) #Store the manual matches for use in analysis of other similar files      
        for p in peaks:
            if p[0] in self.manualAnalysisResults["D"].keys():
                self.manualAnalysisResults["D"][p[0]].append(tuple(p[1:]))
            else:
                self.manualAnalysisResults["D"][p[0]] = [tuple(p[1:])]
                
    def add_auto_peaks(self, peaks, variances):
        """"This function is called when the user selects the auto-matching option, and the peaks are simply added to a list where they will be auto matched by the analyzer later."""
        curReg = self.fitRegions[self.peakCounter.get()]
        curReg[0] = self.file1Energies[curReg[0]]
        curReg[1] = self.file1Energies[curReg[1]]
        self.manualChoices.append((curReg, ["A"]))
        self.manualAnalysisResults["ND"] += peaks
        self.manualAnalysisResults["V"] += variances
    
    def run_analysis(self):
        self.file1Energies = self.fileInfo[0][0]
        self.file1CPS = self.fileInfo[0][1]
        fitRegions, variances = self.do_peak_fitting(self.file1Energies, self.file1CPS, self.ROIs)
        self.fitRegions = fitRegions
        while self.peakCounter.get() < len(fitRegions):
            regionData = fitRegions[self.peakCounter.get()]
            enerRegion = self.file1Energies[regionData[0]:regionData[1]+1]
            self.show_frame(ReviewFitPage)
            self.frames[ReviewFitPage].populate_values(regionData[0], regionData[1], regionData[2:], variances[self.peakCounter.get()])
            self.wait_variable(self.peakCounter)
        foundElements, disregardedPeaks = self.find_elements(self.manualAnalysisResults, self.targets)
        finalMasses = self.get_masses(foundElements)
        self.filesDict[self.filesList[0]] = [foundElements, finalMasses, disregardedPeaks]
        for i in range(1, len(self.fileInfo)):
            energies, cps = self.fileInfo[i]
            foundElements, finalMasses, disregardedPeaks = self.reanalyze(energies, cps, self.manualChoices, self.targets)
            self.filesDict[self.filesList[i]] = [foundElements, finalMasses, disregardedPeaks]
        self.show_frame(ResultsViewer)
        self.frames[ResultsViewer].show_files(self.filesDict)

    
    def multiple_gaussian_and_secant(self, x, *params): 
        """This function models a multiple gaussian curve with a linear model for background noise."""
        y = np.zeros_like(x)
        y += x * float(params[0])
        y += float(params[1])
        for i in range(2, len(params), 3):
            ctr = float(params[i])
            amp = float(params[i+1])
            wid = float(params[i+2])
            y = y + amp * np.exp( -((x - ctr)/wid)**2)
        return y
    def multiple_gaussian(self, x, *params): 
        """This function models a multiple gaussian curve."""
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            ctr = float(params[i])
            amp = float(params[i+1])
            wid = float(params[i+2])
            y = y + amp * np.exp( -((x - ctr)/wid)**2)
        return y
    
    def multiple_peak_and_background(self, peaks, background, x, params):
        y = np.zeros_like(x)
        leftCounter = 0
        rightCounter = background.get_num_params()
        
        background.update_params(params[leftCounter:rightCounter])
        y += background.get_ydata(x)
        
        for peak in peaks:
            leftCounter = rightCounter
            rightCounter = leftCounter + peak.get_num_params()
            peak.update_params(params[leftCounter:rightCounter])
            y += peak.get_ydata(x)
            
        return y
    
    def mpb_no_update(self, peaks, background, x):
        y = np.zeros_like(x)
        y += background.get_ydata(x)
        for peak in peaks:
            y += peak.get_ydata(x)
        return y
    def calculate_ROIs(self, centers, energies, cps):
        for center in centers:
            searchLeft = center
            searchRight = center + 1
            v1 = 25
            v2 = .05
            firstPass = True
            while searchRight - searchLeft > 50 or firstPass:
                while not (max(abs(max(cps[searchLeft-4:searchLeft]) - cps[searchLeft]), abs(min(cps[searchLeft-4:searchLeft]) - cps[searchLeft])) > max(cps[searchLeft]/v1, v2)):
                      searchLeft -= 1
                while not (max(abs(max(cps[searchRight:searchRight+4]) - cps[searchRight]), abs(min(cps[searchRight:searchRight+4]) - cps[searchRight])) > max(cps[searchRight]/v1, v2)):
                      searchRight += 1
                v1 /= 2
                v2 *= 2
                firstPass = False
            while searchRight - searchLeft < 15:
                searchLeft -= 1
                searchRight += 1
    def fit_region(self, regionGuess, energies, cps):
        pass
    def do_peak_fitting(self,energies, cps, ROIs):
        """This function takes a lot of guesses at initial parameters, then uses a curve fitter to find the best-fitting peak equation"""
        peaks_to_return = []
        variances = []
        for region in ROIs:
            i = 0
            while energies[i] < region[0]:
                i+=1
            lowerBound = energies[i-1]
            lowerIndex = i-1
            while energies[i] < region[1]:
                i+=1
            upperBound = energies[i]
            upperIndex = i
            peaksInRange = []
            i=0
            while float(self.all_peaks_sens[i][1]) < lowerBound:
                i+=1
            while i < len(self.all_peaks_sens) and float(self.all_peaks_sens[i][1]) < upperBound:
                peaksInRange.append(self.all_peaks_sens[i])
                i+=1
                
            #Background line guessing algorthm (finds flat areas on each side of peak)
            thresh = 5000
            searchLeft = (lowerIndex + upperIndex)//2
            searchRight = (lowerIndex + upperIndex)//2 + 1
            if lowerBound < thresh:
                v1=25
                v2=10
                v3=.05
            else:
                v1=100
                v2=50
                v3=.015
            while not ((abs(sum(cps[searchLeft-4:searchLeft])/4-cps[searchLeft]) < cps[searchLeft]/v1 and abs(cps[searchLeft-4]-cps[searchLeft]) < cps[searchLeft]/v2 and abs(cps[searchLeft-3]-cps[searchLeft]) < cps[searchLeft]/v2 and abs(cps[searchLeft-2]-cps[searchLeft]) < cps[searchLeft]/v2 and abs(cps[searchLeft-1]-cps[searchLeft]) < cps[searchLeft]/v2) or (sum(cps[searchLeft-4:searchLeft+1])<v3)):
                searchLeft -= 1
            if lowerBound < thresh:
                while not ((abs(sum(cps[searchRight+1:searchRight+5])/4-cps[searchRight]) < cps[searchRight]/v1 and abs(cps[searchRight+1]-cps[searchRight]) < cps[searchRight]/v2 and abs(cps[searchRight+2]-cps[searchRight]) < cps[searchRight]/v2 and abs(cps[searchRight+3]-cps[searchRight]) < cps[searchRight]/v2 and abs(cps[searchRight+4]-cps[searchRight]) < cps[searchRight]/v2) or (sum(cps[searchRight+1:searchRight+5])<v3)):
                    searchRight += 1
            else:
                while cps[searchRight]>cps[searchLeft]:
                    searchRight += 1
            while searchRight - searchLeft > 50: #If our ROI is too wide
                v1/=2
                v2/=2
                v3*=2
                searchLeft = (lowerIndex + upperIndex)//2
                searchRight = (lowerIndex + upperIndex)//2 + 1
                while not ((abs(sum(cps[searchLeft-4:searchLeft])/4-cps[searchLeft]) < cps[searchLeft]/v1 and abs(cps[searchLeft-4]-cps[searchLeft]) < cps[searchLeft]/v2 and abs(cps[searchLeft-3]-cps[searchLeft]) < cps[searchLeft]/v2 and abs(cps[searchLeft-2]-cps[searchLeft]) < cps[searchLeft]/v2 and abs(cps[searchLeft-1]-cps[searchLeft]) < cps[searchLeft]/v2) or (sum(cps[searchLeft-4:searchLeft+1])<v3)):
                    searchLeft -= 1
                if lowerBound < thresh:
                    while not ((abs(sum(cps[searchRight+1:searchRight+5])/4-cps[searchRight]) < cps[searchRight]/v1 and abs(cps[searchRight+1]-cps[searchRight]) < cps[searchRight]/v2 and abs(cps[searchRight+2]-cps[searchRight]) < cps[searchRight]/v2 and abs(cps[searchRight+3]-cps[searchRight]) < cps[searchRight]/v2 and abs(cps[searchRight+4]-cps[searchRight]) < cps[searchRight]/v2) or (sum(cps[searchRight+1:searchRight+5])<v3)):
                        searchRight += 1
                else:
                    while cps[searchRight]>cps[searchLeft]:
                        searchRight += 1
            while searchRight - searchLeft < 15: #If our ROI is too narrow, just widen it a bit
                searchLeft -= 1
                searchRight += 1
                
            enerToFit = energies[searchLeft:searchRight+1]
            cpsToFit = cps[searchLeft:searchRight+1]
            start = cpsToFit[0]
            end = cpsToFit[-1]
            tempData = copy.deepcopy(cpsToFit)
            slope = (end - start) / (enerToFit[-1] - enerToFit[0])
            intercept = start - enerToFit[0] * slope
            tempData = tempData - (slope*enerToFit + np.full_like(tempData, intercept))
            newMin = min(tempData)
            tempData = tempData - np.full_like(tempData, newMin)
            
            #Peak finding
            maxima = peakdet(cpsToFit,max(cpsToFit)/50,enerToFit)
            guessedParams = [slope, intercept]
            for i,j in maxima:
                if j < max(cpsToFit)/25:
                    maxima.remove((i,j))
            maxima = sorted(maxima,key=lambda x:-x[1])
            while len(maxima) * 3 + 2 > searchRight - searchLeft:
                del maxima[-1]
            
            #TODO: get width estimation from FWHM input
            for m,_ in maxima:
                guessedParams += [m, cpsToFit[binary_search_find_nearest(enerToFit, m)], 1.00]
            
            #Actual curve fitting
            if len(guessedParams) >= 5:
                try:
                    popt, pcov = curve_fit(self.multiple_gaussian_and_secant, enerToFit, cpsToFit, p0=guessedParams)
                    peaks_to_return.append([searchLeft, searchRight] + list(popt))
                    pvar = list(np.diag(pcov))[2:]
                    if len(pvar) == 0:
                        pvar = [.1,.1,.1]#maybe change
                    variances.append(pvar)
                except RuntimeError:
                    try:
                        #fallack, use our guess for background line and give the algorithm less params to work with
                        popt, pcov = curve_fit(self.multiple_gaussian, enerToFit, tempData, p0=guessedParams[2:])
                        peaks_to_return.append([searchLeft, searchRight] + [slope, intercept] + list(popt))
                        pvar = list(np.diag(pcov))
                        if len(pvar) == 0:
                            pvar = [.1,.1,.1]#maybe change
                        variances.append(pvar)
                    except RuntimeError:
                        peaks_to_return.append([searchLeft, searchRight])
                        variances.append([])
                        #TODO: as a catch-all for special cases, just integrate under the data. Present the graph of the data(with maybe even more context) to the researcher and let them choose the bounds
                        """self.show_frame(ManualIntegrationPage)
                        self.frames[ManualIntegrationPage].populate_values(enerToFit,cpsToFit)"""
                        #TODO: wait until we get a callback
                        #TODO: figure out variance for this case
                        #TODO: figure out how to differentiate integral data from fit data in find_elements
                        pass
            else:
                peaks_to_return.append([searchLeft, searchRight])
                variances.append([])
        """peaks_to_return = [[left bound, right bound, slope, intercept, ctr1, amp1, wid1...,ctrn,ampn,widn]]"""
        return peaks_to_return, variances
    
    def get_possibilites_list(self, ctr):
        """Get potential peaks within a certain range of an energy value"""
        j = 0
        poss = []
        while float(self.all_peaks_sens[j][1]) < ctr - 10:
            j+=1
        while float(self.all_peaks_sens[j][1]) < ctr + 10:
            poss.append(self.all_peaks_sens[j])
            j+=1
        while len(poss)>5:
            if len(poss) == 6:
                poss = poss[:-1]
            else:
                poss = poss[1:-1]
        return poss
    def find_elements(self, manualAnalysisResults, targets):
        """Match peaks to their respective elements"""
        #run through once, find possibilites for each peak (all within a 2-3keV range)
        peaksWithPossibilites = []
        finalPredictions = copy.deepcopy(manualAnalysisResults["D"])
        """
        finalPreditcions = {"Al-28":[(observed peak energy, theorhetical peak energy, peak sensitivity, mass estimate, stddev for mass estimate)]}
        """
        potentialPredictions = dict()
        """
        potentialPredictions = {observed peak energy:{"Al-28":[theorhetical peak energy, peak sensitivity, mass estimeate, stddev for mass estimate]}}
        """
        foundPeaks = manualAnalysisResults["ND"]
        errors = manualAnalysisResults["V"]
        #Finds possible sources for each peak
        for i in range(0,len(foundPeaks),3):
            ctr = foundPeaks[i]
            amp = foundPeaks[i+1]
            ampVar = errors[i+1] #this is a variance
            wid = foundPeaks[i+2]
            widVar = errors[i+2] #so is this
            poss = []
            area = amp * abs(wid) * math.sqrt(2*math.pi)
            areaErr = math.sqrt(var_mul(amp, ampVar, wid, widVar)) * math.sqrt(2*math.pi) #this is now standard deviation
            j = 0
            #TODO: binsearch the first one
            while float(self.all_peaks_sens[j][1]) < ctr + 2:
                j+=1
            j -= 1
            while float(self.all_peaks_sens[j][1]) >= ctr - 2:
                poss.append(self.all_peaks_sens[j])
                j-=1
            #TODO: actually  figure out what to do here
            if len(poss) == 0:
                while float(self.all_peaks_sens[j][1]) < ctr + 7:
                    j+=1
                j-=1
                while float(self.all_peaks_sens[j][1]) >= ctr - 7:
                    poss.append(self.all_peaks_sens[j])
                    j-=1
                
            peaksWithPossibilites.append([[ctr,amp,wid,area,areaErr],poss]) 
        for peak,poss in peaksWithPossibilites:
            if len(poss) == 1:
                sym = poss[0][0]
                if sym in finalPredictions.keys():
                    #check if the predicted peak is already used, if so assign the closer observation to it
                    done = False
                    for pred in finalPredictions[sym]:
                        if float(poss[0][1]) == pred[1]:
                            if abs(peak[0] - pred[0][0]) < 1:
                                pass
                            elif abs(float(poss[0][1]) - peak[0]) < abs(pred[1] - pred[0][0]):
                                finalPredictions[sym].remove(pred)
                            else:
                                done = True
                    if not done:        
                        finalPredictions[sym] += [(peak, float(poss[0][1]), float(poss[0][2]), peak[3]/float(poss[0][2]), peak[4]/float(poss[0][2]))]
                else:
                    finalPredictions[sym] = [(peak, float(poss[0][1]), float(poss[0][2]), peak[3]/float(poss[0][2]), peak[4]/float(poss[0][2]))]
            else:
                peakString = ",".join(map(str,peak))
                potentialPredictions[peakString] = dict()
                for potential_peak in poss:
                   sym = potential_peak[0]
                   potentialPredictions[peakString][sym] = [float(potential_peak[1]),potential_peak[2],peak[3]/float(potential_peak[2]),peak[4]/float(potential_peak[2])]
        for peak in potentialPredictions.keys():
            for sym in potentialPredictions[peak].keys():
                if sym in finalPredictions.keys():
                    bestIndex = 0
                    bestSens = 0
                    for i in range(len(finalPredictions[sym])):
                        if finalPredictions[sym][i][2] > bestSens:
                            bestSens = finalPredictions[sym][i][2]
                            bestIndex = i
                    bestGuess = finalPredictions[sym][bestIndex][3]
                    r = bestGuess / potentialPredictions[peak][sym][2]
                    if r>1:
                        r = 1/r
                    potentialPredictions[peak][sym].append(r)
                else:
                    potentialPredictions[peak][sym].append("N/A")
                        
        #by cross-referencing the predicted mass from our peak and each possible element with the mass that other peaks of each element predict, determine which element has created the peak of interest
        for peak in potentialPredictions.keys():
            ratio = 0
            bestMatch = ""
            for k in potentialPredictions[peak].keys():
                if potentialPredictions[peak][k][-1] != "N/A" and potentialPredictions[peak][k][-1] > ratio:
                    bestMatch = k
                    ratio = potentialPredictions[peak][k][-1]
            if bestMatch in targets: #TODO: Handle case of different isotopes
                finalPredictions[bestMatch].append(([list(map(float, peak.split(",")))] + potentialPredictions[peak][bestMatch][0:4]))
        toDelete = []
        disregardedPeaks = dict()
        for k in finalPredictions.keys():
            if k not in targets:
                toDelete.append(k)
        for k in toDelete:
            del finalPredictions[k]
        for k in finalPredictions.keys(): 
            toDelete = []
            sensList = [float(l[2]) for l in finalPredictions[k]]
            maxSensLoc = sensList.index(max(sensList))
            bestGuess = finalPredictions[k][maxSensLoc][3]
            bestStdev = finalPredictions[k][maxSensLoc][4]
            for peak in finalPredictions[k]:
                r = abs(peak[3] - bestGuess)
                if r > 2*bestStdev:
                    toDelete.append(peak)
            for peak in toDelete:
                finalPredictions[k].remove(peak)
                if k in disregardedPeaks.keys():
                    disregardedPeaks[k].append(peak)
                else:
                    disregardedPeaks[k] = [peak]
        return finalPredictions, disregardedPeaks
    def get_masses(self, finalPredictions):
        """Gets the masses of elements given their peak intensities, using the provided sensitivity file"""
        output = dict()
        for element in finalPredictions.keys():
            masses = [peak[3] for peak in finalPredictions[element]]
            stddevs = [peak[4] for peak in finalPredictions[element]]
            weights = [1/s for s in stddevs]
            mass = sum([a*b for a,b in zip(masses,weights)])/sum(weights)
            finalDev = math.sqrt(sum([(a*b)**2 for a,b in zip(weights, stddevs)]))/sum(weights)
            scale = int(element.split("-")[1]) - 1 #TODO: this should actually be the element's atomic mass, just approx. for now
            output[element] = [mass, finalDev, mass/scale, finalDev/scale]
        return output
    def reanalyze(self,energies, cps, manualChoices, targets):
        """Reruns analysis on files similar to the one that was manually corrected"""
        newAnalysisResults = {"D":{},"ND":[],"V":[]} #same format as manualAnalysisResults
        for region, choice in manualChoices:
            if choice[0] == "MI":
                i = 0
                data = []
                while energies[i] < choice[1] - 5:
                    i += 1
                while energies[i] < choice[1]:
                    data.append(cps[i])
                    i += 1
                j=i
                while energies[j] < choice[2]:
                    j += 1
                k=j
                while energies[k] < choice[2] + 5:
                    data.append(cps[k])
                    k += 1
                n = j-i+1
                area = (energies[j] - energies[i]) / (2 * n) * (sum(cps[i:j+1])+sum(cps[i+1:j]))
                mean = sum(data) / len(data)
                stdev = (sum([((x - mean) ** 2) for x in data]) / len(data)) ** 0.5
                areaStdev = (energies[j] - energies[i]) / n * (2*n-2) * stdev
                mass = area/float(choice[3][2])
                massStdev = areaStdev/float(choice[3][2])
                if choice[3][0] not in newAnalysisResults["D"].keys():
                    newAnalysisResults["D"][choice[3][0]] = [(["MI","MI","MI", area], float(choice[3][1]), float(choice[3][2]), mass, massStdev)]
                else:
                    newAnalysisResults["D"][choice[3][0]].append((["MI","MI","MI", area], float(choice[3][1]), float(choice[3][2]), mass, massStdev))
            else:
                i = 0
                while energies[i] < region[0]:
                    i += 1
                j=i
                while energies[j] < region[1]:
                    j += 1
                try:
                    popt, pcov = curve_fit(self.multiple_gaussian_and_secant, energies[i:j], cps[i:j], p0=region[2:])
                    popt = list(popt)
                    pvar = np.diag(pcov)
                    if len(pvar) == 0 or not np.isfinite(pvar).all():
                        print("oof")
                        region[2] = float(popt[0])
                        region[3] = float(popt[1])
                        raise RuntimeError
                except RuntimeError:
                    print("yay recovery")
                    tempData = copy.deepcopy(cps[i:j])
                    for l,m in enumerate(tempData):
                            tempData[l] = m - (region[2] * energies[i+l] + region[3])
                    popt, pcov = curve_fit(self.multiple_gaussian, energies[i:j], tempData, p0=region[4:])
                    popt = [region[2], region[3]] + list(popt)
                    pvar = [.1,.1] + list(np.diag(pcov))
                    if len(pvar) == 0 or not np.isfinite(pvar).all():
                            print("ouch")
                            pvar = [.1] * len(popt)
                if choice[0] == "MS":
                    for i in range(2, len(popt), 3):
                        ctr = popt[i]
                        amp = popt[i+1]
                        ampVar = pvar[i+1]
                        wid = popt[i+2]
                        widVar = pvar[i+2]
                        for p in choice[1]:
                            if abs(ctr - p[1][0]) < .25 or abs(ctr - p[2]) < .25: #TODO: perhaps change these numbers
                                area = amp * wid * math.sqrt(2*math.pi)
                                stdev = math.sqrt(var_mul(amp, ampVar, wid, widVar)) * math.sqrt(2*math.pi)
                                mass = float(area) / float(p[3])
                                massStd = float(stdev) / float(p[3])
                                if p[0] in newAnalysisResults["D"].keys():
                                    newAnalysisResults["D"][p[0]].append([[ctr, amp, wid, area], p[2], p[3], mass, massStd])
                                else:
                                    newAnalysisResults["D"][p[0]] = [[[ctr, amp, wid, area], p[2], p[3], mass, massStd]]
                                break #I hate to do it but it makes sense
                else:
                    newAnalysisResults["ND"] += list(popt[2:])
                    newAnalysisResults["V"] += list(pvar[2:])
        foundElements, disregardedPeaks = self.find_elements(newAnalysisResults, targets)
        return foundElements, self.get_masses(foundElements), disregardedPeaks

