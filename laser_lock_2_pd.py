### Python GUI program for laser locking system
### Built April 2015 by T.J. Procter
### Basic build of gui made with help from http://sebsauvage.net/python/gui/
### Help for the streaming class obtained from tutorials on plot.ly ###

### Change to 3.0 November 2015

### When building with pyinstaller use: pyinstaller "~project.py" --hidden-import=scipy.special._ufuncs_cxx

import wx                               # Import wxpython for GUI
import wx.lib.plot as wxplot            # Import to make gui plots
import visa                             # Import to speak with scope
import time                             # Import to do tricks with time
import datetime                         # Import to use for timestamps
import numpy as np                      # Import for all the numpy goodness
from scipy.optimize import leastsq      # Import for peak fitting/finding
from copy import deepcopy               # Import to allow copies
import threading                        # Import for threads
from ctypes import *                    # C types needed for connection with wavemeter
import plotly                           # Functions needed for streaming program
import plotly.plotly as py   
import plotly.tools as tls    
from plotly.graph_objs import *

global operations, important, useful, graph_axes, graph_axes_2
operations = {}                         # Boolean dictionary of the status of operations
operations["go"] = False
operations["connected"] = False
operations["grabbing"] = False
operations["wave_connected"] = False
operations["streaming"] = False
operations["logging"] = False
operations["displayed"] = False
operations["agilent_grab"] = False

important = {}                          # Dictionary of important variables for program
important["Lock Difference"] = 10.0            # Arbitrary starting values
important["Laser Location"] = 0.0
important["HeNe Left Location"] = 0.0
important["HeNe Right Location"] = 0.0
important["HeNe Sideband Location"] = 0.0
important["HeNe Width"] = 600.0
important["HeNe Center"] = 600.0
important["V Out"] = 0.0
important["Ramp Offset"] = 0.5
important["Ramp Amplitude"] = 2.5
important["Wavenumber"] = 0.0
important["Agilent_1"] = 0.0
important["Agilent_2"] = 0.0
important["Error State"] = 1               # Start in error state
important["Kp"] = 0.0                   # Declare and set feedback constants
important["Ki"] = 0.0
important["Ramp Kp"] = 0.005

global monitor, mon_options
mon_options = []
for key in important:
    mon_options.append(key)
mon_options = sorted(mon_options)
monitor = mon_options[0]                 # Initial monitoring value

useful = {}                             # Dictionary of useful variables
useful["invert"] = 1
useful["invert_2"] = 1
useful["background"] = -0.01         # Initial locking parameters
useful["background_2"] = -0.01
useful["fitwidth"] = 100
useful["lockpos"] = 600.0
useful["view_length"] = 500          # Initial monitoring array length
useful["Kpset"] = 1.5
useful["Kiset"] = 0.0
useful["LockRate"] = 2

graph_axes = [int(0),int(1200),float(-0.1),float(0.0)]            # Graph Properties [xmin, xmax, ymin,ymax]
graph_axes_2 = [int(0),int(1200),float(-0.1),float(0.0)]            # Graph Properties [xmin, xmax, ymin,ymax]

# Create Starting points for fitting peaks: [fwhm, intensity, peak location]
global ptisa0, phenel0, phener0
ptisa0 = np.array([2.5,0.5,600])
phenel0 = np.array([5,0.5,200])
phener0 = np.array([5,0.5,1000])
pside0 = np.array([5,0.5,600])

# Create lists for monitoring information: [Y axis label, List of values, important value to monitor]
global mon_info, wm_mon_info, agilent_mon_info, monitor_xaxis
mon_info, wm_mon_info, agilent_mon_info = {},{},{}
mon_info["Laser Location"] = ["Laser Location (Channel #)",[],"Laser Location"]
mon_info["HeNe Left Location"] = ["Left HeNe Location (Channel #)",[],"HeNe Left Location"]
mon_info["HeNe Right Location"] = ["Right HeNe Location (Channel #)",[],"HeNe Right Location"]
mon_info["HeNe Sideband Location"] = ["HeNe Sideband Location (Channel #)",[],"HeNe Sideband Location"]
mon_info["Lock Difference"] = ["Difference (Channel #)",[],"Lock Difference"]
mon_info["HeNe Width"] = ["Width (Channel #)",[],"HeNe Width"]
mon_info["HeNe Center"] = ["HeNe Center (Channel #)",[],"HeNe Center"]
mon_info["V Out"] = ["Voltage",[],"V Out"]
mon_info["Ramp Offset"] = ["Voltage",[],"Ramp Offset"]
mon_info["Ramp Amplitude"] = ["Voltage",[],"Ramp Amplitude"]
mon_info["Error State"] = ["Error (1/0)",[],"Error State"]
mon_info["Kp"] = ["Kp (arb.)",[],"Kp"]
mon_info["Ki"] = ["Ki (arb.)",[],"Ki"]
mon_info["Ramp Kp"] = ["Ramp Kp (arb.)",[],"Ramp Kp"]

wm_mon_info["Wavenumber"] = ["Wavenumber (cm-1)",[],"Wavenumber"]

agilent_mon_info["Agilent_1"] = ["Reading",[],"Agilent_1"]
agilent_mon_info["Agilent_2"] = ["Reading",[],"Agilent_2"]

mon_length = 1000                           # Total num of values to save
monitor_xaxis = list(xrange(mon_length))    # xAxis list for plotting

def lorentzian(x,p,bkg):                        # Global function that creates lorentzian in range (x) using parameters (p)
    return bkg + (-abs(p[1])) * (((p[0]**2)/4.0)/(((x-p[2])**2)+(p[0]**2/4.0)))

class grabbing_thread(threading.Thread):    # Creating threading class for grabbing scope data
    def __init__(self):
        threading.Thread.__init__(self)     # Initiate thread
    def run(self):
        print "Starting to grab data"
        global xaxis1, waveform1, canvas1, xaxis2, waveform2, canvas2
        while operations["grabbing"] == True:
            xaxis1, waveform1 = self.grabdata(pd_channel.GetValue(),useful["invert"])                   # Grab waveform data whilst required
            xaxis2, waveform2 = self.grabdata(pd_2_channel.GetValue(),useful["invert_2"])                   # Grab waveform data whilst required
            time.sleep(0.05)                                    # After grab and sleep create traces for scope view
            mainline = wxplot.PolyLine(zip(xaxis1,waveform1),colour=wx.BLUE)
            exline = wxplot.PolyLine([(graph_axes[0],useful["background"]),(graph_axes[1],useful["background"])])
            henelline = wxplot.PolyLine([(phenel0[2],graph_axes[2]),(phenel0[2],graph_axes[3])],colour=wx.RED)
            henellower = wxplot.PolyLine([(phenel0[2]-useful["fitwidth"],graph_axes[2]),(phenel0[2]-useful["fitwidth"],graph_axes[3])],colour=wx.RED,style=wx.DOT)
            henelupper = wxplot.PolyLine([(phenel0[2]+useful["fitwidth"],graph_axes[2]),(phenel0[2]+useful["fitwidth"],graph_axes[3])],colour=wx.RED,style=wx.DOT)
            henerline = wxplot.PolyLine([(phener0[2],graph_axes[2]),(phener0[2],graph_axes[3])],colour=wx.RED)
            henerlower = wxplot.PolyLine([(phener0[2]-useful["fitwidth"],graph_axes[2]),(phener0[2]-useful["fitwidth"],graph_axes[3])],colour=wx.RED,style=wx.DOT)
            henerupper = wxplot.PolyLine([(phener0[2]+useful["fitwidth"],graph_axes[2]),(phener0[2]+useful["fitwidth"],graph_axes[3])],colour=wx.RED,style=wx.DOT)
            sideline = wxplot.PolyLine([(pside0[2],graph_axes[2]),(pside0[2],graph_axes[3])],colour=wx.RED)
            sidelower = wxplot.PolyLine([(pside0[2]-useful["fitwidth"],graph_axes[2]),(pside0[2]-useful["fitwidth"],graph_axes[3])],colour=wx.RED,style=wx.DOT)
            sideupper = wxplot.PolyLine([(pside0[2]+useful["fitwidth"],graph_axes[2]),(pside0[2]+useful["fitwidth"],graph_axes[3])],colour=wx.RED,style=wx.DOT)
            gc = wxplot.PlotGraphics([mainline,exline,henelline,henellower,henelupper,
                                      henerline,henerlower,henerupper,
                                      sideline,sideupper,sidelower],"HeNe PD","Channel","Voltage")
            canvas1.Draw(gc,xAxis=(graph_axes[0],graph_axes[1]),yAxis=(graph_axes[2],graph_axes[3]))       # Update the graph

            ###### Graph 2

            mainline2 = wxplot.PolyLine(zip(xaxis2,waveform2),colour=wx.BLUE)
            exline = wxplot.PolyLine([(graph_axes_2[0],useful["background_2"]),(graph_axes_2[1],useful["background_2"])])
            laserline = wxplot.PolyLine([(ptisa0[2],graph_axes_2[2]),(ptisa0[2],graph_axes_2[3])],colour='FOREST GREEN')
            laserlower = wxplot.PolyLine([(ptisa0[2]-useful["fitwidth"],graph_axes_2[2]),(ptisa0[2]-useful["fitwidth"],graph_axes_2[3])],colour='FOREST GREEN',style=wx.DOT)
            laserupper = wxplot.PolyLine([(ptisa0[2]+useful["fitwidth"],graph_axes_2[2]),(ptisa0[2]+useful["fitwidth"],graph_axes_2[3])],colour='FOREST GREEN',style=wx.DOT)
            setline2 = wxplot.PolyLine([(useful["lockpos"],graph_axes_2[2]),(useful["lockpos"],graph_axes_2[3])],colour=wx.RED,width=2)
            gc2 = wxplot.PlotGraphics([mainline2,exline,laserline,laserlower,laserupper,setline2],"Laser PD","Channel","Voltage")
            canvas2.Draw(gc2,xAxis=(graph_axes_2[0],graph_axes_2[1]),yAxis=(graph_axes_2[2],graph_axes_2[3]))       # Update the graph
        print "Stopped grabbing scope data"
        
    def grabdata(self,channel,invert):                                         # Function to grab information from the scope
        while True:                                             # Sometimes waveform has errors so this loop catches them
            try:
                rigol.write(":WAVeform:SOURce CHANnel"+str(int(channel)))
                waveform = rigol.query(":WAVeform:DATA?")
                waveform = waveform.encode('ascii','ignore')
                waveform = waveform.split(',')
                waveform.pop(0)
                waveform = np.array(waveform)                   # Make data suitable for use in program
                waveform = waveform.astype(np.float)
                x = [i for i in xrange(len(waveform))]          # Create x axis of values for length of waveform
                return x, waveform * invert                        # Return x axis and waveform
            except Exception,e:
                print str(e)
                print "Error importing values, will try again"

class lock_thread(threading.Thread):        # Create Threading Class for locking routine
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        global important, operations, useful, phenel0, ptisa0, phener0, pside0
        while True:
            try:
                self.timescale = float(rigol.query("TIM:MAIN:SCAL?"))                           # Find the time scale
                self.chan1scale = float(rigol.query(":CHAN1:SCAL?"))                            # Find the voltage scale
                important["Ramp Offset"] = round(float(rigol.query("SOURCE1:VOLTAGE:OFFSET?")),4)    # Find current ramp offset
                important["Ramp Amplitude"] = round(float(rigol.query("SOURCE1:VOLTAGE?")),4)    # Find current ramp amplitude
                important["V Out"] = round(float(rigol.query("SOURCE2:VOLTAGE:OFFSET?")),4)      # Find current vout
                break
            except Exception,e:
                print str(e)
                print "Problem getting scope info, will try again"
        
        # Create Starting points for fitting peaks: [fwhm, intensity, peak location]
        ptisa0 = np.array([2500.*self.timescale,3.*self.chan1scale,float(laserentry.GetValue())])
        phenel0 = np.array([5000.*self.timescale,3.*self.chan1scale,float(henelentry.GetValue())])
        phener0 = np.array([5000.*self.timescale,3.*self.chan1scale,float(henerentry.GetValue())])
        pside0 = np.array([5000.*self.timescale,3.*self.chan1scale,float(sidebandentry.GetValue())])
        important["HeNe Center"] = (phenel0[2]+phener0[2])/2                   # Calcualte center of hene peaks
                         
        useful["lockpos"] = float(phenel0[2])                   # Set lock position to left HeNe, this can be changed later
        important["Kp"] = 0.0                                               # Set starting feedback constants
        important["Ki"] = 0.0
        useful["Kpset"] = float(feedback.GetValue())                     # Set the target feedback constants
        useful["Kiset"] = float(integral.GetValue())
        important["Ramp Kp"] = float(rampfeed.GetValue())                    # Set the feedback constant for the ramp offset

        self.offsetwrite = important["Ramp Offset"]                              # Create var for writing ramp offset
        self.amplitudewrite = important["Ramp Amplitude"]                              # Create var for writing ramp amplitude
        self.voutwrite = important["V Out"]                                  # Create var for writing vout
        self.prevoffs = self.offsetwrite                                    # Set previous write values to know when change has occured
        self.prevamp = self.amplitudewrite
        self.prevvout = self.voutwrite

        self.errors = [0] * 200                                             # Create list for storing errors

        while operations["go"] == True:                                     # Main loop for locking routine
            while True:                                                     # Loop to fit the peaks
                try:
                    ptisa0,self.fitt = self.fitpeak(ptisa0,waveform2,xaxis2,useful["fitwidth"],useful["background_2"])
                    phenel0,self.fithl = self.fitpeak(phenel0,waveform1,xaxis1,useful["fitwidth"],useful["background"])
                    phener0,self.fithr = self.fitpeak(phener0,waveform1,xaxis1,useful["fitwidth"],useful["background"])
                    pside0,self.fits = self.fitpeak(pside0,waveform1,xaxis1,useful["fitwidth"],useful["background"])
                    important["HeNe Center"] = (phenel0[2]+phener0[2])/2       # Set new hene center
                    important["HeNe Width"] = phener0[2] - phenel0[2]        # Set new hene width
                    important["Laser Location"] = ptisa0[2]
                    important["HeNe Left Location"] = phenel0[2]
                    important["HeNe Right Location"] = phener0[2]
                    important["HeNe Sideband Location"] = pside0[2]
                    break
                except Exception,e:
                    print str(e)    
                    print "Likely something wrong with the fitting routine"
            
            # Check hene peaks have been correctly fitted
            if self.fithl == " Good " and self.fithr == " Good " and round(phenel0[2],1) != round(phener0[2],1):
                self.middle = (xaxis1[0] + xaxis1[-1])/2                      # Find centre of x axis
                self.centdiff = important["HeNe Center"] - self.middle         # Find how far away hene peaks are from center

                self.optwidth = (xaxis1[0] + xaxis1[-1])*2/3                      # Determine optimum HeNe width (usually 800)
                self.widthdiff = important["HeNe Width"] - self.optwidth         # Find how far away hene peaks are from optimum width

                # Change ramp offset and amplitude depending on distance from center of scope.  Normalised by"HeNe Center"                      
                important["Ramp Offset"] += important["Ramp Kp"] * (self.centdiff/important["HeNe Width"])
                self.offsetwrite = round(important["Ramp Offset"],3)             # Round to 3 decimal places
                if self.offsetwrite < (self.prevoffs) or self.offsetwrite > (self.prevoffs+0.002):    # Check value has changed
                    try:
                        rigol.write("SOURCE1:VOLTAGE:OFFSET "+str(self.offsetwrite))    # Write new ramp offset value
                        rigol.write("TRIG:EDG:LEV "+str(self.offsetwrite))              # Write new trigger value
                        rigol.write(":CHAN1:OFFS "+str(-1*self.offsetwrite))            # Write new view offset
                        rampoff.SetValue(self.offsetwrite)                              # Update GUI
                        self.prevoffs = self.offsetwrite                                # Save old value for next check
                    except:
                        pass

                important["Ramp Amplitude"] += important["Ramp Kp"] * (self.widthdiff/important["HeNe Width"])
                self.amplitudewrite = round(important["Ramp Amplitude"],3)             # Round to 3 decimal places
                if self.amplitudewrite < (self.prevamp) or self.amplitudewrite > (self.prevamp+0.002):    # Check value has changed
                    try:
                        rigol.write("SOURCE1:VOLTAGE "+str(self.amplitudewrite))    # Write new ramp amplitude value
                        rampamp.SetValue(self.amplitudewrite)                              # Update GUI
                        self.prevamp = self.amplitudewrite                                # Save old value for next check
                    except:
                        pass

            if abs(important["Kp"]) < abs(useful["Kpset"]):                                 
                important["Kp"] += np.sign(useful["Kpset"])*0.05                     # If not at set feedback value increment up
            if abs(important["Kp"]) > abs(useful["Kpset"]):
                important["Kp"] = useful["Kpset"]                                    # If above the set feedback value change to the set value

            if abs(important["Ki"]) < abs(useful["Kiset"]):                                 
                important["Ki"] += np.sign(useful["Kiset"])*0.002                    # If not at set feedback value increment up
            if abs(important["Ki"]) > abs(useful["Kiset"]):
                important["Ki"] = useful["Kiset"]                                    # If above the set feedback value change to the set value

            if freerun.GetValue():                                                              # If freerun is checked then feedback constants are set to 0
                important["Kp"] = 0.0
                important["Ki"] = 0.0
                
            if (self.fitt == " Good " and self.fithl == " Good " and                       # Check fits are good and make sure two peaks haven't fit in the same place
                self.fithr == " Good " and                                                
                round(phenel0[2],1) != round(phener0[2],1)):
                if important["Error State"] == 1:                                          # Update gui if previously had an error
                    main_label.SetLabel("Back Running")
                important["Error State"] = 0                                               # Set error to 0
                # Find laser lock difference from setpoint (relative to hene center)
                if setpoint.GetValue() == "Left HeNe":          
                    useful["lockpos"] = important["HeNe Left Location"]                   # Change set point and set Kp to 0 for smoother shift
                elif setpoint.GetValue() == "Right HeNe":
                    useful["lockpos"] = important["HeNe Right Location"] 
                elif setpoint.GetValue() == "Sideband":
                    useful["lockpos"] = important["HeNe Sideband Location"]
                elif setpoint.GetValue() == "HeNe Center":
                    useful["lockpos"] = important["HeNe Center"]

                important["Lock Difference"] = ptisa0[2] - useful["lockpos"]
                self.error = (important["Lock Difference"]/important["HeNe Width"])             # Define error value. Normalise to hene width.
                self.errors.append(self.error)                                          # Add error to errors
                self.errors.pop(0)                                                      # Remove old value
                 # Set vout depending on distance from set point.  (Proportional) + (Integral)
                important["V Out"] -= (important["Kp"] * self.error) + (important["Ki"] * (sum(self.errors)/len(self.errors)))        
                if abs(important["V Out"]) > 5.0:
                    important["V Out"] = np.sign(important["V Out"])*5.0                  # Limit vout to +-5.0
                if abs(important["V Out"]) > 2.0:
                    main_label.SetLabel("WARNING! Approaching Vout Limit!")
                self.voutwrite = round(important["V Out"],3)                             # Round write voltage to 3 dp
                if self.voutwrite != self.prevvout:
                    try:
                        rigol.write("SOURCE2:VOLTAGE:OFFSET "+str(self.voutwrite))      # If different write the new vout
                        dcoff.SetValue(self.voutwrite)                                  # Update the GUI
                        self.prevvout = self.voutwrite                                  # Save old value for next check
                    except:
                        pass
            else:                                               # If any of the above fails...
                important["Error State"] = 1                       # Set error state as 1
                important["Kp"] = 0.0                           # Set feedback constant to 0
                important["Ki"] = 0.0
                time.sleep(0.5)                                 # Wait 0.5 seconds
                main_label.SetLabel("Locking Error...")         # Update GUI and try again 

            laserentry.SetValue(round(ptisa0[2],1))             # Update laser location in GUI
            henelentry.SetValue(round(phenel0[2],1))            # Update hene left location in GUI                    
            henerentry.SetValue(round(phener0[2],1))            # Update hene right location in GUI
            sidebandentry.SetValue(round(pside0[2],1))
            currentfeedback.SetLabel(str(important["Kp"]))      # Update current feedback constants in GUI
            currentintegral.SetLabel(str(important["Ki"]))
            self.fit_update(henelfit,self.fithl)                # Update laser fit state in GUI
            self.fit_update(henerfit,self.fithr)                # Update hene left fit state in GUI
            self.fit_update(laserfit,self.fitt)                 # Update hene right fit state in GUI
            self.fit_update(sidebandfit,self.fits)                 # Update hene right fit state in GUI
            laser_diff.SetLabel(str(round(important["Lock Difference"],3)))         # Update laser difference from lock position in GUI
            if (abs(important["Kp"]) == 0 or abs(important["Kp"]) == 0.0) and (abs(important["Ki"]) == 0 or abs(important["Ki"]) == 0.0):
                second_label.SetLabel(" Free Running ")           # Update GUI as free running if Kpset = 0
                second_label.SetBackgroundColour(wx.YELLOW)
            elif abs(np.average(mon_info["Lock Difference"][1][-10:])) < 1.5:
                second_label.SetLabel("Locked")                 # If lock difference is below 1.5 update GUI as locked
                second_label.SetBackgroundColour(wx.GREEN)
            elif abs(np.average(mon_info["Lock Difference"][1][-10:])) < 2.5:
                second_label.SetLabel("Almost Locked")          # If lock difference is below 2.5 update GUI as almost locked
                second_label.SetBackgroundColour(wx.YELLOW)
            else:
                second_label.SetLabel("Trying to lock...")      # Otherwise update GUI as trying to lock
                second_label.SetBackgroundColour(wx.RED)
            
            for key in mon_info:                                # Update monitoring info
                mon_info[key][1].append(important[mon_info[key][2]])
                if len(mon_info[key][1]) > mon_length:
                    mon_info[key][1].pop(0)
            if monitor != "Wavenumber" and monitor != "Agilent_1" and monitor != "Agilent_2":          # If not requesting wavenumber or agilent reading, update the monitor graph
                monitor_data = mon_info[monitor][1]
                monitor_plot = zip(monitor_xaxis[:useful["view_length"]],monitor_data[-useful["view_length"]:])
                monitor_line = wxplot.PolyLine(monitor_plot)
                mon_gc = wxplot.PlotGraphics([monitor_line],"Monitor","Entry",mon_info[monitor][0])
                mon_canvas.Draw(mon_gc,yAxis=(min(monitor_data[-useful["view_length"]:])-0.0001,max(monitor_data[-useful["view_length"]:])+0.0001))
            try:
                time.sleep(useful["LockRate"])                                 # Sleep time before starting loop again (i.e crude refresh rate)
            except:
                pass

    def fit_update(self,label,condition):                       # Function to update GUI goodness of fit
        label.SetLabel(condition)
        if condition == " Good ":
            label.SetBackgroundColour(wx.GREEN)
        else:
            label.SetBackgroundColour(wx.RED)

    def residuals(self,p,y,x,bkg):                                  # Function that takes fit parameters, waveform data (y) and fit range (x)
        return np.subtract(y,lorentzian(x,p,bkg))                   # Returns difference between data (y) and fit
    
    def fitpeak(self,pstart,waveform,xaxis,fitwidth,bkg):           # Function to fit peaks and determine if fit is good or not
        lowbound = int(pstart[2])-fitwidth                      # Set lower boundary to search for peak
        if lowbound < 0:                                        # If this value is below zero set boundary to zero
            lowbound = int(0)
        highbound = int(pstart[2])+fitwidth                     # Set higher boundary (this can go beyond scope range)
        fity = deepcopy(waveform[lowbound:highbound])           # Set data range within boundaries
        fity[fity > bkg] = bkg    # Set values above background (exclusion level) to background
        fitx = deepcopy(xaxis[lowbound:highbound])              # Set the x axis range for creating peak
        try:
            p0 = leastsq(self.residuals,pstart,args=(fity,fitx,bkg))# Least squares minimising routine to fit peak
            pend0 = deepcopy(p0[0])                             # Create array of fitted parameters
        except Exception,e:
            print str (e)
            pend0 = [0.0,0.0,0.0]                               # If error set all values to 0 so a bad fit is declared
        
        if (abs(pend0[0]) > 0.1*5000.*self.timescale and        # Judge fit condition. Fwhm and intensities within expected values
            abs(pend0[0]) < 10.*5000.*self.timescale and
            abs(pend0[1]) > 0.05*3.*self.chan1scale and
            abs(pend0[1]) < 10.*3.*self.chan1scale):            
            fit = " Good "                                        # Declare the fit as good
        else:
            fit = " Bad "                                         # Otherwise declare fit as bad
            pend0[0] = 5000.*self.timescale                     # Set the FWHM to a reasonable value
            pend0[1] = 3.*self.chan1scale                       # Set the intensity to a reasonable value
            pend0[2] = pstart[2]                                # Set the location back to where it started looking
        return pend0, fit                                       # Return fitted parameters and goodness of fit

class agilent_thread(threading.Thread):                         # Creating thread for grabing agilent info
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        global important, operations, agilent
        print "Starting grabbing info from agilent"

        while operations["agilent_grab"] == True:
            time.sleep(1)
            try:
                data = agilent.query("READ?")
                data = data.split(",")
                important["Agilent_1"] = float(data[0].split(" ")[0])
                print important["Agilent_1"]
                important["Agilent_2"] = float(data[4].split(" ")[0])
                print important["Agilent_2"]
            except Exception, e:
                print str(e)
            for key in agilent_mon_info:                                             # Update monitoring info
                    agilent_mon_info[key][1].append(important[agilent_mon_info[key][2]])  
            if len(agilent_mon_info[key][1]) > mon_length:
                agilent_mon_info[key][1].pop(0)
            if monitor == "Agilent_1" or monitor == "Agilent_2":                                         # If requesting wavenumber, update the monitor graph
                monitor_data = agilent_mon_info[monitor][1]
                monitor_plot = zip(monitor_xaxis[:useful["view_length"]],monitor_data[-useful["view_length"]:])
                monitor_line = wxplot.PolyLine(monitor_plot)
                mon_gc = wxplot.PlotGraphics([monitor_line],"Agilent","Entry",agilent_mon_info[monitor][0])
                mon_canvas.Draw(mon_gc,yAxis=(min(monitor_data[-useful["view_length"]:]),max(monitor_data[-useful["view_length"]:])))
        print "Closing Agilent grabbing"

class wavemeter_thread(threading.Thread):                       # Creating threading class for commumation with wavemeter
    def __init__(self):
        threading.Thread.__init__(self)                         # Initiate thread
    def run(self):
        global important, operations
        print "Starting Connection with Wavemeter"
        self.handles = CDLL('bristol_wavemeter\\CLDevIFace.dll')        # Import Handles to cummunicate with the wavemeter
        self.handles.CLGetLambdaReading.restype = c_double              # Set 'GetLambdaReading' to a c type double
        self.handles.CLGetPowerReading.restype = c_double               # Set 'GetLambdaReading' to a c type double

        self.opendevice = self.handles.CLOpenUSBSerialDevice            # Load open device funstion
        self.closedevice = self.handles.CLCloseDevice                   # Load close device function
        self.getlambda = self.handles.CLGetLambdaReading                # Load get lamda function
        self.getpower = self.handles.CLGetPowerReading                  # Load get power function

        self.CommPort = int(comm_port.GetValue())                                        # Define the CommPort the wavemeter is on
        print "Using Value:", str(self.CommPort)                        # Print to confirm
        self.DevHandle = self.opendevice(self.CommPort)                 # Open communication with wavemeter on port number
        if self.DevHandle == 1:
            print "Connected to Device"                                 # Connect to device and print values
            comm_port.Disable()
            wavemeter_button.Disable()
            important["Wavenumber"] = self.getlambda(self.DevHandle)
            self.power = self.getpower(self.DevHandle)
            print "Wavenumber: ", important["Wavenumber"]
            print "Power: ", self.power
        else:
            print "Failed to connect to device"                         # Print warning if connection fails
            operations["wave_connected"] = False
            comm_port.Enable()
            wavemeter_button.Enable()

        while operations["wave_connected"] == True:                     # Initiate loop for grabbing data
            time.sleep(0.1)
            try:
                important["Wavenumber"] = self.getlambda(self.DevHandle)# Grab wavenumber
                self.power = self.getpower(self.DevHandle)              # Grab power
            except Exception, e:
                print str(e)
            if operations["wave_connected"] == True:                    # Annoyingly this has to be included in the loop or else it hangs up on writing to the GUI
                wn_label.SetLabel(str(round(important["Wavenumber"],3))+" cm-1")    # Set the wavenumber in the GUI
                power_label.SetLabel(str(round(self.power,3))+" mW")                # Set the power in the GUI

                for key in wm_mon_info:                                             # Update monitoring info
                    wm_mon_info[key][1].append(important[wm_mon_info[key][2]])  
                if len(wm_mon_info[key][1]) > mon_length:
                    wm_mon_info[key][1].pop(0)
                if monitor == "Wavenumber":                                         # If requesting wavenumber, update the monitor graph
                    monitor_data = wm_mon_info[monitor][1]
                    monitor_plot = zip(monitor_xaxis[:useful["view_length"]],monitor_data[-useful["view_length"]:])
                    monitor_line = wxplot.PolyLine(monitor_plot)
                    mon_gc = wxplot.PlotGraphics([monitor_line],"Monitor","Entry",wm_mon_info[monitor][0])
                    mon_canvas.Draw(mon_gc,yAxis=(min(monitor_data[-useful["view_length"]:]),max(monitor_data[-useful["view_length"]:])))
        print "Closing Wavemeter Device"
        comm_port.Enable()
        wavemeter_button.Enable()
        self.DevHandle = self.closedevice(self.CommPort)                # Close connection with device

class streaming_thread(threading.Thread):                       # Creating threading class for streaming program
    def __init__(self):
        threading.Thread.__init__(self)                         # Initiate thread
    def run(self): 
        print plotly.__version__                                # Check version
        stream_ids = tls.get_credentials_file()['stream_ids']   # Get stream token ids
        print "Window for stream will open.\n If needed, log in using laser as username and the usual password."

        trace1 = Scatter(x=[],y=[],xaxis='x1',yaxis='y1',mode='lines',
            line=Line(width=0.75),name='Photodiode',stream = Stream(token=stream_ids[1]))       # Initialize trace of oscilloscope signal

        trace2 = Scatter(x=[],y=[],mode='markers',name='Lock Diff',xaxis='x3',yaxis='y3',
            marker=Marker(color=[]),stream = Stream(token=stream_ids[2],maxpoints = 1000))      # Initialize trace of lock error(difference)

        trace3 = Scatter(x=[600,600],y=[-0.02,-0.4],xaxis='x1',yaxis='y1',
            text=["Laser Locking Stream!","written by T.J. PROCTER 2015"],
            textfont=Font(family='sans serif',size=15,color='red'),
            textposition='bottom',mode='markers+text',name='Error Messages',
            marker=Marker(opacity=0.001),stream = Stream(token=stream_ids[3],maxpoints = 2))    # Initialize trace of error messages to overlay on scope signal

        trace4 = Scatter(x=[0,0],y=[0,1],xaxis='x2',yaxis='y2',text=["Start","Wavenumber"],
            textfont=Font(family='sans serif',size=20,color='green'),
            textposition='bottom',mode='markers+text',name='Error Messages',marker=Marker(opacity=0.001),
            stream = Stream(token=stream_ids[4],maxpoints = 2))                                 # Initialize trace of lock status and wavenumber messages

        # Trace for histogramming lock errors.  Histogram not working with empty array streaming as of 17 March 2015
        trace5 = Histogram(y=[],name='Lock Err Histogram',marker=Marker(color='rgb(148,103,189)'),
            xaxis='x4',yaxis='y4',stream = Stream(token=stream_ids[5],maxpoints = 1000))

        data = Data([trace1,trace2,trace3,trace4,trace5])       # Add traces together

        # Add title to layout object and set axes
        layout = Layout(title='Laser Lock',showlegend=False,
                        xaxis=XAxis(showgrid=False,showticklabels=False,zeroline=False,domain=[0,0.7],anchor='y1'),
                        xaxis2=XAxis(showgrid=False,zeroline=False,domain=[0.7,1],showticklabels=False,anchor='y2'),
                        xaxis3=XAxis(showgrid=False,zeroline=False,domain=[0,0.75],anchor='y3'),
                        xaxis4=XAxis(showgrid=False,zeroline=False,domain=[0.8,1],anchor='y4'),
                        yaxis=YAxis(title='Photodiode Voltage!!',showticklabels=False,showgrid=False,zeroline=False,domain=[0,0.5],anchor='x1'),
                        yaxis2=YAxis(showgrid=False,zeroline=False,showticklabels=False,domain=[0,0.5],range=[-1.5,2.0],anchor='x2'),
                        yaxis3=YAxis(title='Lock Diff kHz (wrt FSR HeNe)',titlefont=Font(color='rgb(148,103,189)'),showticklabels=False,showgrid=False,
                                     tickfont=Font(color='rgb(148,103,189)'),domain=[0.6,1],anchor='x3'),
                        yaxis4=YAxis(showgrid=False,showticklabels=True,domain=[0.6,1],anchor='x4')                                
                        )

        fig = Figure(data=data, layout=layout)                  # Make a figure object
        unique_url = py.plot(fig, filename='s7_laser_lock')     # Send fig to Plotly, initialize streaming plot, open new tab

        s1 = py.Stream(stream_ids[1])                           # Make instances of the Stream link objects
        s2 = py.Stream(stream_ids[2])
        s3 = py.Stream(stream_ids[3])
        s4 = py.Stream(stream_ids[4])
        s5 = py.Stream(stream_ids[5])

        s1.open()                                   # Open the streams
        s2.open()
        s3.open()
        s4.open()
        s5.open()

        time.sleep(5)                               # Delay start of stream by 5 sec (time to switch tabs)

        s1_data = {}                                # Empty dictionary for stream 1 data
        sendytisa = 0.0
        lockoldold = 10.0                           # Memory of lock values
        lockold = 10.0                              # More recent memory of lock values
        print "Starting STreaming"
        while operations["streaming"] == True:      # Streaming loop
            try:
                s1_data['x'] = np.array(xaxis1)      # Set s1 x data to scope xaxis1
                s1_data['y'] = np.array(waveform1)   # Set s1 y data to scope trace
                xtime = datetime.datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')
                s1.write(s1_data)                   # Write scope trace
                if abs(important["Lock Difference"]) < 1.0:
                    colorset = 'lightgreen'         # Set colour scheme depending on lock difference
                else:
                    colorset = 'red'
                if important["Error State"] == 0:      # If no error update lock difference in kHz using FSR of etalon = 300 MHz
                    sendytisa = important["Lock Difference"]*300000/important["HeNe Width"]
                s2.write(dict(x=xtime,y=sendytisa,marker=Marker(color=colorset)))
                s5.write(dict(y=sendytisa))
                if important["Error State"] == 1:      # If error show messages...
                    s3.write(dict(text=["LOCKING ERROR"," "]))
                    s4.write(dict(text=[" ",str(round(important["Wavenumber"],4))+" cm-1"]))
                elif important["Lock Difference"] == lockold and important["Lock Difference"] == lockoldold:
                    s3.write(dict(text=["LOCKING ERROR","Locking Program Might Not Be Running"]))
                    s4.write(dict(text=[" ",str(round(important["Wavenumber"],4))+" cm-1"]))
                elif (abs(important["Kp"]) == 0 or abs(important["Kp"]) == 0.0) and (abs(important["Ki"]) == 0 or abs(important["Ki"]) == 0.0):
                    s3.write(dict(text=[" "," "]))
                    s4.write(dict(text=["Free Running",str(round(important["Wavenumber"],4))+" cm-1"]))
                elif abs(np.average(mon_info["Lock Difference"][1][-10:])) < 1.5:
                    s3.write(dict(text=[" "," "]))
                    s4.write(dict(text=["Locked",str(round(important["Wavenumber"],4))+" cm-1"]))
                else:
                    s3.write(dict(text=[" "," "]))
                    s4.write(dict(text=[" ",str(round(important["Wavenumber"],4))+" cm-1"]))
                lockoldold = lockold                # Pass on old value
                lockold = important["Lock Difference"]     # Save value
                stream_label.SetLabel("Data Sent: "+str(datetime.datetime.now().strftime('%H:%M:%S')))      # Update GUI of time sent
            except Exception,e:
                print str(e)
            time.sleep(4)
        s3.write(dict(text=["Program Stopped","Program Stopped"]))
        s4.write(dict(text=["Stopped","Stopped"]))
        s1.close()                                  # Close the streams
        s2.close()
        s3.close()
        s4.close()
        s5.close()
        print "Streams closed"

class logging_thread(threading.Thread):                        # Creating threading class for logging data
    def __init__(self):
        threading.Thread.__init__(self)                        # Initiate thread
    def run(self):
        print "Logging Data"
        # Create files for logging data
        self.filename1 = "logs\\" + str(time.asctime()).replace(':','') + ".csv"
        self.filename2 = "M:\Laser_Spec\Laser_Logs\\" + str(time.asctime()).replace(':','') + ".csv"
        self.f1 = open(self.filename1,"a")
        self.f2 = open(self.filename2,"a")
        self.txt = str()
        for key in mon_options:
            self.txt += key
            self.txt += ','
        self.txt += "Time Stamp\n"
        self.f1.write(self.txt)
        self.f2.write(self.txt)
        log_status.SetLabel(self.filename1)                     # Show local filename in GUI
        while operations["logging"] == True:
            self.txt = str()
            for key in mon_options:
                self.txt += str(important[key])
                self.txt += ','
            self.txt += str(datetime.datetime.now().strftime('%c')) + "\n"
            if operations["go"] == True:
                self.f1.write(self.txt)
                self.f2.write(self.txt)
                log_label.SetLabel("Saving Log: "+str(datetime.datetime.now().strftime('%H:%M:%S')))
            time.sleep(5)
        print "Stopping data logging"
        self.f1.close()
        self.f2.close()

class PlotGraph(wxplot.PlotCanvas):                     # Class to create canvas for showing scope
    def __init__(self,parent,id):
        wxplot.PlotCanvas.__init__(self,parent,id,style=wx.BORDER_NONE)
        self.data = ([0,0],[0,1])
        line = wxplot.PolyLine(self.data)
        gc = wxplot.PlotGraphics([line],"Nothing To Display","Channel","Voltage")
        self.Draw(gc,xAxis=(0,1200),yAxis=(-1,0))
        #self.Draw(gc,xAxis=(0,1200),yAxis=(0,6))       # Added for other etalon

class PlotMonitor(wxplot.PlotCanvas):                   # Class to create canvas for showing monitor info
    def __init__(self,parent,id):
        wxplot.PlotCanvas.__init__(self,parent,id,style=wx.BORDER_NONE)
        monitor_line = wxplot.PolyLine([(0,0),(1,1)])
        mon_gc = wxplot.PlotGraphics([monitor_line],"Monitor","Entry","Value")
        self.Draw(mon_gc)

class laser_lock_wx(wx.Frame):                                  # Create Class for GUI
    def __init__(self,parent,id,title):                         # Initialize with info on parent, id and name
        wx.Frame.__init__(self,parent,id,title)                 # Create a frame for the GUI
        self.parent = parent                                    # Keep track of parent if needed
        self.initialize(id)                                     # Start initialize function

    ####### Declare globals

    global important, useful, operations

    ####### Main operation functions

    def start_lock(self,event):                                 # Function to start locking thread
        if operations["connected"] == False:                    # Check that the scope has been connected
            main_label.SetLabel("Scope not connected!")
        elif operations["go"] == True:
            main_label.SetLabel("Program Already Running")      # Check program isn't already running
        else:
            operations["go"] = True                             # If not then set go to true
            main_label.SetLabel("Starting Locking Program")     # Update GUI
            self.reset_button.Disable()
            self.start_button.Disable()
            self.stop_button.Enable()
            rampoff.Disable()
            rampamp.Disable()
            rampfeed.Disable()
            dcoff.Disable()
            feedback.Disable()
            integral.Disable()
            excludeentry.Disable()
            excludeentry2.Disable()
            thread1 = lock_thread()                             # Create thread
            thread1.start()                                     # Start thread

    def stop_lock(self,event):                                  # Function to stop locking
        if operations["connected"] == False:                    # Check that the scope has been connected
            main_label.SetLabel("Scope not connected!")
        elif operations["go"] == False:
            main_label.SetLabel("No Program Running")           # Check program is running
        else:
            operations["go"] = False                            # Set go to False
            important["Kp"] = 0.0                               # Set feedbacks to 0.0
            important["Ki"] = 0.0
            main_label.SetLabel("Stopped!")
            second_label.SetLabel("Stopped!")
            self.reset_button.Enable()
            self.start_button.Enable()
            self.stop_button.Disable()
            rampoff.Enable()
            rampamp.Enable()
            rampfeed.Enable()
            dcoff.Enable()
            feedback.Enable()
            integral.Enable()
            excludeentry.Enable()
            excludeentry2.Enable()
            
    def start_wavemeter(self,event):                            # Function to start wavemeter thread
        if operations["wave_connected"] == True:
            main_label.SetLabel("Wavemeter already connected")
        else:
            operations["wave_connected"] = True
            main_label.SetLabel("Starting Wavemeter Program" )
            self.wave_label.SetLabel("Connected")
            comm_port.Disable()
            wavemeter_button.Disable()
            thread2 = wavemeter_thread()                
            thread2.start()  

    def start_stream(self,event):                               # Function to start streaming thread
        if operations["streaming"] == True:
            main_label.SetLabel("Already streaming data" )
        else:
            operations["streaming"] = True
            main_label.SetLabel("Starting Streaming Program" )
            self.stream_status.SetLabel("Streaming")
            self.stream_button.Disable()
            thread3 = streaming_thread()
            thread3.start()

    def start_logging(self,event):                              # Function to start logging data thread
        if operations["logging"] == True:
            main_label.SetLabel("Already logging data")
        else:
            operations["logging"] = True
            main_label.SetLabel("Logging Data")
            self.logging_button.Disable()
            thread4 = logging_thread()                    
            thread4.start()       

    def quit_program(self,event):                               # Function to quit whole program safely
        print "Closing all connections, will quit in 5 seconds!"
        for key in operations:  operations[key] = False         # Set all operations to false
        time.sleep(5)                                           # Give time for programs to stop
        try:
            rigol.close()                                       # Close connection to rigol if connected
            print "Closed connection with rigol"
        except:
            print "Rigol not currently connected"
            pass                                       
        self.Close(True)                                        # Destroy window    

    ####### Functions dealing with the scope

    def display_devices(self,event):                            # Function to display connected devices in terminal
        global rm,devices
        main_label.SetLabel("See Terminal for Device List")
        rm = visa.ResourceManager()                                                 # Create a resource manager to open up communications with the Rigol Device
        devices = rm.list_resources()                                               # Find resources to connect to
        print "These are the devices connected to the computer:\n"
        for i in xrange(len(devices)):
            print i, devices[i]
        operations["displayed"] = True
        if operations["connected"] == False:  
            self.connect_button.Enable()
            self.scope_num.Enable()
            pd_channel.Enable()
            pd_2_channel.Enable()
        if operations["agilent_grab"] == False:
            self.agilent_start_button.Enable()
            self.agilent_num.Enable()
            channel1_num.Enable()
            channel2_num.Enable()
            

    def connect_scope(self,event):                              # Function to connect scope
        global rigol,rm,devices,thread5
        if operations["connected"] == True:                     # Check scope is not already connected
            main_label.SetLabel("Scope is already connected")
        elif operations["displayed"] == False:
            main_label.SetLabel("Please Display Devices First")
        else:
            try:
                main_label.SetLabel("Connecting...")                                            
                rigol = rm.open_resource(devices[int(self.scope_num.GetValue())])           # Connect to device number stated in GUI
                print "\nConnected with device", rigol.query("*IDN?")
                rigol.write(":WAVeform:SOURce CHANnel"+str(int(pd_channel.GetValue())))            # Make sure scope is set up to outfut waveform from Channel 3
                rampoff.SetValue(round(float(rigol.query("SOURCE1:VOLTAGE:OFFSET?")),4))
                rampamp.SetValue(round(float(rigol.query("SOURCE1:VOLTAGE?")),4))
                dcoff.SetValue(round(float(rigol.query("SOURCE2:VOLTAGE:OFFSET?")),4))      # Set ramp and vout values in GUI depending on how the scope is set
                self.scope_label.SetLabel(" Connected with device ")                          # Change GUI to show scope is connected
                self.scope_label.SetBackgroundColour(wx.GREEN)
                main_label.SetLabel("Connected")
                operations["connected"] = True
                operations["grabbing"] = True
                self.connect_button.Disable()
                self.scope_num.Disable()
                self.reset_button.Enable()
                self.start_button.Enable()
                self.stop_button.Disable()
                pd_channel.Disable()
                pd_2_channel.Disable()
                rampoff.Enable()
                rampamp.Enable()
                rampfeed.Enable()
                feedback.Enable()
                integral.Enable()
                excludeentry.Enable()
                excludeentry2.Enable()
                freerun.Enable()
                dcoff.Enable()
                setpoint.Enable()
                thread5 = grabbing_thread()                                                 # Create thread to grab and plot scope data
                thread5.start()                                                             # Start thread
            except:
                self.scope_label.SetLabel("Failed to connect with device")                  # Connection may fail
                try:
                    rigol.close()                                                           # Close connection whatever device connected
                    print "Closed connection with Device"
                except:
                    print "No Device connected"
                    pass
                pass

    def reset_scope(self,event):                                # Function to reset scope
        if operations["displayed"] == False:
            main_label.SetLabel("Please Display Devices First")
        elif operations["connected"] == False:                    # Check that the scope has been connected
            main_label.SetLabel("Scope not connected") 
        elif operations["go"] == True:
            main_label.SetLabel("Not when locking!")            # Check program isn't already running
        else:
            operations["grabbing"] = False
            time.sleep(3)
            main_label.SetLabel("Resetting scope")                               
            rigol.write("*RST")                                 # Reset the scope so that conditions for locking can be set
            time.sleep(5)                                       # Sleep for 5 seconds to allow scope to restart before set up
            rigol.write(":WAVeform:FORMat ASC")                 # Set up scope
            rigol.write(":CHAN1:DISP ON")
            rigol.write(":CHAN1:PROB 1")
            rigol.write(":CHAN1:SCAL 0.1")
            rigol.write(":CHAN1:OFFS -0.4")
            rigol.write(":CHAN2:DISP ON")
            rigol.write(":CHAN2:PROB 1")
            rigol.write(":CHAN2:SCAL 0.5")
            rigol.write(":CHAN3:DISP ON")
            rigol.write(":CHAN3:PROB 1")
            rigol.write(":CHAN3:SCAL 0.01")
            rigol.write(":CHAN3:OFFS 0.03")
            rigol.write(":CHAN3:INV ON")
            rigol.write(":CHAN4:DISP ON")
            rigol.write(":CHAN4:PROB 1")
            rigol.write(":CHAN4:SCAL 0.01")
            rigol.write(":CHAN4:OFFS 0.02")
            rigol.write(":CHAN4:INV ON")
            rigol.write("TIM:MAIN:SCAL 0.002")
            rigol.write("TRIG:EDG:LEV 0.5")
            rigol.write("TRIG:EDG:SLOP POS")
            
            rigol.write(":SOUR1:FUNC RAMP")                     # Set up the Ramp
            rigol.write(":SOUR1:VOLT:OFFS 0.09")            
            rigol.write(":SOUR1:VOLT 2.49")
            rigol.write(":SOUR1:FREQ 16")
            rigol.write(":SOUR1:OUTP 1")
            
            time.sleep(1)                                       # Trick to get 50% symmetry ramp, scope doesnt do it if set straight at 50
            rigol.write(":SOUR1:FUNC:RAMP:SYMM 49")
            time.sleep(1)
            rigol.write(":SOUR1:FUNC:RAMP:SYMM 50")
            
            rigol.write(":SOUR2:FUNC DC")                       # Set up the TiSa voltage feedback
            rigol.write(":SOUR2:VOLT:OFFS 0.5")
            dcoff.SetValue(0.5)
            rigol.write(":SOUR2:OUTP 1")
            rigol.write(":WAVeform:SOURce CHANnel"+str(int(pd_channel.GetValue())))

            time.sleep(1)
            rigol.write(":ACQ:TYPE AVER")
            rigol.write(":ACQ:AVER 16")

            time.sleep(0.5)
            rigol.write(":SOUR1:VOLT:OFFS 0.40")
            rigol.write(":SOUR1:VOLT 2.30")
            rampoff.SetValue(0.4)
            rampamp.SetValue(2.3)

            print "\nThe scope should have reset and set conditions for the locking system"
            print "\nThe scope should show and be triggered off the ramp, show the photodiode signal and the output voltage to the laser."
            print "*************\nAdjust the bias on the HV amplifier to centralise the two HeNe peaks"
            print "If no photodiode signal the HV amplifier may need to be reset."
            print "Adjust the DC offset of Source 2 to expose the TiSa peak on the screen"
            print "*************\n"
            print "\nThe ramp should be generated from Source 1, viewed on Channel 1 and sent to the HV amplifier."
            print "The TiSa feedback signal should be generated from Source 2, viewed on Channel 2 and sent to the TiSa."
            print "The photodiode signal should be viewed on Channel 3 via the Coherent Analyzer."
            print "Remember:  The HV amplifier can take up to 20 minutes to settle."
            time.sleep(1)
            operations["grabbing"] = True
            time.sleep(1)
            thread10 = grabbing_thread()                         # Create thread
            thread10.start()                                     # Start thread
            main_label.SetLabel("Ready to go!") 

    def move_ramp(self,event):
        if operations["connected"] == False:                                 # Check that the scope has been connected
            main_label.SetLabel("Scope not connected!")
        else:
            if operations["go"] == False:
                rigol.write("TRIG:EDG:LEV "+str(rampoff.GetValue()))         # Set trigger level on scope
                rigol.write(":SOUR1:VOLT:OFFS "+str(rampoff.GetValue()))     # Set ramp offset on scope
                rigol.write(":SOUR1:VOLT "+str(rampamp.GetValue()))     # Set ramp amplitude on scope
                rigol.write(":CHAN1:OFFS "+str(-1*rampoff.GetValue()))
            important["Ramp Offset"] = float(rampoff.GetValue())                  # Change offset for locking program
            important["Ramp Amplitude"] = float(rampamp.GetValue())                  # Change amplitude for locking program

    def move_tisa(self,event):
        if operations["connected"] == False:                                 # Check that the scope has been connected
            main_label.SetLabel("Scope not connected!")
        else:
            if operations["go"] == False:
                rigol.write(":SOUR2:VOLT:OFFS "+str(dcoff.GetValue()))       # Set locking voltage on scope
            important["V Out"] = float(dcoff.GetValue())                      # Change vout for locking program

    def print_help(self,event):
        print "If list of network devices does not load there is a communications error.  try these steps in order before contacting T Procter for help:"
        print "Even if you are not using these units, if they are part of the system they can still hold up Ni-VISA"
        print "This can be checked by opening up NI MAX" 
        print "\n1. Restart Erthernet/GPIB adapter by DAQ computers"
        print "2. Restart Agilent unit"
        print "3. On Rigol scope: Go to Utility menue and change USB device from Computer to PictBridge and back again"
        print "\n4. Restart this Windows Machine"

    ####### Functions dealing with the agilent unit

    def start_agilent_grab(self,event):                              # Function to connect scope
        global agilent,rm,devices
        if operations["agilent_grab"] == True:                     # Check agilent is not already connected
            main_label.SetLabel("Agilent is already logging")
        elif operations["displayed"] == False:
            main_label.SetLabel("Please Display Devices First")
        else:
            try:
                main_label.SetLabel("Connecting...")                                            
                agilent = rm.open_resource(devices[int(self.agilent_num.GetValue())])           # Connect to device number stated in GUI
                print "\nConnected with device", agilent.query("*IDN?")
                agilent.write("ROUT:SCAN (@"+str(channel1_num.GetValue())+","+str(channel2_num.GetValue())+")")
                agilent.write("TRIG:SOUR IMM")
                agilent.write("TRIG:COUNT 1")               
                print agilent.query("READ?")                                                 # Acquire first data point
                operations["agilent_grab"] = True
                self.agilent_start_button.Disable()
                self.agilent_stop_button.Enable()
                self.agilent_num.Disable()
                channel1_num.Disable()
                channel2_num.Disable()
                main_label.SetLabel("Grabbing Agilent Data")
                thread6 = agilent_thread()                                                 # Create thread to grab and plot scope data
                thread6.start()  
            except:
                print "Failed to connect with device"                  # Connection may fail
                pass

    def stop_agilent_grab(self,event):
        if operations["agilent_grab"] == True:
            operations["agilent_grab"] = False
            main_label.SetLabel("Stopping Agilent Data")
            self.agilent_start_button.Enable()
            self.agilent_stop_button.Disable()
            self.agilent_num.Enable()
            channel1_num.Enable()
            channel2_num.Enable()

    ####### Functions for changing variables with widgets

    def change_feedback(self,event):
        useful["Kpset"] = float(event.GetValue())                        # Change the set feedback constant

    def change_rampfeed(self,event):
        important["Ramp Kp"] = float(rampfeed.GetValue())                    # Change the ramp feedback constant

    def change_integral(self,event):         
        useful["Kiset"] = float(integral.GetValue())                     # Change Ki

    def change_background(self,event):
        useful["background"] = float(excludeentry.GetValue())            # Change exclusion level

    def change_background_2(self,event):
        useful["background_2"] = float(excludeentry2.GetValue())            # Change exclusion level

    def reset_constants(self,event):
        important["Kp"] = 0.0                               # Set feedbacks to 0.0
        important["Ki"] = 0.0

    def change_graph(self,event):
        global graph_axes
        if (int(gr_x_max.GetValue()) - 5) > int(gr_x_min.GetValue()) and float(gr_y_max.GetValue()) > float(gr_y_min.GetValue()):
            graph_axes[0] = int(gr_x_min.GetValue())                            # Change graph axes
            graph_axes[1] = int(gr_x_max.GetValue())
            graph_axes[2] = float(gr_y_min.GetValue())
            graph_axes[3] = float(gr_y_max.GetValue())

    def change_graph_2(self,event):
        global graph_axes_2
        if (int(gr_x_max_2.GetValue()) - 5) > int(gr_x_min_2.GetValue()) and float(gr_y_max_2.GetValue()) > float(gr_y_min_2.GetValue()):
            graph_axes_2[0] = int(gr_x_min_2.GetValue())                            # Change graph axes
            graph_axes_2[1] = int(gr_x_max_2.GetValue())
            graph_axes_2[2] = float(gr_y_min_2.GetValue())
            graph_axes_2[3] = float(gr_y_max_2.GetValue())

    def update_params(self,event):
        global phenel0, ptisa0, phener0                         # Update fitting parameters with values from GUI
        phenel0[2] = float(henelentry.GetValue())
        ptisa0[2] = float(laserentry.GetValue())
        phener0[2] = float(henerentry.GetValue())
        pside0[2] = float(sidebandentry.GetValue())
        important["HeNe Center"]=(phener0[2]+phenel0[2])/2

    def OnSelect(self, e):                                                  # Change what monitor graph is showing
        global monitor
        monitor = e.GetString()

    def View_Select(self,e):
        useful["view_length"] = int(e.GetString())

    def InvertData(self,e):
        if self.invert_box.GetValue():
            useful["invert"] = -1
        else:
            useful["invert"] = 1

    def InvertData_2(self,e):
        if self.invert_box_2.GetValue():
            useful["invert_2"] = -1
        else:
            useful["invert_2"] = 1

    ####### Functions to create widgets for GUI

    def create_label(self,id,sizer,text,bkgcol,forecol,x,y,h,w,fontsize):   # Function to create text label
        font = wx.Font(fontsize,wx.FONTFAMILY_DECORATIVE,                   # Define font
                       wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        label = wx.StaticText(panel,id,label=text)                          # Create label and set the text
        label.SetFont(font)                                                 # Set font
        label.SetBackgroundColour(bkgcol)                                   # Change background colour
        label.SetForegroundColour(forecol)                                  # Change foreground colour
        sizer.Add(label,(x,y),(h,w),wx.ALIGN_CENTER,3)                      # Add label to grid
        return label                                                        # Return label for use

    def create_slider(self,id,sizer,text,x,y,h,w,function,low,high):        # Function to create text entry
        slide = wx.Slider(panel,id,value=text,style=wx.SL_HORIZONTAL,minValue=low,maxValue=high)                 # Create entry and set text
        sizer.Add(slide,(x,y),(h,w),wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND,3)                                 # Add entry to grid
        if function:
            self.Bind(wx.EVT_SCROLL,function,slide)                         # Create function action if one is sent
        return slide                                                        # Return entry for use

    def create_spinctrl(self,id,sizer,text,x,y,h,w,function,inc,low,high):  # Function to create value control
        spin = wx.SpinCtrlDouble(panel,id,value=text)                       # Create control and set text
        spin.SetRange(low,high)
        spin.SetIncrement(inc)
        sizer.Add(spin,(x,y),(h,w),wx.ALL|wx.ALIGN_CENTER_VERTICAL,3)     # Add control to grid
        if function:
            spin.Bind(wx.EVT_SPINCTRLDOUBLE,function)                       # Create function action if one is sent
        return spin                                                         # Return control for use

    def create_button(self,id,sizer,text,x,y,h,w,fontsize,function):        # Function to create button
        font = wx.Font(fontsize,wx.FONTFAMILY_DECORATIVE,                   # Define font
                       wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        button = wx.Button(panel,id,label=text)                             # Create button and set text
        button.SetFont(font)                                                # Set font
        sizer.Add(button,(x,y),(h,w),wx.ALL|wx.ALIGN_CENTER,3)              # Add button to grid
        self.Bind(wx.EVT_BUTTON,function,button)                            # Assign function to button
        return button                                                       # Return button for use

    def create_combo(self,id,sizer,val,x,y,h,w,function,options):  # Function to create value control
        combo = wx.ComboBox(panel,id,value=val,choices=options,style=wx.CB_READONLY,size=(120,30))                       # Create control and set text
        sizer.Add(combo,(x,y),(h,w),wx.ALL|wx.ALIGN_CENTER,3)     # Add control to grid
        if function:
            combo.Bind(wx.EVT_COMBOBOX,function)                       # Create function action if one is sent
        return combo 

    def create_invert(self,id,sizer,val,x,y,h,w):  # Function to create value control
        invert_box = wx.CheckBox(panel,id,val)
        invert_box.SetValue(False)
        sizer.Add(invert_box,(x,y),(h,w),wx.ALL|wx.ALIGN_CENTER,3)     # Add control to grid
        return invert_box

    ####### Create GUI Function

    def initialize(self,id):
        global panel
        panel = wx.Panel(self,id)

        # Create Sizers
        self.mainSizer = wx.FlexGridSizer(1,2,1,1)
        self.mainSizer.AddGrowableCol(1,1)
        self.graphSizer = wx.GridSizer(3,1,1,1)
        self.graph_propertiesSizer = wx.FlexGridSizer(2,1,1,1)
        self.graph_propertiesSizer.AddGrowableRow(1,1)
        self.graph_propertiesSizer.AddGrowableCol(0,1)
        self.graph_properties_2_Sizer = wx.FlexGridSizer(2,1,1,1)
        self.graph_properties_2_Sizer.AddGrowableRow(1,1)
        self.graph_properties_2_Sizer.AddGrowableCol(0,1)
        self.monitorSizer = wx.FlexGridSizer(2,1,1,1)
        self.monitorSizer.AddGrowableRow(1,1)
        self.monitorSizer.AddGrowableCol(0,1)
        self.topSizer = wx.BoxSizer(wx.VERTICAL)
        self.titleSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.statusSizer = wx.GridSizer(1,2,1,1)
        self.doubleSizer = wx.GridSizer(1,2,1,1)
        self.controlSizer = wx.GridSizer(1,2,1,1)
        self.waveSizer = wx.GridSizer(2,2,1,1)
        self.logSizer = wx.GridSizer(1,3,1,1)
        self.streamSizer = wx.GridSizer(1,3,1,1)
        self.endSizer = wx.GridSizer(1,2,1,1)
        
        self.scopebox = wx.StaticBox(panel, label='VISA Devices')
        self.scopeBoxSizer = wx.StaticBoxSizer(self.scopebox,wx.VERTICAL)
        self.scopeSizer = wx.GridBagSizer(hgap=1, vgap=1)

        self.settingsbox = wx.StaticBox(panel, label='Locking Settings')
        self.settingsBoxSizer = wx.StaticBoxSizer(self.settingsbox,wx.VERTICAL)
        self.settingsSizer = wx.GridBagSizer(hgap=1, vgap=1)

        self.locksbox = wx.StaticBox(panel, label='Locking Information')
        self.lockBoxSizer = wx.StaticBoxSizer(self.locksbox,wx.VERTICAL)
        self.laserSizer = wx.GridBagSizer(hgap=1, vgap=1)
        self.lockSizer = wx.GridBagSizer(hgap=1, vgap=1)

        # Create Title for title sizer
        font = wx.Font(16,wx.DECORATIVE,wx.NORMAL,wx.NORMAL)
        title = wx.StaticText(panel,id,'Polarizer Locking Program')
        title.SetFont(font)
        self.titleSizer.Add(title,0,wx.ALL,5)

        # Create labels fo status sizer
        global main_label, second_label
        main_label = wx.StaticText(panel,id,'Not Doing Anything')
        main_label.SetBackgroundColour(wx.NullColour)                       # Change background colour
        main_label.SetForegroundColour(wx.BLACK)                      # Change foreground colour
        main_label.SetFont(font)

        second_label = wx.StaticText(panel,id,'No Status')
        second_label.SetBackgroundColour(wx.NullColour)                       # Change background colour
        second_label.SetForegroundColour(wx.BLACK)                      # Change foreground colour
        second_label.SetFont(font)

        self.statusSizer.AddMany([(main_label,0,wx.ALIGN_CENTER),(second_label,0,wx.ALIGN_CENTER)])

        # Fill Scope Sizer
        global channel1_num,channel2_num
        self.display_button = self.create_button(id,self.scopeSizer,'  Display Devices  ',0,0,1,2,12,self.display_devices)
        self.scope_label = self.create_label(id,self.scopeSizer,'   Scope Not Connected   ',wx.RED,wx.BLACK,1,0,1,2,12)
        self.device_label = self.create_label(id,self.scopeSizer,'Scope Device Num:',wx.NullColour,wx.BLACK,2,0,1,1,10)
        self.scope_num = self.create_spinctrl(id,self.scopeSizer,'0',2,1,1,1,None,1,0,10)
        self.scope_num.Disable()
        self.connect_button = self.create_button(id,self.scopeSizer,'  Connect Scope  ',3,0,1,1,12,self.connect_scope)
        self.connect_button.Disable()
        self.reset_button = self.create_button(id,self.scopeSizer,'  Reset Scope  ',3,1,1,1,10,self.reset_scope)
        self.reset_button.Disable()
        self.agilent_label = self.create_label(id,self.scopeSizer,'Agilent Device Num:',wx.NullColour,wx.BLACK,4,0,1,1,10)
        self.agilent_num = self.create_spinctrl(id,self.scopeSizer,'3',4,1,1,1,None,1,0,10)
        self.channel1_label = self.create_label(id,self.scopeSizer,'Agilent_1 Channel Num:',wx.NullColour,wx.BLACK,5,0,1,1,10)
        channel1_num = self.create_spinctrl(id,self.scopeSizer,'205',5,1,1,1,None,1,0,1000)
        self.channel2_label = self.create_label(id,self.scopeSizer,'Agilent_2 Channel Num:',wx.NullColour,wx.BLACK,6,0,1,1,10)
        channel2_num = self.create_spinctrl(id,self.scopeSizer,'206',6,1,1,1,None,1,0,1000)
        self.agilent_num.Disable()
        channel1_num.Disable()
        channel2_num.Disable()
        self.agilent_start_button = self.create_button(id,self.scopeSizer,'Start Agilent Grab',7,0,1,1,10,self.start_agilent_grab)
        self.agilent_stop_button = self.create_button(id,self.scopeSizer,'Stop Agilent Grab',7,1,1,1,10,self.stop_agilent_grab)
        self.agilent_start_button.Disable()
        self.agilent_stop_button.Disable()
        self.scopeBoxSizer.Add(self.scopeSizer,0,wx.CENTER)

        # Fill Settings Sizer
        global rampoff, rampamp, rampfeed, feedback, currentfeedback, integral, currentintegral, pd_channel, pd_2_channel
        self.pd_label = self.create_label(id,self.settingsSizer,'HeNe PD Channel:',wx.NullColour,wx.BLACK,0,0,1,1,10)
        pd_channel = self.create_spinctrl(id,self.settingsSizer,str(3),0,1,1,1,None,1,1,4)
        pd_channel.Disable()
        self.pd_2_label = self.create_label(id,self.settingsSizer,'Laser PD Channel:',wx.NullColour,wx.BLACK,1,0,1,1,10)
        pd_2_channel = self.create_spinctrl(id,self.settingsSizer,str(4),1,1,1,1,None,1,1,4)
        pd_2_channel.Disable()
        self.ramplabel = self.create_label(id,self.settingsSizer,'Ramp Offset:',wx.NullColour,wx.BLACK,2,0,1,1,10)
        rampoff = self.create_spinctrl(id,self.settingsSizer,str(important["Ramp Offset"]),2,1,1,1,self.move_ramp,0.01,-2,2)
        rampoff.Disable()
        self.rampamplabel = self.create_label(id,self.settingsSizer,'Ramp Amplitude:',wx.NullColour,wx.BLACK,3,0,1,1,10)
        rampamp = self.create_spinctrl(id,self.settingsSizer,str(important["Ramp Amplitude"]),3,1,1,1,self.move_ramp,0.01,-3,3)
        rampamp.Disable()
        self.rampfeed_label = self.create_label(id,self.settingsSizer,'Ramp Feedback:',wx.NullColour,wx.BLACK,4,0,1,1,10)
        rampfeed = self.create_spinctrl(id,self.settingsSizer,str(important["Ramp Kp"]),4,1,1,1,self.change_rampfeed,0.001,0,0.5)
        rampfeed.Disable()
        self.feedbacklabel = self.create_label(id,self.settingsSizer,'Kp Set:',wx.NullColour,wx.BLACK,5,0,1,1,10)
        feedback = self.create_spinctrl(id,self.settingsSizer,str(useful["Kpset"]),5,1,1,1,self.change_feedback,0.1,0,2)
        self.feedbackcurrentlabel = self.create_label(id,self.settingsSizer,'Current Kp:',wx.NullColour,wx.BLACK,6,0,1,1,10)
        currentfeedback = self.create_label(id,self.settingsSizer,str(important["Kp"]),wx.NullColour,wx.BLACK,6,1,1,1,10)
        self.integral_label = self.create_label(id,self.settingsSizer,'Ki Set:',wx.NullColour,wx.BLACK,7,0,1,1,10)
        integral = self.create_spinctrl(id,self.settingsSizer,str(useful["Kiset"]),7,1,1,1,self.change_integral,0.1,0,1)
        feedback.Disable()
        integral.Disable()
        self.currentintegral_label = self.create_label(id,self.settingsSizer,'Current Ki:',wx.NullColour,wx.BLACK,8,0,1,1,10)
        currentintegral = self.create_label(id,self.settingsSizer,str(important["Ki"]),wx.NullColour,wx.BLACK,8,1,1,1,10)      
        self.settingsBoxSizer.Add(self.settingsSizer,0,wx.CENTER)

        # Fill Laser & Lock Sizer
        global henelentry, sidebandentry, laserentry, henerentry, dcoff, setpoint, laser_diff, henelfit, laserfit, henerfit, sidebandfit, excludeentry, excludeentry2, freerun
        self.henellabel = self.create_label(id,self.laserSizer,'Left HeNe Loc:',wx.NullColour,wx.BLACK,0,0,1,2,10)
        self.sidebandlabel = self.create_label(id,self.laserSizer,'Sideband Loc:',wx.NullColour,wx.BLACK,0,2,1,2,10)
        self.henerlabel = self.create_label(id,self.laserSizer,'Right HeNe Loc:',wx.NullColour,wx.BLACK,0,4,1,2,10)
        self.laserlabel = self.create_label(id,self.laserSizer,'Laser Loc:',wx.NullColour,wx.BLACK,0,8,1,2,10)

        henelentry = self.create_slider(id,self.laserSizer,phenel0[2],1,0,1,2,self.update_params,0,1200)
        sidebandentry = self.create_slider(id,self.laserSizer,pside0[2],1,2,1,2,self.update_params,0,1200)
        henerentry = self.create_slider(id,self.laserSizer,phener0[2],1,4,1,2,self.update_params,0,1200)
        laserentry = self.create_slider(id,self.laserSizer,ptisa0[2],1,8,1,2,self.update_params,0,1200)

        henelfit = self.create_label(id,self.laserSizer,' Not Yet ',wx.YELLOW,wx.BLACK,2,0,1,2,12)
        sidebandfit = self.create_label(id,self.laserSizer,' Not Yet ',wx.YELLOW,wx.BLACK,2,2,1,2,12)
        henerfit = self.create_label(id,self.laserSizer,' Not Yet ',wx.YELLOW,wx.BLACK,2,4,1,2,12)
        laserfit = self.create_label(id,self.laserSizer,' Not Yet ',wx.YELLOW,wx.BLACK,2,8,1,2,12)

        self.lockBoxSizer.Add(self.laserSizer,0,wx.CENTER)

        self.excludelabel = self.create_label(id,self.lockSizer,'HeNe Cutoff:',wx.NullColour,wx.BLACK,1,0,1,1,10)
        excludeentry = self.create_spinctrl(id,self.lockSizer,str(useful["background"]),1,1,1,1,self.change_background,0.01,-50,0)
        excludeentry.Disable()
        self.excludelabel2 = self.create_label(id,self.lockSizer,'Laser Cutoff:',wx.NullColour,wx.BLACK,1,2,1,1,10)
        excludeentry2 = self.create_spinctrl(id,self.lockSizer,str(useful["background_2"]),1,3,1,1,self.change_background_2,0.01,-50,0)
        excludeentry2.Disable()
        freerun = self.create_invert(id,self.lockSizer,"Free Run",1,5,1,1)
        freerun.Disable()

        self.dclabel = self.create_label(id,self.lockSizer,'V Out:',wx.NullColour,wx.BLACK,2,0,1,1,10)
        dcoff = self.create_spinctrl(id,self.lockSizer,str(important["V Out"]),2,1,1,1,self.move_tisa,0.01,-3,3)
        self.setpointlabel = self.create_label(id,self.lockSizer,'Set Point:',wx.NullColour,wx.BLACK,2,2,1,1,10)
        setpoint = self.create_combo(id,self.lockSizer,"Left HeNe",2,3,1,1,self.reset_constants,["Left HeNe","Right HeNe","Sideband", "HeNe Center"])
        dcoff.Disable()
        setpoint.Disable()
        self.diff_label = self.create_label(id,self.lockSizer,'Lock Diff:',wx.NullColour,wx.BLACK,2,4,1,1,10)
        laser_diff = self.create_label(id,self.lockSizer,str(important["Lock Difference"]),wx.NullColour,wx.BLACK,2,5,1,1,10)
        self.lockBoxSizer.Add(self.lockSizer,0,wx.CENTER)

        # Fill Control Sizer
        self.start_button = wx.Button(panel,id,label='START LOCK')
        self.Bind(wx.EVT_BUTTON,self.start_lock,self.start_button)
        self.stop_button = wx.Button(panel,id,label='STOP LOCK')
        self.Bind(wx.EVT_BUTTON,self.stop_lock,self.stop_button)
        self.start_button.Disable()
        self.stop_button.Disable()
        self.controlSizer.AddMany([(self.start_button,0,wx.EXPAND),(self.stop_button,0,wx.EXPAND)])

        # Fill Wavemeter Sizer
        global wn_label, power_label, wavemeter_button, comm_port
        wavemeter_button = wx.Button(panel,id,label='Connect Wavemeter')
        wavemeter_button.Enable()
        self.Bind(wx.EVT_BUTTON,self.start_wavemeter,wavemeter_button)
        font = wx.Font(20,wx.DECORATIVE,wx.NORMAL,wx.NORMAL)
        wn_label = wx.StaticText(panel,id,label='Wave Number:')
        wn_label.SetFont(font)
        self.miniwavesizer = wx.GridSizer(1,3,1,1)
        self.comm_label = wx.StaticText(panel,id,label='CommPort Num:')
        self.miniwavesizer.Add(self.comm_label,0,wx.LEFT,5)
        comm_port = wx.SpinCtrlDouble(panel,id,value=str(5),size=(60,20))
        comm_port.SetRange(0,10)
        comm_port.SetIncrement(1)
        comm_port.Enable()
        self.miniwavesizer.Add(comm_port,0,wx.LEFT,5)
        self.wave_label = wx.StaticText(panel,id,label='Not Connected')
        self.miniwavesizer.Add(self.wave_label,0,wx.LEFT,5)
        power_label = wx.StaticText(panel,id,label='Power:')
        power_label.SetFont(font)
        self.waveSizer.AddMany([(wavemeter_button,0,wx.ALIGN_CENTER|wx.EXPAND),(wn_label,0,wx.ALIGN_CENTER),
                                (self.miniwavesizer,0,wx.ALIGN_CENTER),(power_label,0,wx.ALIGN_CENTER)])

        # Fill Log Sizer
        global log_status, log_label
        self.logging_button = wx.Button(panel,id,label="Log Data")
        self.Bind(wx.EVT_BUTTON,self.start_logging,self.logging_button)
        self.logging_button.Enable()
        log_status = wx.StaticText(panel,id,label='Not Logging')
        log_label = wx.StaticText(panel,id,label='Not Logging')
        self.logSizer.AddMany([(self.logging_button,0,wx.ALIGN_CENTER|wx.EXPAND),(log_status,0,wx.LEFT|wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL,10),
                                (log_label,0,wx.ALIGN_CENTER)])

        # Fill Stream Sizer
        global stream_label
        self.stream_button = wx.Button(panel,id,label="Stream to Web")
        self.Bind(wx.EVT_BUTTON,self.start_stream,self.stream_button)
        self.stream_button.Enable()
        self.stream_status = wx.StaticText(panel,id,label='Not Streaming')
        stream_label = wx.StaticText(panel,id,label='Not Streaming')
        self.streamSizer.AddMany([(self.stream_button,0,wx.ALIGN_CENTER|wx.EXPAND),(self.stream_status,0,wx.LEFT|wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL,10),
                                (stream_label,0,wx.ALIGN_CENTER)])

        # Fill End Sizer
        self.help_button = wx.Button(panel,id,label="Help")
        self.Bind(wx.EVT_BUTTON,self.print_help,self.help_button)
        self.quit_button = wx.Button(panel,id,label="Quit")
        self.Bind(wx.EVT_BUTTON,self.quit_program,self.quit_button)
        self.endSizer.AddMany([(self.help_button,0,wx.EXPAND),(self.quit_button,0,wx.EXPAND)])
                
        # Add sub-sizers and lines to topSizer
        self.topSizer.Add(self.titleSizer,0,wx.CENTER)
        self.topSizer.Add(wx.StaticLine(panel),0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(self.statusSizer,0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(wx.StaticLine(panel),0,wx.ALL|wx.EXPAND,5)
        self.doubleSizer.Add(self.scopeBoxSizer,0,wx.EXPAND)
        self.doubleSizer.Add(self.settingsBoxSizer,0,wx.EXPAND)
        self.topSizer.Add(self.doubleSizer,0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(self.lockBoxSizer,0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(self.controlSizer,0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(wx.StaticLine(panel),0,wx.ALL|wx.EXPAND,3)
        self.topSizer.Add(self.waveSizer,0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(wx.StaticLine(panel),0,wx.ALL|wx.EXPAND,3)
        self.topSizer.Add(self.logSizer,0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(self.streamSizer,0,wx.ALL|wx.EXPAND,5)
        self.topSizer.Add(wx.StaticLine(panel),0,wx.ALL|wx.EXPAND,3)
        self.topSizer.Add(self.endSizer,0,wx.ALL|wx.EXPAND,5)

        # Fill graph properties sizer

        global gr_x_min, gr_x_max, gr_y_min, gr_y_max
        self.graphlabelsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.x_min_label = wx.StaticText(panel,id,label='X axes min:')
        self.graphlabelsizer.Add(self.x_min_label,0,wx.LEFT,5)
        gr_x_min = wx.SpinCtrlDouble(panel,id,value=str(graph_axes[0]),size=(60,20))
        gr_x_min.SetRange(graph_axes[0]-1000,graph_axes[0]+1000)
        gr_x_min.SetIncrement(1)
        gr_x_min.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph)
        self.graphlabelsizer.Add(gr_x_min,0,wx.LEFT,5)

        self.x_max_label = wx.StaticText(panel,id,label='X axes max:')
        self.graphlabelsizer.Add(self.x_max_label,0,wx.LEFT,5)
        gr_x_max = wx.SpinCtrlDouble(panel,id,value=str(graph_axes[1]),size=(60,20))
        gr_x_max.SetRange(graph_axes[1]-1000,graph_axes[1]+1000)
        gr_x_max.SetIncrement(1)
        gr_x_max.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph)
        self.graphlabelsizer.Add(gr_x_max,0,wx.LEFT,5)

        self.y_min_label = wx.StaticText(panel,id,label='Y axes min:')
        self.graphlabelsizer.Add(self.y_min_label,0,wx.LEFT,5)        
        gr_y_min = wx.SpinCtrlDouble(panel,id,value=str(graph_axes[2]),size=(60,20))
        gr_y_min.SetRange(graph_axes[2]-1000,graph_axes[2]+1000)
        gr_y_min.SetIncrement(0.1)
        gr_y_min.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph)
        self.graphlabelsizer.Add(gr_y_min,0,wx.LEFT,5)

        self.y_max_label = wx.StaticText(panel,id,label='Y axes max:')
        self.graphlabelsizer.Add(self.y_max_label,0,wx.LEFT,5)        
        gr_y_max = wx.SpinCtrlDouble(panel,id,value=str(graph_axes[3]),size=(60,20))
        gr_y_max.SetRange(graph_axes[3]-1000,graph_axes[3]+1000)
        gr_y_max.SetIncrement(0.1)
        gr_y_max.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph)
        self.graphlabelsizer.Add(gr_y_max,0,wx.LEFT,5)
 
        self.invert_box = wx.CheckBox(panel,id,'Invert')
        self.invert_box.SetValue(False)
        wx.EVT_CHECKBOX(self, self.invert_box.GetId(),self.InvertData)
        self.graphlabelsizer.Add(self.invert_box,0,wx.LEFT,15)

        # Fill graph properties 2 sizer

        global gr_x_min_2, gr_x_max_2, gr_y_min_2, gr_y_max_2
        self.graphlabel_2_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.x_min_label_2 = wx.StaticText(panel,id,label='X axes min:')
        self.graphlabel_2_sizer.Add(self.x_min_label_2,0,wx.LEFT,5)
        gr_x_min_2 = wx.SpinCtrlDouble(panel,id,value=str(graph_axes_2[0]),size=(60,20))
        gr_x_min_2.SetRange(graph_axes_2[0]-1000,graph_axes_2[0]+1000)
        gr_x_min_2.SetIncrement(1)
        gr_x_min_2.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph_2)
        self.graphlabel_2_sizer.Add(gr_x_min_2,0,wx.LEFT,5)

        self.x_max_label_2 = wx.StaticText(panel,id,label='X axes max:')
        self.graphlabel_2_sizer.Add(self.x_max_label_2,0,wx.LEFT,5)
        gr_x_max_2 = wx.SpinCtrlDouble(panel,id,value=str(graph_axes_2[1]),size=(60,20))
        gr_x_max_2.SetRange(graph_axes_2[1]-1000,graph_axes_2[1]+1000)
        gr_x_max_2.SetIncrement(1)
        gr_x_max_2.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph_2)
        self.graphlabel_2_sizer.Add(gr_x_max_2,0,wx.LEFT,5)

        self.y_min_label_2 = wx.StaticText(panel,id,label='Y axes min:')
        self.graphlabel_2_sizer.Add(self.y_min_label_2,0,wx.LEFT,5)        
        gr_y_min_2 = wx.SpinCtrlDouble(panel,id,value=str(graph_axes_2[2]),size=(60,20))
        gr_y_min_2.SetRange(graph_axes_2[2]-1000,graph_axes_2[2]+1000)
        gr_y_min_2.SetIncrement(0.1)
        gr_y_min_2.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph_2)
        self.graphlabel_2_sizer.Add(gr_y_min_2,0,wx.LEFT,5)

        self.y_max_label_2 = wx.StaticText(panel,id,label='Y axes max:')
        self.graphlabel_2_sizer.Add(self.y_max_label_2,0,wx.LEFT,5)        
        gr_y_max_2 = wx.SpinCtrlDouble(panel,id,value=str(graph_axes_2[3]),size=(60,20))
        gr_y_max_2.SetRange(graph_axes_2[3]-1000,graph_axes_2[3]+1000)
        gr_y_max_2.SetIncrement(0.1)
        gr_y_max_2.Bind(wx.EVT_SPINCTRLDOUBLE,self.change_graph_2)
        self.graphlabel_2_sizer.Add(gr_y_max_2,0,wx.LEFT,5)
 
        self.invert_box_2 = wx.CheckBox(panel,id,'Invert')
        self.invert_box_2.SetValue(False)
        wx.EVT_CHECKBOX(self, self.invert_box_2.GetId(),self.InvertData_2)
        self.graphlabel_2_sizer.Add(self.invert_box_2,0,wx.LEFT,15)

        # Fill Monitor sizer
        view_options = ['100','500','1000']
        self.monlabelSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.mon_label = wx.StaticText(panel,id,label='Currently Monitoring:')
        self.monlabelSizer.Add(self.mon_label,0,wx.LEFT,10)
        self.combo = wx.ComboBox(panel,id,value=monitor,choices=mon_options,style=wx.CB_READONLY,size=(200,30))
        self.combo.Bind(wx.EVT_COMBOBOX,self.OnSelect)
        self.monlabelSizer.Add(self.combo,0,wx.LEFT,10)
        self.view_length = wx.ComboBox(panel,id,value=view_options[1],choices=view_options,style=wx.CB_READONLY)
        self.view_length.Bind(wx.EVT_COMBOBOX,self.View_Select)
        self.monlabelSizer.Add(self.view_length,0,wx.ALIGN_RIGHT|wx.LEFT,100)

        global mon_canvas
        mon_canvas = PlotMonitor(panel,id)

        self.monitorSizer.Add(self.monlabelSizer,0,wx.ALL|wx.LEFT,5)
        self.monitorSizer.Add(mon_canvas,0,wx.ALL|wx.EXPAND,5)

        # Fill Graph sizer
        global canvas1, canvas2
        canvas1 = PlotGraph(panel,id)
        canvas2 = PlotGraph(panel,id)

        self.graph_propertiesSizer.Add(self.graphlabelsizer,0,wx.ALL|wx.LEFT,5)
        self.graph_propertiesSizer.Add(canvas1,0,wx.ALL|wx.EXPAND,5)

        self.graph_properties_2_Sizer.Add(self.graphlabel_2_sizer,0,wx.ALL|wx.LEFT,5)
        self.graph_properties_2_Sizer.Add(canvas2,0,wx.ALL|wx.EXPAND,5)

        self.graphSizer.Add(self.graph_propertiesSizer,0,wx.ALL|wx.EXPAND,5)
        self.graphSizer.Add(self.graph_properties_2_Sizer,0,wx.ALL|wx.EXPAND,5)
        self.graphSizer.Add(self.monitorSizer,0,wx.ALL|wx.EXPAND,5)

        # Fill Main sizer
        self.mainSizer.Add(self.topSizer,0,wx.ALL|wx.EXPAND,5)
        self.mainSizer.Add(self.graphSizer,0,wx.ALL|wx.EXPAND,5)

        # Set panel attributes
        panel.SetSizer(self.mainSizer)
        self.SetSizeHints(1300,880,2000,880)                #(Width Min, Height Min, Width Max, Height Max)
        self.topSizer.Fit(self)
        
if __name__  == "__main__":
    app = wx.App()                                              # Creat wx python App
    frame = laser_lock_wx(None,-1,"Laser Lock").Show()          # Create frame using class built for GUI
    print "Welcome to the 2-photodiode laser locking program"
    app.MainLoop()                                              # Start loop for GUI interaction
