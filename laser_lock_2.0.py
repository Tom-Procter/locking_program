### Python GUI program for laser locking system
### Built April 2015 by T.J. Procter
### Basic build of gui made with help from http://sebsauvage.net/python/gui/
### Help for the streaming class obtained from tutorials on plot.ly ###

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

global operations, important, xaxis, waveform
operations = {}                         # Boolean dictionary of the status of operations
operations["go"] = False
operations["connected"] = False
operations["grabbing"] = False
operations["wave_connected"] = False
operations["streaming"] = False
operations["logging"] = False

important = {}                          # Dictionary of important variables for program
important["tisadiff"] = 10.0            # Arbitrary starting values
important["henewidth"] = 600.0
important["henecent"] = 600.0
important["vout"] = 0.0
important["offset"] = 0.5
important["wavenumber"] = 0.0
important["view_length"] = 500          # Initial monitoring array length
important["errstate"] = 1               # Start in error state
important["Kp"] = 0.0                   # Declare and set feedback constants
important["Kpset"] = 0.5
important["Ki"] = 0.0
important["Kiset"] = 0.0
important["rampKp"] = 0.01
important["background"] = -0.12         # Initial locking parameters
important["fitwidth"] = 100
important["lockpos"] = 0.0

xaxis = [1,2,3,4,5]                     # Arbitraty starting trace
waveform = [1,1,1,1,1]

# Create Starting points for fitting peaks: [fwhm, intensity, peak location]
global ptisa0, phenel0, phener0
ptisa0 = np.array([2.5,0.5,600])
phenel0 = np.array([5,0.5,200])
phener0 = np.array([5,0.5,1000])

# Create lists for monitoring information: [Y axis label, List of values, important value to monitor]
global mon_info, wm_mon_info, monitor, monitor_xaxis
mon_info, wm_mon_info = {},{}
mon_info["Lock Difference"] = ['Difference (Channel #)',[],"tisadiff"]
mon_info["HeNe Width"] = ['Width (Channel #)',[],"henewidth"]
mon_info["HeNe Center"] = ['HeNe Center (Channel #)',[],"henecent"]
mon_info["V Out"] = ['Voltage',[],"vout"]
mon_info["Ramp Offset"] = ['Voltage',[],"offset"]
mon_info["Error State"] = ['Error (1/0)',[],"errstate"]
mon_info["Kp"] = ['Kp (arb.)',[],"Kp"]
mon_info["Ki"] = ['Ki (arb.)',[],"Ki"]
wm_mon_info["Wavenumber"] = ['Wavenumber (cm-1)',[],"wavenumber"]

monitor = "Lock Difference"                 # Initial monitoring value
mon_length = 1000                           # Total num of values to save
monitor_xaxis = list(xrange(mon_length))    # xAxis list for plotting

def lorentzian(x,p):                        # Global function that creates lorentzian in range (x) using parameters (p)
    return important["background"] + (-abs(p[1])) * (((p[0]**2)/4.0)/(((x-p[2])**2)+(p[0]**2/4.0)))

class grabbing_thread(threading.Thread):    # Creating threading class for grabbing scope data
    def __init__(self):
        threading.Thread.__init__(self)     # Initiate thread
    def run(self):
        print "Starting to grab data"
        global xaxis, waveform, canvas
        while operations["grabbing"] == True:
            xaxis, waveform = self.grabdata()                   # Grab waveform data whilst required
            time.sleep(0.05)                                    # After grab and sleep create traces for scope view
            mainline = wxplot.PolyLine(zip(xaxis,waveform),colour=wx.BLUE)
            exline = wxplot.PolyLine([(0,important["background"]),(1200,important["background"])])
            henelline = wxplot.PolyLine([(phenel0[2],-0.8),(phenel0[2],0.0)],colour=wx.RED)
            henellower = wxplot.PolyLine([(phenel0[2]-important["fitwidth"],-0.8),(phenel0[2]-important["fitwidth"],0.0)],colour=wx.RED,style=wx.DOT)
            henelupper = wxplot.PolyLine([(phenel0[2]+important["fitwidth"],-0.8),(phenel0[2]+important["fitwidth"],0.0)],colour=wx.RED,style=wx.DOT)
            laserline = wxplot.PolyLine([(ptisa0[2],-0.8),(ptisa0[2],0.0)],colour='FOREST GREEN')
            laserlower = wxplot.PolyLine([(ptisa0[2]-important["fitwidth"],-0.8),(ptisa0[2]-important["fitwidth"],0.0)],colour='FOREST GREEN',style=wx.DOT)
            laserupper = wxplot.PolyLine([(ptisa0[2]+important["fitwidth"],-0.8),(ptisa0[2]+important["fitwidth"],0.0)],colour='FOREST GREEN',style=wx.DOT)
            henerline = wxplot.PolyLine([(phener0[2],-0.8),(phener0[2],0.0)],colour=wx.RED)
            henerlower = wxplot.PolyLine([(phener0[2]-important["fitwidth"],-0.8),(phener0[2]-important["fitwidth"],0.0)],colour=wx.RED,style=wx.DOT)
            henerupper = wxplot.PolyLine([(phener0[2]+important["fitwidth"],-0.8),(phener0[2]+important["fitwidth"],0.0)],colour=wx.RED,style=wx.DOT)
            setline = wxplot.PolyLine([(important['henecent']+(important["lockpos"]*important['henewidth']/2),-0.8),(important['henecent']+(important["lockpos"]*important['henewidth']/2),-0.0)],colour=wx.BLACK,style=wx.DOT)
            gc = wxplot.PlotGraphics([mainline,exline,henelline,henellower,henelupper,
                                      laserline,laserlower,laserupper,
                                      henerline,henerlower,henerupper,setline],"Scope","Channel","Voltage")
            canvas.Draw(gc,xAxis=(0,1200),yAxis=(-0.8,0))       # Update the graph
        print "Stopped grabbing scope data"
        
    def grabdata(self):                                         # Function to grab information from the scope
        while True:                                             # Sometimes waveform has errors so this loop catches them
            try:
                waveform = rigol.query(":WAVeform:DATA?")
                waveform = waveform.encode('ascii','ignore')
                waveform = waveform.split(',')
                waveform.pop(0)
                waveform = np.array(waveform)                   # Make data suitable for use in program
                waveform = waveform.astype(np.float)
                x = [i for i in xrange(len(waveform))]          # Create x axis of values for length of waveform
                return x, waveform                              # Return x axis and waveform
            except Exception,e:
                print str(e)
                print "Error importing values, will try again"

class lock_thread(threading.Thread):        # Create Threading Class for locking routine
    def __init__(self):
        threading.Thread.__init__(self)
    def run(self):
        global important, operations, phenel0, ptisa0, phener0
        while True:
            try:
                self.timescale = float(rigol.query("TIM:MAIN:SCAL?"))                           # Find the time scale
                self.chan1scale = float(rigol.query(":CHAN1:SCAL?"))                            # Find the voltage scale
                important["offset"] = round(float(rigol.query("SOURCE1:VOLTAGE:OFFSET?")),4)    # Find current offset
                important["vout"] = round(float(rigol.query("SOURCE2:VOLTAGE:OFFSET?")),4)      # Find current vout
                break
            except Exception,e:
                print str(e)
                print "Problem getting scope info, will try again"
        
        # Create Starting points for fitting peaks: [fwhm, intensity, peak location]
        ptisa0 = np.array([2500.*self.timescale,3.*self.chan1scale,float(laserentry.GetValue())])
        phenel0 = np.array([5000.*self.timescale,3.*self.chan1scale,float(henelentry.GetValue())])
        phener0 = np.array([5000.*self.timescale,3.*self.chan1scale,float(henerentry.GetValue())])
        important["henecent"] = (phenel0[2]+phener0[2])/2                   # Calcualte center of hene peaks
                         
        important["lockpos"] = float(setpoint.GetValue())                   # Set lock position 
        important["Kp"] = 0.0                                               # Set starting feedback constants
        important["Ki"] = 0.0
        important["Kpset"] = float(feedback.GetValue())                     # Set the target feedback constants
        important["Kiset"] = float(integral.GetValue())
        important["rampKp"] = float(rampfeed.GetValue())                    # Set the feedback constant for the ramp offset

        self.offsetwrite = important["offset"]                              # Create var for writing offset
        self.voutwrite = important["vout"]                                  # Create var for writing vout
        self.prevoffs = self.offsetwrite                                    # Set previous write values to know when change has occured
        self.prevvout = self.voutwrite

        self.errors = [0] * 200                                             # Create list for storing errors

        while operations["go"] == True:                                     # Main loop for locking routine
            while True:                                                     # Loop to fit the peaks
                try:
                    ptisa0,self.fitt = self.fitpeak(ptisa0,waveform,xaxis,important["fitwidth"])
                    phenel0,self.fithl = self.fitpeak(phenel0,waveform,xaxis,important["fitwidth"])
                    phener0,self.fithr = self.fitpeak(phener0,waveform,xaxis,important["fitwidth"])
                    important["henecent"] = (phenel0[2]+phener0[2])/2       # Set new hene center
                    important["henewidth"] = phener0[2] - phenel0[2]        # Set new hene width
                    break
                except Exception,e:
                    print str(e)    
                    print "Likely something wrong with the fitting routine"
            
            # Check hene peaks have been correctly fitted
            if self.fithl == 'good' and self.fithr == 'good' and round(phenel0[2],1) != round(phener0[2],1):
                self.middle = (xaxis[0] + xaxis[-1])/2                      # Find centre of x axis
                self.centdiff = important["henecent"] - self.middle         # Find how far away hene peaks are from center
                # Change offset value depending on distance from center of scope.  Normalised by henewidth                       
                important["offset"] += important["rampKp"] * (self.centdiff/important["henewidth"])
                self.offsetwrite = round(important["offset"],3)             # Round to 3 decimal places
                if self.offsetwrite < (self.prevoffs) or self.offsetwrite > (self.prevoffs+0.002):    # Check value has changed
                    try:
                        rigol.write("SOURCE1:VOLTAGE:OFFSET "+str(self.offsetwrite))    # Write new ramp offset value
                        rigol.write("TRIG:EDG:LEV "+str(self.offsetwrite))              # Write new trigger value
                        rigol.write(":CHAN1:OFFS "+str(-1*self.offsetwrite))            # Write new view offset
                        rampoff.SetValue(self.offsetwrite)                              # Update GUI
                        self.prevoffs = self.offsetwrite                                # Save old value for next check
                    except:
                        pass

            if abs(important["Kp"]) < abs(important["Kpset"]):                                 
                important["Kp"] += np.sign(important["Kpset"])*0.05                     # If not at set feedback value increment up
            if abs(important["Kp"]) > abs(important["Kpset"]):
                important["Kp"] = important["Kpset"]                                    # If above the set feedback value change to the set value

            if abs(important["Ki"]) < abs(important["Kiset"]):                                 
                important["Ki"] += np.sign(important["Kiset"])*0.002                    # If not at set feedback value increment up
            if abs(important["Ki"]) > abs(important["Kiset"]):
                important["Ki"] = important["Kiset"]                                    # If above the set feedback value change to the set value
                
            if (self.fitt == 'good' and self.fithl == 'good' and   
                self.fithr == 'good' and                                                # Check fits are good
                round(phenel0[2],1) != round(ptisa0[2],1) and
                round(phenel0[2],1) != round(phener0[2],1) and                          # Check to make sure two peaks haven't fit in the same place
                round(ptisa0[2],1) != round(phener0[2],1)):
                if important["errstate"] == 1:                                          # Update gui if previously had an error
                    main_label.SetLabel("Back Running")
                important["errstate"] = 0                                               # Set error to 0
                # Find laser lock difference from setpoint (relative to hene center)
                important["tisadiff"] = ptisa0[2] - (important["lockpos"]*important["henewidth"]/2) - important["henecent"]
                self.error = (important["tisadiff"]/important["henewidth"])             # Define error value. Normalise to hene width.
                self.errors.append(self.error)                                          # Add error to errors
                self.errors.pop(0)                                                      # Remove old value
                 # Set vout depending on distance from set point.  (Proportional) + (Integral)
                important["vout"] -= (important["Kp"] * self.error) + (important["Ki"] * (sum(self.errors)/len(self.errors)))        
                if abs(important["vout"]) > 5.0:
                    important["vout"] = np.sign(important["vout"])*5.0                  # Limit vout to +-5.0
                self.voutwrite = round(important["vout"],3)                             # Round write voltage to 3 dp
                if self.voutwrite != self.prevvout:
                    try:
                        rigol.write("SOURCE2:VOLTAGE:OFFSET "+str(self.voutwrite))      # If different write the new vout
                        dcoff.SetValue(self.voutwrite)                                  # Update the GUI
                        self.prevvout = self.voutwrite                                  # Save old value for next check
                    except:
                        pass
            else:                                               # If any of the above fails...
                important["errstate"] = 1                       # Set error state as 1
                important["Kp"] = 0.0                           # Set feedback constant to 0
                important["Ki"] = 0.0
                time.sleep(0.5)                                 # Wait 0.5 seconds
                main_label.SetLabel("Locking Error...")         # Update GUI and try again 

            laserentry.SetValue(round(ptisa0[2],1))             # Update laser location in GUI
            henelentry.SetValue(round(phenel0[2],1))            # Update hene left location in GUI                    
            henerentry.SetValue(round(phener0[2],1))            # Update hene right location in GUI
            currentfeedback.SetLabel(str(important["Kp"]))      # Update current feedback constants in GUI
            currentintegral.SetLabel(str(important["Ki"]))
            self.fit_update(henelfit,self.fithl)                # Update laser fit state in GUI
            self.fit_update(henerfit,self.fithr)                # Update hene left fit state in GUI
            self.fit_update(laserfit,self.fitt)                 # Update hene right fit state in GUI
            laser_diff.SetLabel(str(round(important["tisadiff"],3)))         # Update laser difference from lock position in GUI
            if (abs(important["Kp"]) == 0 or abs(important["Kp"]) == 0.0) and (abs(important["Ki"]) == 0 or abs(important["Ki"]) == 0.0):
                second_label.SetLabel("Free Running")           # Update GUI as free running if Kpset = 0
                second_label.SetBackgroundColour(wx.YELLOW)
            elif abs(important["tisadiff"]) < 1.5:
                second_label.SetLabel("Locked")                 # If lock difference is below 1.5 update GUI as locked
                second_label.SetBackgroundColour(wx.GREEN)
            elif abs(important["tisadiff"]) < 2.5:
                second_label.SetLabel("Almost Locked")          # If lock difference is below 2.5 update GUI as almost locked
                second_label.SetBackgroundColour(wx.YELLOW)
            else:
                second_label.SetLabel("Trying to lock...")      # Otherwise update GUI as trying to lock
                second_label.SetBackgroundColour(wx.RED)
            
            for key in mon_info:                                # Update monitoring info
                mon_info[key][1].append(important[mon_info[key][2]])
                if len(mon_info[key][1]) > mon_length:
                    mon_info[key][1].pop(0)
            if monitor != "Wavenumber":                         # If not requesting wavenumber, update the monitor graph
                monitor_data = mon_info[monitor][1]
                monitor_plot = zip(monitor_xaxis[:important["view_length"]],monitor_data[-important["view_length"]:])
                monitor_line = wxplot.PolyLine(monitor_plot)
                mon_gc = wxplot.PlotGraphics([monitor_line],"Monitor","Entry",mon_info[monitor][0])
                mon_canvas.Draw(mon_gc,yAxis=(min(monitor_data[-important["view_length"]:])-0.0001,max(monitor_data[-important["view_length"]:])+0.0001))
            try:
                time.sleep(0.5)                                 # Sleep time before starting loop again (i.e crude refresh rate)
            except:
                pass

    def fit_update(self,label,condition):                       # Function to update GUI goodness of fit
        label.SetLabel(condition)
        if condition == 'good':
            label.SetBackgroundColour(wx.GREEN)
        else:
            label.SetBackgroundColour(wx.RED)

    def residuals(self,p,y,x):                                  # Function that takes fit parameters, waveform data (y) and fit range (x)
        return np.subtract(y,lorentzian(x,p))                   # Returns difference between data (y) and fit
    
    def fitpeak(self,pstart,waveform,xaxis,fitwidth):           # Function to fit peaks and determine if fit is good or not
        lowbound = int(pstart[2])-fitwidth                      # Set lower boundary to search for peak
        if lowbound < 0:                                        # If this value is below zero set boundary to zero
            lowbound = int(0)
        highbound = int(pstart[2])+fitwidth                     # Set higher boundary (this can go beyond scope range)
        fity = deepcopy(waveform[lowbound:highbound])           # Set data range within boundaries
        fity[fity > important["background"]] = important["background"]    # Set values above background (exclusion level) to background
        fitx = deepcopy(xaxis[lowbound:highbound])              # Set the x axis range for creating peak
        try:
            p0 = leastsq(self.residuals,pstart,args=(fity,fitx))# Least squares minimising routine to fit peak
            pend0 = deepcopy(p0[0])                             # Create array of fitted parameters
        except Exception,e:
            print str (e)
            pend0 = [0.0,0.0,0.0]                               # If error set all values to 0 so a bad fit is declared
        
        if (abs(pend0[0]) > 0.1*5000.*self.timescale and        # Judge fit condition. Fwhm and intensities within expected values
            abs(pend0[0]) < 10.*5000.*self.timescale and
            abs(pend0[1]) > 0.1*3.*self.chan1scale and
            abs(pend0[1]) < 10.*3.*self.chan1scale):            
            fit = 'good'                                        # Declare the fit as good
        else:
            fit = 'bad'                                         # Otherwise declare fit as bad
            pend0[0] = 5000.*self.timescale                     # Set the FWHM to a reasonable value
            pend0[1] = 3.*self.chan1scale                       # Set the intensity to a reasonable value
            pend0[2] = pstart[2]                                # Set the location back to where it started looking
        return pend0, fit                                       # Return fitted parameters and goodness of fit

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

        self.CommPort = int(5)                                          # Define the CommPort the wavemeter is on
        print "Using Value:", str(self.CommPort)                        # Print to confirm
        self.DevHandle = self.opendevice(self.CommPort)                 # Open communication with wavemeter on port number
        if self.DevHandle == 1:
            print "Connected to Device"                                 # Connect to device and print values
            important["wavenumber"] = self.getlambda(self.DevHandle)
            self.power = self.getpower(self.DevHandle)
            print "Wavenumber: ", important["wavenumber"]
            print "Power: ", self.power
        else:
            print "Failed to connect to device"                         # Print warning if connection fails
            operations["wave_connected"] = False
            
        while operations["wave_connected"] == True:                     # Initiate loop for grabbing data
            time.sleep(0.1)
            try:
                important["wavenumber"] = self.getlambda(self.DevHandle)# Grab wavenumber
                self.power = self.getpower(self.DevHandle)              # Grab power
            except Exception, e:
                print str(e)
            if operations["wave_connected"] == True:                    # Annoyingly this has to be included in the loop or else it hangs up on writing to the GUI
                wn_label.SetLabel(str(round(important["wavenumber"],3))+" cm-1")    # Set the wavenumber in the GUI
                power_label.SetLabel(str(round(self.power,3))+" mW")                # Set the power in the GUI

                for key in wm_mon_info:                                             # Update monitoring info
                    wm_mon_info[key][1].append(important[wm_mon_info[key][2]])  
                if len(wm_mon_info[key][1]) > mon_length:
                    wm_mon_info[key][1].pop(0)
                if monitor == "Wavenumber":                                         # If requesting wavenumber, update the monitor graph
                    monitor_data = wm_mon_info[monitor][1]
                    monitor_plot = zip(monitor_xaxis[:important["view_length"]],monitor_data[-important["view_length"]:])
                    monitor_line = wxplot.PolyLine(monitor_plot)
                    mon_gc = wxplot.PlotGraphics([monitor_line],"Monitor","Entry",wm_mon_info[monitor][0])
                    mon_canvas.Draw(mon_gc,yAxis=(min(monitor_data[-important["view_length"]:]),max(monitor_data[-important["view_length"]:])))
        print "Closing Wavemeter Device"
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
        trace5 = Histogram(y=[5],name='y density',marker=Marker(color='rgb(148,103,189)'),
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

        while operations["streaming"] == True:      # Streaming loop
            try:
                s1_data['x'] = np.array(xaxis)      # Set s1 x data to scope xaxis
                s1_data['y'] = np.array(waveform)   # Set s1 y data to scope trace
                xtime = datetime.datetime.now().strftime('%y-%m-%d %H:%M:%S.%f')
                s1.write(s1_data)                   # Write scope trace
                if abs(important["tisadiff"]) < 1.0:
                    colorset = 'lightgreen'         # Set colour scheme depending on lock difference
                else:
                    colorset = 'red'
                if important["errstate"] == 0:      # If no error update lock difference in kHz using FSR of etalon = 300 MHz
                    sendytisa = important["tisadiff"]*300000/important["henewidth"]
                s2.write(dict(x=xtime,y=sendytisa,marker=Marker(color=colorset)))
                if important["errstate"] == 1:      # If error show messages...
                    s3.write(dict(text=["LOCKING ERROR",""]))
                    s4.write(dict(text=["",str(round(important["wavenumber"],4))+" cm-1"]))
                elif important["tisadiff"] == lockold and important["tisadiff"] == lockoldold:
                    s3.write(dict(text=["LOCKING ERROR","Locking Program Might Not Be Running"]))
                    s4.write(dict(text=["",str(round(important["wavenumber"],4))+" cm-1"]))
                elif (abs(important["Kp"]) == 0 or abs(important["Kp"]) == 0.0) and (abs(important["Ki"]) == 0 or abs(important["Ki"]) == 0.0):
                    s3.write(dict(text=["",""]))
                    s4.write(dict(text=["Free Running",str(round(important["wavenumber"],4))+" cm-1"]))
                    s5.write(dict(y=sendytisa))
                elif abs(important["tisadiff"]) < 1.0:
                    s3.write(dict(text=["",""]))
                    s4.write(dict(text=["Locked",str(round(important["wavenumber"],4))+" cm-1"]))
                    s5.write(dict(y=sendytisa))
                else:
                    s3.write(dict(text=["",""]))
                    s4.write(dict(text=["",str(round(important["wavenumber"],4))+" cm-1"]))
                    s5.write(dict(y=sendytisa))
                lockoldold = lockold                # Pass on old value
                lockold = important["tisadiff"]     # Save value
                stream_label.SetLabel("Data Sent: "+str(datetime.datetime.now().strftime('%H:%M:%S')))      # Update GUI of time sent
            except Exception,e:
                print str(e)
            time.sleep(1)
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
        self.txt = "Laser Loc,Left Hene Loc,Right Hene Loc,Lock Difference,Locking Voltage,Kp,Lock Error State,Wavenumber,Time\n"
        self.f1.write(self.txt)
        self.f2.write(self.txt)
        log_status.SetLabel(self.filename1)                     # Show local filename in GUI
        while operations["logging"] == True:
            self.txt = (str(ptisa0[2]) + ',' + str(phenel0[2]) + ','
                            + str(phener0[2]) + ',' + str(important["tisadiff"]) + ','
                            + str(important["vout"]) + ',' + str(important["Kp"]) + ','
                            + str(important["errstate"]) + ',' + str(important["wavenumber"]) + ','
                            + str(datetime.datetime.now().strftime('%H:%M:%S')) + "\n")
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
        self.data = zip(xaxis,waveform)
        line = wxplot.PolyLine(self.data)
        gc = wxplot.PlotGraphics([line],"Scope","Channel","Voltage")
        self.Draw(gc,xAxis=(0,1200),yAxis=(-1,0))

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

    ####### Main operation functions

    def start_lock(self,event):                                 # Function to start locking thread
        global operations
        if operations["connected"] == False:                    # Check that the scope has been connected
            main_label.SetLabel("Scope not connected!")
        elif operations["go"] == True:
            main_label.SetLabel("Program Already Running")      # Check program isn't already running
        else:
            operations["go"] = True                             # If not then set go to true
            main_label.SetLabel("Starting Locking Program")     # Update GUI
            thread1 = lock_thread()                             # Create thread
            thread1.start()                                     # Start thread

    def start_wavemeter(self,event):                            # Function to start wavemeter thread
        global operations
        if operations["wave_connected"] == True:
            main_label.SetLabel("Wavemeter already connected")
        else:
            operations["wave_connected"] = True
            main_label.SetLabel("Starting Wavemeter Program" )
            self.wave_label.SetLabel("Connected")
            thread2 = wavemeter_thread()                
            thread2.start()  

    def start_stream(self,event):                               # Function to start streaming thread
        global operations
        if operations["streaming"] == True:
            main_label.SetLabel("Already streaming data" )
        else:
            operations["streaming"] = True
            main_label.SetLabel("Starting Streaming Program" )
            self.stream_status.SetLabel("Streaming")
            thread3 = streaming_thread()
            thread3.start()

    def start_logging(self,event):                              # Function to start logging data thread
        global operations
        if operations["logging"] == True:
            main_label.SetLabel("Already logging data")
        else:
            operations["logging"] = True
            main_label.SetLabel("Logging Data")
            thread4 = logging_thread()                    
            thread4.start()       

    def stop_lock(self,event):                                  # Function to stop locking
        global operations, important
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

    def quit_program(self,event):                               # Function to quit whole program safely
        print "Closing all connections, will quit in 5 seconds!"
        global operations
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

    def connect_scope(self,event):                              # Function to connect scope
        global rigol, operations
        if operations["connected"] == True:                     # Check scope is not already connected
            main_label.SetLabel("Scope is already connected")
        else:
            try:
                main_label.SetLabel("Connecting...")
                rm = visa.ResourceManager()                                                 # Create a resource manager to open up communications with the Rigol Device
                devices = rm.list_resources()                                               # Find resources to connect to
                print "These are the devices connected to the computer:\n"
                for i in xrange(len(devices)):
                    print i, devices[i]                                                     # Print list of resources
                rigol = rm.open_resource(devices[int(self.scope_num.GetValue())])           # Connect to device number stated in GUI
                print "\nConnected with device", rigol.query("*IDN?")
                rigol.write(":WAVeform:SOURce CHANnel3")                                    # Make sure scope is set up to outfut waveform from Channel 3
                rampoff.SetValue(round(float(rigol.query("SOURCE1:VOLTAGE:OFFSET?")),4))
                dcoff.SetValue(round(float(rigol.query("SOURCE2:VOLTAGE:OFFSET?")),4))      # Set ramp and vout values in GUI depending on how the scope is set
                self.scope_label.SetLabel("Connected with device")                          # Change GUI to show scope is connected
                self.scope_label.SetBackgroundColour(wx.GREEN)
                main_label.SetLabel("Connected")
                operations["connected"] = True
                operations["grabbing"] = True
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
        global operations
        if operations["connected"] == False:                    # Check that the scope has been connected
            main_label.SetLabel("Scope not connected") 
        elif operations["go"] == True:
            main_label.SetLabel("Not when locking!")            # Check program isn't already running
        else:
            operations["grabbing"] = False
            main_label.SetLabel("Resetting scope")                               
            rigol.write("*RST")                                 # Reset the scope so that conditions for locking can be set
            time.sleep(5)                                       # Sleep for 5 seconds to allow scope to restart before set up
            rigol.write(":WAVeform:FORMat ASC")                 # Set up scope
            rigol.write(":CHAN1:DISP ON")
            rigol.write(":CHAN1:PROB 1")
            rigol.write(":CHAN1:SCAL 0.1")
            rigol.write(":CHAN1:OFFS -0.5")
            rigol.write(":CHAN2:DISP ON")
            rigol.write(":CHAN2:PROB 1")
            rigol.write(":CHAN2:SCAL 0.5")
            rigol.write(":CHAN3:DISP ON")
            rigol.write(":CHAN3:PROB 1")
            rigol.write(":CHAN3:SCAL 0.075")
            rigol.write(":CHAN3:OFFS 0.3")
            rigol.write("TIM:MAIN:SCAL 0.002")
            rigol.write("TRIG:EDG:LEV 0.5")
            rigol.write("TRIG:EDG:SLOP POS")
            
            rigol.write(":SOUR1:FUNC RAMP")                     # Set up the Ramp
            rigol.write(":SOUR1:VOLT 2.5")
            rigol.write(":SOUR1:VOLT:OFFS 0.5")
            rampoff.SetValue(0.5)
            rigol.write(":SOUR1:FREQ 15")
            rigol.write(":SOUR1:OUTP 1")
            
            time.sleep(1)                                       # Trick to get 50% symmetry ramp, scope doesnt do it if set straight at 50
            rigol.write(":SOUR1:FUNC:RAMP:SYMM 49")
            time.sleep(1)
            rigol.write(":SOUR1:FUNC:RAMP:SYMM 50")
            
            rigol.write(":SOUR2:FUNC DC")                       # Set up the TiSa voltage feedback
            rigol.write(":SOUR2:VOLT:OFFS 0.0")
            dcoff.SetValue(0.0)
            rigol.write(":SOUR2:OUTP 1")
            rigol.write(":WAVeform:SOURce CHANnel3")            # Make sure scope is set up to outfut waveform from Channel 3

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
            operations["grabbing"] = 1
            thread5 = grabbing_thread()                         # Create thread
            thread5.start()                                     # Start thread
            main_label.SetLabel("Ready to go!") 

    def move_ramp(self,event):
        global important
        if operations["connected"] == False:                                 # Check that the scope has been connected
            main_label.SetLabel("Scope not connected!")
        else:
            if operations["go"] == False:
                rigol.write("TRIG:EDG:LEV "+str(rampoff.GetValue()))         # Set trigger level on scope
                rigol.write(":SOUR1:VOLT:OFFS "+str(rampoff.GetValue()))     # Set ramp offset on scope
                rigol.write(":CHAN1:OFFS "+str(-1*rampoff.GetValue()))
            important["offset"] = float(rampoff.GetValue())                  # Change offset for locking program

    def move_tisa(self,event):
        global important
        if operations["connected"] == False:                                 # Check that the scope has been connected
            main_label.SetLabel("Scope not connected!")
        else:
            if operations["go"] == False:
                rigol.write(":SOUR2:VOLT:OFFS "+str(dcoff.GetValue()))       # Set locking voltage on scope
            important["vout"] = float(dcoff.GetValue())                      # Change vout for locking program

    def print_help(self,event):
        print "If list of network devices does not load there is a communications error.  try these steps in order before contacting T Procter for help:"
        print "Even if you are not using these units, if they are part of the system they can still hold up Ni-VISA"
        print "This can be checked by opening up NI MAX" 
        print "\n1. Restart Erthernet/GPIB adapter by DAQ computers"
        print "2. Restart Agilent unit"
        print "3. On Rigol scope: Go to Utility menue and change USB device from Computer to PictBridge and back again"
        print "\n4. Restart this Windows Machine"

    ####### Functions for changing variables with widgets

    def change_feedback(self,event):
        global important
        important["Kpset"] = float(event.GetValue())                        # Change the set feedback constant

    def change_rampfeed(self,event):
        global important
        important["rampKp"] = float(rampfeed.GetValue())                    # Change the ramp feedback constant

    def change_setpoint(self,event):
        global important
        important["Kp"], important["Ki"] = 0.0, 0.0           
        important["lockpos"] = float(setpoint.GetValue())                   # Change set point and set Kp to 0 for smoother shift

    def change_integral(self,event):
        global important          
        important["Kiset"] = float(integral.GetValue())                     # Change Ki

    def change_background(self,event):
        global important
        important["background"] = float(excludeentry.GetValue())            # Change exclusion level

    def update_params(self,event):
        global phenel0, ptisa0, phener0                         # Update fitting parameters with values from GUI
        phenel0[2] = float(henelentry.GetValue())
        ptisa0[2] = float(laserentry.GetValue())
        phener0[2] = float(henerentry.GetValue())
        important['henewidth']=phener0[2]-phenel0[2]
        important['henecent']=(phener0[2]+phenel0[2])/2

    def OnSelect(self, e):                                                  # Change what monitor graph is showing
        global monitor
        monitor = e.GetString()
        self.mon_name.SetLabel(monitor)

    def View_Select(self,e):
        global important
        important["view_length"] = int(e.GetString())

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
        sizer.Add(spin,(x,y),(h,w),wx.ALL|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND,3)     # Add control to grid
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

    ####### Create GUI Function

    def initialize(self,id):
        global panel, important
        panel = wx.Panel(self,id)

        # Create Sizers
        self.mainSizer = wx.FlexGridSizer(1,2,1,1)
        self.mainSizer.AddGrowableCol(1,1)
        self.graphSizer = wx.GridSizer(2,1,1,1)
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
        
        self.scopebox = wx.StaticBox(panel, label='Rigol Scope')
        self.scopeBoxSizer = wx.StaticBoxSizer(self.scopebox,wx.VERTICAL)
        self.scopeSizer = wx.GridBagSizer(hgap=1, vgap=1)

        self.settingsbox = wx.StaticBox(panel, label='Locking Settings')
        self.settingsBoxSizer = wx.StaticBoxSizer(self.settingsbox,wx.VERTICAL)
        self.settingsSizer = wx.GridBagSizer(hgap=1, vgap=1)

        self.locksbox = wx.StaticBox(panel, label='Locking Information')
        self.lockBoxSizer = wx.StaticBoxSizer(self.locksbox,wx.VERTICAL)
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
        self.scope_label = self.create_label(id,self.scopeSizer,'   Scope Not Connected   ',wx.RED,wx.BLACK,0,0,1,2,12)
        self.device_label = self.create_label(id,self.scopeSizer,'VISA Device Number:',wx.NullColour,wx.BLACK,1,0,1,1,10)
        self.scope_num = self.create_spinctrl(id,self.scopeSizer,'0',1,1,1,1,None,1,0,10)
        self.connect_button = self.create_button(id,self.scopeSizer,'  Connect Scope  ',2,0,1,2,12,self.connect_scope)
        self.reset_button = self.create_button(id,self.scopeSizer,'  Reset Scope  ',4,0,1,2,10,self.reset_scope)
        self.scopeBoxSizer.Add(self.scopeSizer,0,wx.CENTER)

        # Fill Settings Sizer
        global rampoff, rampfeed, excludeentry, feedback, currentfeedback, integral, currentintegral
        self.ramplabel = self.create_label(id,self.settingsSizer,'Ramp Offset:',wx.NullColour,wx.BLACK,0,0,1,1,10)
        rampoff = self.create_spinctrl(id,self.settingsSizer,str(important["offset"]),0,1,1,1,self.move_ramp,0.01,-2,2)
        self.rampfeed_label = self.create_label(id,self.settingsSizer,'Ramp Feedback:',wx.NullColour,wx.BLACK,1,0,1,1,10)
        rampfeed = self.create_spinctrl(id,self.settingsSizer,str(important["rampKp"]),1,1,1,1,self.change_rampfeed,0.001,0,0.5)
        self.excludelabel = self.create_label(id,self.settingsSizer,'Exclusion Level:',wx.NullColour,wx.BLACK,2,0,1,1,10)
        excludeentry = self.create_spinctrl(id,self.settingsSizer,str(important["background"]),2,1,1,1,self.change_background,0.01,-1,0)
        self.feedbacklabel = self.create_label(id,self.settingsSizer,'Kp Set:',wx.NullColour,wx.BLACK,3,0,1,1,10)
        feedback = self.create_spinctrl(id,self.settingsSizer,str(important["Kpset"]),3,1,1,1,self.change_feedback,0.1,0,2)
        self.feedbackcurrentlabel = self.create_label(id,self.settingsSizer,'Current Kp:',wx.NullColour,wx.BLACK,4,0,1,1,10)
        currentfeedback = self.create_label(id,self.settingsSizer,str(important["Kp"]),wx.NullColour,wx.BLACK,4,1,1,1,10)
        self.integral_label = self.create_label(id,self.settingsSizer,'Ki Set:',wx.NullColour,wx.BLACK,5,0,1,1,10)
        integral = self.create_spinctrl(id,self.settingsSizer,str(important["Kiset"]),5,1,1,1,self.change_integral,0.1,0,1)
        self.currentintegral_label = self.create_label(id,self.settingsSizer,'Current Ki:',wx.NullColour,wx.BLACK,6,0,1,1,10)
        currentintegral = self.create_label(id,self.settingsSizer,str(important["Ki"]),wx.NullColour,wx.BLACK,6,1,1,1,10)      
        self.settingsBoxSizer.Add(self.settingsSizer,0,wx.CENTER)

        # Fill Lock Sizer
        global henelentry, laserentry, henerentry, dcoff, setpoint, laser_diff, henelfit, laserfit, henerfit
        self.henellabel = self.create_label(id,self.lockSizer,'Left Hene Location:',wx.NullColour,wx.BLACK,0,0,1,2,10)
        self.laserlabel = self.create_label(id,self.lockSizer,'Laser Location:',wx.NullColour,wx.BLACK,0,2,1,2,10)
        self.henerlabel = self.create_label(id,self.lockSizer,'Right Hene Location:',wx.NullColour,wx.BLACK,0,4,1,2,10)
        henelentry = self.create_slider(id,self.lockSizer,phenel0[2],1,0,1,2,self.update_params,0,1200)
        laserentry = self.create_slider(id,self.lockSizer,ptisa0[2],1,2,1,2,self.update_params,0,1200)
        henerentry = self.create_slider(id,self.lockSizer,phener0[2],1,4,1,2,self.update_params,0,1200)
        henelfit = self.create_label(id,self.lockSizer,'Not Yet',wx.YELLOW,wx.BLACK,2,0,1,2,12)
        laserfit = self.create_label(id,self.lockSizer,'Not Yet',wx.YELLOW,wx.BLACK,2,2,1,2,12)
        henerfit = self.create_label(id,self.lockSizer,'Not Yet',wx.YELLOW,wx.BLACK,2,4,1,2,12)

        self.dclabel = self.create_label(id,self.lockSizer,'V Out:',wx.NullColour,wx.BLACK,4,0,1,1,10)
        dcoff = self.create_spinctrl(id,self.lockSizer,str(important["vout"]),4,1,1,1,self.move_tisa,0.01,-2,2)
        self.setpointlabel = self.create_label(id,self.lockSizer,'Set Point:',wx.NullColour,wx.BLACK,4,2,1,1,10)
        setpoint = self.create_spinctrl(id,self.lockSizer,str(important["lockpos"]),4,3,1,1,self.change_setpoint,0.01,-0.9,0.9)
        self.diff_label = self.create_label(id,self.lockSizer,'Lock Diff:',wx.NullColour,wx.BLACK,4,4,1,1,10)
        laser_diff = self.create_label(id,self.lockSizer,str(important["tisadiff"]),wx.NullColour,wx.BLACK,4,5,1,1,10)
        self.lockBoxSizer.Add(self.lockSizer,0,wx.CENTER)

        # Fill Control Sizer
        self.start_button = wx.Button(panel,id,label='START LOCK')
        self.Bind(wx.EVT_BUTTON,self.start_lock,self.start_button)
        self.stop_button = wx.Button(panel,id,label='STOP LOCK')
        self.Bind(wx.EVT_BUTTON,self.stop_lock,self.stop_button)
        self.controlSizer.AddMany([(self.start_button,0,wx.EXPAND),(self.stop_button,0,wx.EXPAND)])

        # Fill Wavemeter Sizer
        global wn_label, power_label
        self.wavemeter_button = wx.Button(panel,id,label='Connect Wavemeter')
        self.Bind(wx.EVT_BUTTON,self.start_wavemeter,self.wavemeter_button)
        font = wx.Font(20,wx.DECORATIVE,wx.NORMAL,wx.NORMAL)
        wn_label = wx.StaticText(panel,id,label='Wave Number:')
        wn_label.SetFont(font)
        self.wave_label = wx.StaticText(panel,id,label='Not Connected')
        power_label = wx.StaticText(panel,id,label='Power:')
        power_label.SetFont(font)
        self.waveSizer.AddMany([(self.wavemeter_button,0,wx.ALIGN_CENTER|wx.EXPAND),(wn_label,0,wx.ALIGN_CENTER),
                                (self.wave_label,0,wx.ALIGN_CENTER),(power_label,0,wx.ALIGN_CENTER)])

        # Fill Log Sizer
        global log_status, log_label
        self.logging_button = wx.Button(panel,id,label="Log Data")
        self.Bind(wx.EVT_BUTTON,self.start_logging,self.logging_button)
        log_status = wx.StaticText(panel,id,label='Not Logging')
        log_label = wx.StaticText(panel,id,label='Not Logging')
        self.logSizer.AddMany([(self.logging_button,0,wx.ALIGN_CENTER|wx.EXPAND),(log_status,0,wx.LEFT|wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL,10),
                                (log_label,0,wx.ALIGN_CENTER)])

        # Fill Stream Sizer
        global stream_label
        self.stream_button = wx.Button(panel,id,label="Stream to Web")
        self.Bind(wx.EVT_BUTTON,self.start_stream,self.stream_button)
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

        # Fill Monitor sizer
        mon_options = ['Lock Difference','HeNe Width',
        'HeNe Center','V Out','Ramp Offset','Error State','Kp','Ki','Wavenumber']
        view_options = ['100','500','1000']
        self.monlabelSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.mon_label = wx.StaticText(panel,id,label='Currently Monitoring:')
        self.monlabelSizer.Add(self.mon_label,0,wx.LEFT,10)
        self.mon_name = wx.StaticText(panel,id,label='Lock Difference')
        self.monlabelSizer.Add(self.mon_name,0,wx.LEFT,10)
        self.combo = wx.ComboBox(panel,id,value=mon_options[0],choices=mon_options,style=wx.CB_READONLY)
        self.combo.Bind(wx.EVT_COMBOBOX,self.OnSelect)
        self.monlabelSizer.Add(self.combo,0,wx.ALIGN_RIGHT|wx.LEFT,100)
        self.view_length = wx.ComboBox(panel,id,value=view_options[1],choices=view_options,style=wx.CB_READONLY)
        self.view_length.Bind(wx.EVT_COMBOBOX,self.View_Select)
        self.monlabelSizer.Add(self.view_length,0,wx.ALIGN_RIGHT|wx.LEFT,100)
        global mon_canvas
        mon_canvas = PlotMonitor(panel,id)

        self.monitorSizer.Add(self.monlabelSizer,0,wx.ALL|wx.LEFT,5)
        self.monitorSizer.Add(mon_canvas,0,wx.ALL|wx.EXPAND,5)

        # Fill Graph sizer
        global canvas
        canvas = PlotGraph(panel,id)
        self.graphSizer.Add(canvas,0,wx.ALL|wx.EXPAND,5)
        self.graphSizer.Add(self.monitorSizer,0,wx.ALL|wx.EXPAND,5)

        # Fill Main sizer
        self.mainSizer.Add(self.topSizer,0,wx.ALL|wx.EXPAND,5)
        self.mainSizer.Add(self.graphSizer,0,wx.ALL|wx.EXPAND,5)

        # Set panel attributes
        panel.SetSizer(self.mainSizer)
        self.SetSizeHints(1200,770,2000,770)
        self.topSizer.Fit(self)
        
if __name__  == "__main__":
    app = wx.App()                                              # Creat wx python App
    frame = laser_lock_wx(None,-1,"Laser Lock").Show()          # Create frame using class built for GUI
    app.MainLoop()                                              # Start loop for GUI interaction
