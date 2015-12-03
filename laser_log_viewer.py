import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
plt.style.use('ggplot')

xaxis = []
data1 = []
data2 = []
data3 = []
data4 = []

datanames = ["Agilent_1","Agilent_2",							#0,1
			"Error State","HeNe Center",						#2,3
			"HeNe Left Location","HeNe Right Location",			#4,5
			"HeNe Sideband Location","HeNe Width",				#6,7
			"Ki","Kp",											#8,9
			"Laser Location","Lock Difference",					#10,11
			"Ramp Amplitude","Ramp Kp",							#12,13
			"Ramp Offset","V Out",								#14,15
			"Wavenumber"]										#16

data1name = datanames[12]
data2name = datanames[14]
data3name = datanames[15]
data4name = datanames[16]

allFiles = glob.glob("Data/Nov_2015/*.csv")
for file_ in allFiles:
	print "Loading File:"
	print file_
	df = pd.read_csv(file_,index_col=None, header=0)
	df["Time Stamp"] = pd.to_datetime(df["Time Stamp"])
	df["month"] = df["Time Stamp"].apply(lambda x: x.month)
	df["day"] = df["Time Stamp"].apply(lambda x: x.day)
	df["year"] = df["Time Stamp"].apply(lambda x: x.year)
	df["hour"] = df["Time Stamp"].apply(lambda x: x.hour)
	df["minute"] = df["Time Stamp"].apply(lambda x: x.minute)
	df["second"] = df["Time Stamp"].apply(lambda x: x.second)

	for i,val in enumerate(df["year"]):
		xaxis.append(datetime.datetime(int(df["year"][i]),int(df["month"][i]),int(df["day"][i]),int(df["hour"][i]),int(df["minute"][i]),int(df["second"][i])))
		try:
			data1.append(float(df[data1name][i]))
		except:
			data1.append(0)
		try:
			data2.append(float(df[data2name][i]))
		except:
			data2.append(0)		
		try:
			data3.append(float(df[data3name][i]))
		except:
			data3.append(0)
		try:
			data4.append(float(df[data4name][i]))
		except:
			data4.append(0)

f, axarr = plt.subplots(4, sharex=True)
axarr[0].plot_date(xaxis, data1,'b-')
axarr[0].set_title(data1name)
axarr[1].plot_date(xaxis, data2,'b-')
axarr[1].set_title(data2name)
axarr[2].plot_date(xaxis, data3,'b-')
axarr[2].set_title(data3name)
axarr[3].plot_date(xaxis, data4,'b-')
axarr[3].set_title(data4name)

plt.show()

