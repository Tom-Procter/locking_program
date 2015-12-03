# Python GUI for viewing log data

import wx
import wx.lib.plot as wxplot

import pandas

def get_path(wildcard):
	app = wx.App(None)
	style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
	dialog = wx.FileDialog(None, 'Open', wildcard=wildcard, style=style)
	if dialog.ShowModal() == wx.ID_OK:
		path = dialog.GetPath()
	else:
		path=None
	dialog.Destroy()
	return path

filename = get_path('*.csv')

global info
info = pandas.read_csv(filename)
options = info.columns

class Plot(wxplot.PlotCanvas):
    def __init__(self,parent,id):
        wxplot.PlotCanvas.__init__(self,parent,id,style=wx.BORDER_NONE,size=(400,400))
    	tag = options[0]
    	data = info[tag]
    	xaxis = list(xrange(len(data)))
    	trace = zip(xaxis,data)
        plot_line = wxplot.PolyLine(trace)
        mon_gc = wxplot.PlotGraphics([plot_line],tag,"Entry","Value")
        self.Draw(mon_gc)

class log_plot_wx(wx.Frame):
    def __init__(self,parent,id,title):
        wx.Frame.__init__(self,parent,id,title)
        self.parent = parent
        self.initialize(id)   

    def initialize(self,id): 
        global panel, canvas
        panel = wx.Panel(self,id)

        self.mainSizer = wx.FlexGridSizer(3,1,1,1)
        self.mainSizer.AddGrowableCol(0,1)
        self.mainSizer.AddGrowableRow(1,1)
        self.settingsSizer = wx.BoxSizer(wx.VERTICAL)


        # Create Title
        font = wx.Font(16,wx.DECORATIVE,wx.NORMAL,wx.NORMAL)
        title = wx.StaticText(panel,id,'Laser Log Viewer')
        title.SetFont(font)

        canvas = Plot(panel,id)

        self.combo = wx.ComboBox(panel,id,value=options[0],choices=options,style=wx.CB_READONLY)
        self.combo.Bind(wx.EVT_COMBOBOX,self.update_plot)

        self.settingsSizer.Add(self.combo,0,wx.ALIGN_CENTER,5)

        self.mainSizer.Add(title,0,wx.ALIGN_CENTER)
        self.mainSizer.Add(canvas,0,wx.EXPAND,5)
        self.mainSizer.Add(self.settingsSizer,0,wx.ALL|wx.EXPAND,5)

        panel.SetSizer(self.mainSizer)

        self.SetSizeHints(700,600,1000,1000)
        self.mainSizer.Fit(self)

    def update_plot(self,e):
    	tag = e.GetString()
    	data = info[tag]
    	xaxis = list(xrange(len(data)))
    	trace = zip(xaxis,data)
    	plot_line = wxplot.PolyLine(trace)
    	gc = wxplot.PlotGraphics([plot_line],tag,"Entry","Value")
    	canvas.Draw(gc)

if __name__  == "__main__":
    app = wx.App()                
    frame = log_plot_wx(None,-1,"Laser Log Plot").Show()  
    app.MainLoop()   