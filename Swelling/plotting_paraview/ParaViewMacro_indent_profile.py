# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import numpy as np

#### USER PARAMETERS ####
l0 = 1.4
TimeStep = 20
SampleRate = 15
FolderName = "Data"

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
penaltyxdmf = GetActiveSource()

# set active source
SetActiveSource(penaltyxdmf)

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(Input=penaltyxdmf)
warpByVector1.Vectors = ['POINTS', 'Displacement']

# Set active Source
SetActiveSource(warpByVector1)

# Get active view
renderView1 = GetActiveViewOrCreate("RenderView")

# Set how many samples we want across the indentation profile
sample = np.linspace(0, l0, num=SampleRate)

# Used for file names to track data files
count = 1
for val in sample:

	# create a new 'Plot Over Line'
	PlotOverLineVal = PlotOverLine(Input=warpByVector1, Source='High Resolution Line Source')
	# init the 'High Resolution Line Source' selected for 'Source'
	PlotOverLineVal.Source.Point1 = [0.7, 0.0, val]
	PlotOverLineVal.Source.Point2 = [0.7, 1.5, val]

	# Show data in view
	plotOverLineDisplay = Show(PlotOverLineVal, renderView1)

	# Create a new 'Line Chart View'
	lineChartView = CreateView('XYChartView')

	# get layout
	layout = GetLayout()

	# place view in the layout
	layout.AssignView(2, lineChartView)

	# show data in view
	plotOverLine1Display_1 = Show(PlotOverLineVal, lineChartView)

	# update the view to ensure updated data information
	lineChartView.Update()

	# set active view
	name = "PlotOverLine1"
	plotOverLine = FindSource(name)
	SetActiveSource(plotOverLine)

	# get active source.
	SaveData("./" + FolderName + "/" + str(TimeStep) + "_" + str(count) + ".csv", proxy=plotOverLine)
	#np.savetxt("Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/" + str(count) + ".txt", PlotOverLineVal)

	# Close in Pipeline Browser
	Delete(plotOverLine)
	del plotOverLine

	# Find lineChartView1
	layoutName = "LineChartView1"	# Default name for first LineChart opened
	# get layout
	layout = FindViewOrCreate(layoutName, viewtype='XYChartView')
	SetActiveView(layout)
	# Close lineChartView1
	Delete(layout)
	del layout
	# Close empty frame
	layout = GetLayoutByName("Layout #1")
	layout.Collapse(5)

	count += 1