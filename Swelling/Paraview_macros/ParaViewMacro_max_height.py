# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import numpy as np

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

animationScene = GetAnimationScene()
tk = GetTimeKeeper()
timesteps = tk.TimestepValues

time = []
array = []
for val in timesteps:
	indata = GetActiveSource()
	bounds = indata.GetDataInformation().GetBounds()
	yHeight = bounds[3]
	array.append(yHeight)
	animationScene.GoToNext()

	time.append(val)

data_export = np.column_stack((time, array))
np.savetxt("./output.csv", data_export)
#np.savetxt("Users/idaang/Dropbox/cartilage_mechanics_IDA/FEniCS/Swelling/output.txt", data_export)
