# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import os

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#For test
path_models = ['/scratch/usr/bbkponsm/model/globalscale/sphere3d/R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003/',
        '/scratch/usr/bbkponsm/model/globalscale/sphere3d/R01f_SRW01_1em15_Rodinia_2GPa_Mantle_C10MPa_f005_LR/',
        '/scratch/usr/bbkponsm/model/globalscale/sphere3d/P01a_Pangea_1GPa_Mantle_C40MPa_LR/']

output_model_folder = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/test_sphere_simple/'

try:
    os.makedirs(output_model_folder)
except FileExistsError:
    print("model repository for postprocess already exists")
    # directory already exists
    pass

# create a new 'Sphere'
sphere1 = Sphere(registrationName='Sphere1')

# Properties modified on sphere1
sphere1.Radius = 6371000.0

UpdatePipeline(time=0.0, proxy=sphere1)

# set active source
SetActiveSource(sphere1)


# save data
SaveData(output_model_folder+'sphere_simple.csv', proxy=sphere1, ChooseArraysToWrite=1,
    PointDataArrays=['Normals'])
