# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
solutionpvd = PVDReader(registrationName='solution.pvd', FileName='/Volumes/Jerry/global_models_3d/R01e_Rodinia_2GPa_Mantle_C20MPa_f003_LR/solution.pvd')
solutionpvd.PointArrays = ['velocity', 'stress', 'shear_stress', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'p', 'T', 'continent', 'asthenosphere', 'strain_rate', 'nonadiabatic_temperature', 'depth', 'vertical_heat_flux', 'heat_flux_map', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity', 'thermal_conductivity', 'thermal_diffusivity', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'plastic_yielding', 'v_r', 'v_phi', 'v_theta', 'principal_stress_1', 'principal_stress_2', 'principal_stress_3']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on solutionpvd
solutionpvd.PointArrays = ['T', 'density', 'depth', 'viscosity']

UpdatePipeline(time=0.0, proxy=solutionpvd)

# Properties modified on animationScene1
animationScene1.Stride = 20

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=solutionpvd)
contour1.ContourBy = ['POINTS', 'T']
contour1.Isosurfaces = [2133.2440185546875]
contour1.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour1
contour1.ContourBy = ['POINTS', 'depth']
contour1.Isosurfaces = [100000.0, 298571.4285714286, 497142.85714285716, 695714.2857142857, 894285.7142857143, 1092857.142857143, 1291428.5714285714, 1490000.0, 1688571.4285714286, 1887142.8571428573, 2085714.285714286, 2284285.7142857146, 2482857.1428571427, 2681428.5714285714, 2880000.0]

UpdatePipeline(time=0.0, proxy=contour1)

# save data
SaveData('/Volumes/Jerry/global_models_3d_extract/R01e_Rodinia_2GPa_Mantle_C20MPa_f003_LR/depths_contours.csv', proxy=contour1, WriteTimeSteps=1,
    WriteTimeStepsSeparately=1,
    ChooseArraysToWrite=1,
    PointDataArrays=['T', 'density', 'depth', 'viscosity'],
    AddTime=1)