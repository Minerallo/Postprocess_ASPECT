# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
solution_surfacepvd = PVDReader(registrationName='solution_surface.pvd', FileName='/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/solution_surface.pvd')
solution_surfacepvd.PointArrays = ['surface_strain_rate_tensor', 'surface_stress']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'PVD Reader'
solutionpvd = PVDReader(registrationName='solution.pvd', FileName='/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/solution.pvd')
solutionpvd.PointArrays = ['velocity', 'stress', 'shear_stress', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'p', 'T', 'continent', 'asthenosphere', 'strain_rate', 'nonadiabatic_temperature', 'depth', 'vertical_heat_flux', 'heat_flux_map', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity', 'thermal_conductivity', 'thermal_diffusivity', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'plastic_yielding', 'v_r', 'v_phi', 'v_theta', 'principal_stress_1', 'principal_stress_2', 'principal_stress_3']

# find source
solutionpvd_1 = FindSource('solution.pvd')

# set active source
SetActiveSource(solutionpvd_1)

# set active source
SetActiveSource(solutionpvd)

# find source
solution_surfacepvd_1 = FindSource('solution_surface.pvd')

# Properties modified on solutionpvd
solutionpvd.PointArrays = ['T', 'continent', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'density', 'depth', 'nonadiabatic_temperature', 'p', 'plastic_yielding', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'specific_heat', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity']

UpdatePipeline(time=0.0, proxy=solutionpvd)

# Properties modified on solution_surfacepvd
solution_surfacepvd.PointArrays = []

UpdatePipeline(time=0.0, proxy=solution_surfacepvd)

# set active source
SetActiveSource(solution_surfacepvd)

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=solution_surfacepvd)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = [None, '']

# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=clip2.ClipType)

# Properties modified on clip2
clip2.ClipType = 'Sphere'
clip2.Scalars = ['POINTS', '']
clip2.Invert = 0

# Properties modified on clip2.ClipType
clip2.ClipType.Radius = 5000000.0

UpdatePipeline(time=0.0, proxy=clip2)


# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=clip2.ClipType)



# set active source
SetActiveSource(clip2)


# toggle interactive widget visibility (only when running from the GUI)
ShowInteractiveWidgets(proxy=clip2.ClipType)

# set active source
SetActiveSource(solutionpvd)

# toggle interactive widget visibility (only when running from the GUI)
HideInteractiveWidgets(proxy=clip2.ClipType)

# set active source
SetActiveSource(solutionpvd)

# create a new 'Gradient'
gradient2 = Gradient(registrationName='Gradient2', Input=solutionpvd)
gradient2.ScalarArray = ['POINTS', 'T']

# Properties modified on gradient2
gradient2.ScalarArray = ['POINTS', 'velocity']
gradient2.ComputeDivergence = 1
gradient2.ComputeVorticity = 1
gradient2.ComputeQCriterion = 1

UpdatePipeline(time=0.0, proxy=gradient2)


# set active source
SetActiveSource(gradient2)


# create a new 'Contour'
contour2 = Contour(registrationName='Contour2', Input=gradient2)
contour2.ContourBy = ['POINTS', 'Divergence']
contour2.Isosurfaces = [3.7377972148533445e-06]
contour2.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour2
contour2.ContourBy = ['POINTS', 'T']

UpdatePipeline(time=0.0, proxy=contour2)

# Properties modified on contour2
contour2.Isosurfaces = [1573.0]

# set active source
SetActiveSource(contour2)

# set active source
SetActiveSource(contour2)

# find source
extract2 = FindSource('extract2')

# set active source
SetActiveSource(extract2)

# find source
contour1 = FindSource('Contour1')

# set active source
SetActiveSource(contour1)

# set active source
SetActiveSource(extract2)

# set active source
SetActiveSource(contour2)

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=contour2)
calculator1.Function = ''

# set active source
SetActiveSource(extract2)

# set active source
SetActiveSource(contour2)

# set active source
SetActiveSource(calculator1)

# Properties modified on calculator1
calculator1.ResultArrayName = 'oceanic_age'
calculator1.Function = '((depth / 2.32)^2) / (1e-6 * 31536000)/1e6'

UpdatePipeline(time=0.0, proxy=calculator1)

# rename source object
RenameSource('extract2', calculator1)

# find source
extract3 = FindSource('extract3')

# set active source
SetActiveSource(extract3)

# set active source
SetActiveSource(gradient2)

# set active source
SetActiveSource(extract3)

# set active source
SetActiveSource(gradient2)

# create a new 'Contour'
contour3 = Contour(registrationName='Contour3', Input=gradient2)
contour3.ContourBy = ['POINTS', 'Divergence']
contour3.Isosurfaces = [3.7377972148533445e-06]
contour3.PointMergeMethod = 'Uniform Binning'

# set active source
SetActiveSource(extract3)

# set active source
SetActiveSource(contour3)

# set active source
SetActiveSource(extract3)

# set active source
SetActiveSource(contour3)

# Properties modified on contour3
contour3.ContourBy = ['POINTS', 'depth']
contour3.Isosurfaces = [440000.0, 660000.0, 1000000.0, 2000000.0, 2600000.0]

UpdatePipeline(time=0.0, proxy=contour3)

# rename source object
RenameSource('extract3', contour3)

# find source
extract1 = FindSource('extract1')

# set active source
SetActiveSource(extract1)

# set active source
SetActiveSource(gradient2)

# create a new 'Resample With Dataset'
resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=gradient2,
    DestinationMesh=clip2)
resampleWithDataset1.CellLocator = 'Static Cell Locator'

UpdatePipeline(time=0.0, proxy=resampleWithDataset1)

# set active source
SetActiveSource(resampleWithDataset1)

# set active source
SetActiveSource(resampleWithDataset1)

# rename source object
RenameSource('extract1', resampleWithDataset1)

# save data
SaveData('/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/surface.csv', proxy=resampleWithDataset1, WriteTimeSteps=1,
    WriteTimeStepsSeparately=1,
    Filenamesuffix='_%.4d',
    PointDataArrays=['Divergence', 'Gradient', 'Q Criterion', 'T', 'Vorticity', 'continent', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'density', 'depth', 'nonadiabatic_temperature', 'p', 'plastic_yielding', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'specific_heat', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity', 'vtkValidPointMask'],
    FieldDataArrays=['TIME'],
    AddTime=1)


# set active source
SetActiveSource(contour3)

# set active source
SetActiveSource(calculator1)

# save data
SaveData('/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/lithosphere.csv', proxy=calculator1, WriteTimeSteps=1,
    WriteTimeStepsSeparately=1,
    Filenamesuffix='_%.4d',
    PointDataArrays=['Divergence', 'Gradient', 'Normals', 'Q Criterion', 'T', 'Vorticity', 'continent', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'density', 'depth', 'nonadiabatic_temperature', 'oceanic_age', 'p', 'plastic_yielding', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'specific_heat', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity'],
    FieldDataArrays=['TIME'],
    AddTime=1)

# set active source
SetActiveSource(contour3)

# save data
SaveData('/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/depths.csv', proxy=contour3, WriteTimeSteps=1,
    WriteTimeStepsSeparately=1,
    Filenamesuffix='_%.4d',
    PointDataArrays=['Divergence', 'Gradient', 'Normals', 'Q Criterion', 'T', 'Vorticity', 'continent', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'density', 'depth', 'nonadiabatic_temperature', 'p', 'plastic_yielding', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'specific_heat', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity'],
    FieldDataArrays=['TIME'],
    AddTime=1)



# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()