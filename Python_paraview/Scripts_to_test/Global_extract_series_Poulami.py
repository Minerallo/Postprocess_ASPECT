# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import os

# path_models = [
#     '/scratch/usr/bbkponsm/model/globalscale/sphere3d/P01a_Pangea_1GPa_Mantle_C40MPa_LR/',
#     '/scratch/usr/bbkponsm/model/globalscale/sphere3d/P01b_Pangea_1GPa_Mantle_C40MPa_HR/']
# path_models = ['/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/']

path_models = ['/scratch/projects/bbk00014/model_global_Poulami/12_test_250Ma_3D_test_no_shear_heat/']

for path in path_models:
    filename2=path + 'solution.pvd'

    folder = path.split('/')[-2]
    output_base_folder = '/scratch/projects/bbk00014/model_global_Poulami/post_process/'
    

    output_model_folder = output_base_folder + folder

    # try:
    #     os.makedirs(output_model_folder)
    # except FileExistsError:
    #     print("file already exists")
    #     # directory already exists
    #     pass

    # for i in range(0,len(path_models)+1):

    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # create a new 'PVD Reader'
    solutionpvd = PVDReader(registrationName='solution.pvd', FileName=filename2)
    solutionpvd.PointArrays = [ 'T', 'density', 'depth', 'nonadiabatic_temperature', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'viscosity','llsvps','thermal_conductivity','thermal_expansivity','velocity','p']

    # find source
    solutionpvd_1 = FindSource('solution.pvd')

    # set active source
    SetActiveSource(solutionpvd_1)

    UpdatePipeline(time=0, proxy=solutionpvd)
    
    
    # create a new 'Extract Time Steps'
    extractTimeSteps1 = ExtractTimeSteps(registrationName='ExtractTimeSteps1', Input=solutionpvd)
    extractTimeSteps1.TimeStepIndices = [0]
    extractTimeSteps1.TimeStepRange = [0, 250]

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # Properties modified on extractTimeSteps1
    extractTimeSteps1.SelectionMode = 'Select Time Range'
    extractTimeSteps1.TimeStepInterval = 5

    UpdatePipeline(time=0.0, proxy=extractTimeSteps1)
    
       # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # set active source
    SetActiveSource(extractTimeSteps1)

    
    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=extractTimeSteps1)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=slice1.SliceType)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice1.SliceType)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=slice1.SliceType)

    # Properties modified on slice1
    slice1.SliceType = 'Sphere'

    # Properties modified on slice1.SliceType
    slice1.SliceType.Radius = 6370000.0

    UpdatePipeline(time=0.0, proxy=slice1)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=slice1.SliceType)

    # toggle interactive widget visibility (only when running from the GUI)
    ShowInteractiveWidgets(proxy=slice1.SliceType)
    
    # set active source
    SetActiveSource(extractTimeSteps1)

    # create a new 'Contour'
    contour3 = Contour(registrationName='Contour3', Input=extractTimeSteps1)
    contour3.ContourBy = ['POINTS', 'Divergence']
    contour3.Isosurfaces = [3.7377972148533445e-06]
    contour3.PointMergeMethod = 'Uniform Binning'

    # set active source
    SetActiveSource(contour3)

    # Properties modified on contour3
    contour3.ContourBy = ['POINTS', 'depth']
    contour3.Isosurfaces = [440000.0, 1100000.0, 2600000.0, 2800000.0]

    UpdatePipeline(time=0, proxy=contour3)

    # rename source object
    RenameSource('extract3', contour3)

  # rename source object
    RenameSource('extract1', slice1)
  

    # set active source
    SetActiveSource(slice1)


    #  & os.path.isfile(output_model_folder + '/surface.csv')==False

        # save data
    SaveData(output_model_folder + '/surface.csv', proxy=slice1, WriteTimeSteps=1,
        WriteTimeStepsSeparately=1,
        Filenamesuffix='_%.4d',
        PointDataArrays=['T', 'nonadiabatic_temperature', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'viscosity','llsvps','velocity','p'],
        FieldDataArrays=['TIME'],
        AddTime=1)

    # set active source
    SetActiveSource(contour3)

    # save data
    SaveData(output_model_folder +'/depths.csv', proxy=contour3, WriteTimeSteps=1,
        WriteTimeStepsSeparately=1,
        Filenamesuffix='_%.4d',
        PointDataArrays=[ 'T', 'depth', 'nonadiabatic_temperature', 'v_phi', 'v_r', 'v_theta', 'viscosity','llsvps','thermal_conductivity','thermal_expansivity','velocity','p'],
        FieldDataArrays=['TIME'],
        AddTime=1)




    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    ResetSession()
