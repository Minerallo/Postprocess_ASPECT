#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *

path_models = ['/scratch/usr/bbkponsm/model/globalscale/sphere3d/R01a_Rodinia_1GPa_Mantle_C40MPa_LR/',
    '/scratch/usr/bbkponsm/model/globalscale/sphere3d/S01a_no_continents_C40MPA_TP1700_LR/',
    '/scratch/usr/bbkponsm/model/globalscale/sphere3d/S01h_no_continents_C10MPA_f005_TP1700_LR/']

# filename=path_models[0] + 'solution_surface.pvd'
# filename2=path_models[0] + 'solution.pvd'

# folder = path_models[0].split('/')[-2]

for path in path_models:
    filename=path + 'solution_surface.pvd'
    filename2=path + 'solution.pvd'
    folder = path.split('/')[-2]
    output_base_folder = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/'

    # '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/'
    # '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/'
    # '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/'

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

    # create a new 'PVD Reader'
    solution_surfacepvd = PVDReader(registrationName='solution_surface.pvd', FileName=filename)
    solution_surfacepvd.PointArrays = []

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # create a new 'PVD Reader'
    solutionpvd = PVDReader(registrationName='solution.pvd', FileName=filename2)
    solutionpvd.PointArrays = ['continent','strain_rate', 'v_r', 'v_phi', 'v_theta']

    # find source
    solutionpvd_1 = FindSource('solution.pvd')

    # set active source
    SetActiveSource(solutionpvd_1)

    # set active source
    SetActiveSource(solutionpvd)


    # tsteps = solution_surfacepvd_1.TimestepValues
    # print(tsteps)

    output_directory_full =output_model_folder+'/' 
    # n_files_surface = sum(1 for file in os.listdir(output_directory_full) if 'surface' in file)
    # n_files_lithosphere = sum(1 for file in os.listdir(output_directory_full) if 'lithosphere' in file)
    # n_files_depths = sum(1 for file in os.listdir(output_directory_full) if 'depths' in file)

    # # check if we reached the maximum steps for each files 
    # if n_files_surface<len(tsteps):
    #     n_file_choosen = n_files_surface
    #     index_file = 1
    #     last_time = 0
    # elif n_files_surface==len(tsteps) & n_files_lithosphere < len(tsteps):
    #     n_files_choosen = n_files_lithosphere
    #     index_file = 2
    #     last_time = 0
    # else :
    #     n_files_choosen = n_files_depths
    #     index_file = 3
    #     last_time = tsteps[n_files_depths-1]

    # sum(1 for file in os.listdir('/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/') if 'surface' in file)

    UpdatePipeline(time=0, proxy=solutionpvd)

    UpdatePipeline(time=0, proxy=solution_surfacepvd)

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

    UpdatePipeline(time=0, proxy=clip2)

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

    # create a new 'Resample With Dataset'
    resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=solutionpvd,
        DestinationMesh=clip2)
    resampleWithDataset1.CellLocator = 'Static Cell Locator'

    UpdatePipeline(time=0, proxy=resampleWithDataset1)

    # set active source
    SetActiveSource(resampleWithDataset1)

    # rename source object
    RenameSource('extract1', resampleWithDataset1)

    #  & os.path.isfile(output_model_folder + '/surface.csv')==False
    # if index_file == 1 :
        # save data
    SaveData(output_model_folder + '/surface.csv', proxy=resampleWithDataset1, WriteTimeSteps=1,
        WriteTimeStepsSeparately=1,
        Filenamesuffix='_%.4d',
        PointDataArrays=['continent', 'strain_rate', 'v_phi', 'v_r', 'v_theta','vtkValidPointMask'],
        FieldDataArrays=['TIME'],
        AddTime=1)
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    ResetSession()