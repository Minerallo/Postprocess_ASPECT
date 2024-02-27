# trace generated using paraview version 5.11.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import os

# path_models = /scratch/usr/bbkponsm/model/globalscale/sphere3d/R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003_test/
# path_models = /scratch/usr/bbkponsm/model/globalscale/sphere3d/R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003_test/

path_models = /scratch/usr/bbkponsm/model/globalscale/sphere3d/R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003_test/

#For the purpose of automatizing the postprocess we first loop over the models path to create the repositories so we can track them locally 
# for path in path_models:
#     folder = path.split('/')[-2]
#     check_path_from_sshfs = '/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/sphere3d/'
#     check_model_files = check_path_from_sshfs + folder

#     try:
#         os.makedirs(check_model_files)
#     except FileExistsError:
#         print("model repository for postprocess already exists")
#         # directory already exists
#         pass

#Here is the original loop 
for path in path_models:
    filename=path + 'solution_surface.pvd'
    filename2=path + 'solution.pvd'

    folder = path.split('/')[-2]
    output_base_folder = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/'
    check_path_from_sshfs = '/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/sphere3d/'
    # '/Users/ponsm/Desktop/Model_sphere/Models_HLRN/'

    output_model_folder = output_base_folder + folder
    check_model_files = check_path_from_sshfs + folder
    # path + 'Postprocess'
    # output_base_folder + folder

    try:
        os.makedirs(check_model_files)
    except FileExistsError:
        print("model repository for postprocess already exists")
        # directory already exists
        pass

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
    solutionpvd.PointArrays = ['T', 'continent', 'current_friction_angles', 'current_yield_stresses', 'density', 'depth', 'nonadiabatic_temperature', 'p', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity']

    # find source
    solutionpvd_1 = FindSource('solution.pvd')

    # set active source
    SetActiveSource(solutionpvd_1)

    # set active source
    SetActiveSource(solutionpvd)

  
  # find source
    solution_surfacepvd_1 = FindSource('solution_surface.pvd')

    UpdatePipeline(time=0, proxy=solutionpvd)

    # Properties modified on solution_surfacepvd

    UpdatePipeline(time=0, proxy=solution_surfacepvd)

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

    UpdatePipeline(time=0, proxy=clip2)

    # toggle interactive widget visibility (only when running from the GUI)
    HideInteractiveWidgets(proxy=clip2.ClipType)

    # set active source
    SetActiveSource(clip2)

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

    UpdatePipeline(time=0, proxy=gradient2)


    # set active source
    SetActiveSource(gradient2)


    # create a new 'Contour'
    contour2 = Contour(registrationName='Contour2', Input=gradient2)
    contour2.ContourBy = ['POINTS', 'Divergence']
    contour2.Isosurfaces = [3.7377972148533445e-06]
    contour2.PointMergeMethod = 'Uniform Binning'

    # Properties modified on contour2
    contour2.ContourBy = ['POINTS', 'T']

    UpdatePipeline(time=0, proxy=contour2)

    # Properties modified on contour2
    contour2.Isosurfaces = [1573.0]


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

    UpdatePipeline(time=0, proxy=calculator1)

    # rename source object
    RenameSource('extract2', calculator1)

    # find source
    extract3 = FindSource('extract3')

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

    # Properties modified on contour3
    contour3.ContourBy = ['POINTS', 'depth']
    contour3.Isosurfaces = [440000.0, 2600000.0]

    UpdatePipeline(time=0, proxy=contour3)

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

    UpdatePipeline(time=0, proxy=resampleWithDataset1)

    # set active source
    SetActiveSource(resampleWithDataset1)

    # rename source object
    RenameSource('extract1', resampleWithDataset1)

    tsteps = solution_surfacepvd_1.TimestepValues
    print(tsteps)
    
    check_model_files_full = check_model_files+'/' 
    
    n_files_surface = sum(1 for file in os.listdir(check_model_files_full) if 'surface' in file)
    n_files_lithosphere = sum(1 for file in os.listdir(check_model_files_full) if 'lithosphere' in file)
    n_files_depths = sum(1 for file in os.listdir(check_model_files_full) if 'depths' in file)

    min_files = min(n_files_surface, n_files_lithosphere, n_files_depths)

# Working with local path in sshfs 
        # sum(1 for file in os.listdir('/Users/ponsm/Desktop/Model_sphere/Models_HLRN/P01a_Pangea_1GPa_Mantle_C40MPa_LR/') if 'surface' in file)
# < len(tsteps)
    for i in range(min_files, len(tsteps)):
        time_writing = tsteps[i]
        message2 = f"Current time {time_writing/1e6} My."
        print(message2)

        # Update all the values
        UpdatePipeline(time=tsteps[i], proxy=calculator1)
        UpdatePipeline(time=tsteps[i], proxy=contour3)
        UpdatePipeline(time=tsteps[i], proxy=resampleWithDataset1)

        if i>n_files_surface-1:
            message = f"Writing surface files at time {time_writing}."
            print(message)
            # set active source
            SetActiveSource(resampleWithDataset1)
            SaveData(output_model_folder +'/surface_%.4d.csv' % i, proxy=resampleWithDataset1, 
                PointDataArrays=['Divergence', 'T', 'Vorticity', 'continent', 'current_friction_angles', 'current_yield_stresses', 'density', 'nonadiabatic_temperature', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity'],
                FieldDataArrays=['TIME'],
                AddTime=1)

        if i>n_files_lithosphere-1:
            message = f"Writing lithosphere files at time  {time_writing}."
            print(message)
            # set active source
            SetActiveSource(calculator1)
            # save data
            SaveData(output_model_folder +'/lithosphere_%.4d.csv' % i, proxy=calculator1,
                Filenamesuffix='_%.4d',
                PointDataArrays=['Divergence', 'T', 'Vorticity', 'continent', 'density', 'depth', 'nonadiabatic_temperature', 'oceanic_age', 'p', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity'],
                FieldDataArrays=['TIME'],
                AddTime=1)

        if i>n_files_depths-1:
            message = f"Writing depths files at time {time_writing}."
            print(message)
            # set active source
            SetActiveSource(contour3)
            # save data
            SaveData(output_model_folder +'/depths_%.4d.csv' % i, proxy=contour3,
                PointDataArrays=[ 'T', 'Vorticity', 'density', 'depth', 'nonadiabatic_temperature', 'p', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'shear_stress', 'strain_rate', 'v_phi', 'v_r', 'v_theta', 'velocity', 'viscosity'],
                FieldDataArrays=['TIME'],
                AddTime=1)    


    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    ResetSession()