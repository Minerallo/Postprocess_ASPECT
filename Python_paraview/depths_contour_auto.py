# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
import os

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

path_models = ['/scratch/usr/bbkponsm/model/globalscale/sphere3d/V06c_R01f_Rodinia_2GPa_llsvps_ps_1x50My_init_2Myint_rhoc3160/']


for path in path_models:
    filename2=path + 'solution.pvd'

    folder = path.split('/')[-2]
    output_base_folder = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/'
    check_path_from_sshfs = output_base_folder
    
    output_model_folder = output_base_folder + folder
    check_model_files = check_path_from_sshfs + folder

    try:
        os.makedirs(check_model_files)
    except FileExistsError:
        print("model repository for postprocess already exists")
        # directory already exists
        pass

    # create a new 'PVD Reader'
    solutionpvd = PVDReader(registrationName='solution.pvd', FileName=filename2)
    solutionpvd.PointArrays = ['velocity', 'stress', 'shear_stress', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'p', 'T', 'continent', 'asthenosphere', 'strain_rate', 'nonadiabatic_temperature', 'depth', 'vertical_heat_flux', 'heat_flux_map', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity', 'thermal_conductivity', 'thermal_diffusivity', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'plastic_yielding', 'v_r', 'v_phi', 'v_theta', 'principal_stress_1', 'principal_stress_2', 'principal_stress_3']
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    
    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()
    
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    solutionpvd_1 = FindSource('solution.pvd')

    
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
    contour1.Isosurfaces = [100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2880000]
#     [100000.0, 298571.4285714286, 497142.85714285716, 695714.2857142857, 894285.7142857143, 1092857.142857143, 1291428.5714285714, 1490000.0, 1688571.4285714286, 1887142.8571428573, 2085714.285714286, 2284285.7142857146, 2482857.1428571427, 2681428.5714285714, 2880000.0]
    
    UpdatePipeline(time=0.0, proxy=contour1)
    

    tsteps = solutionpvd_1.TimestepValues
    print(tsteps)
    
    check_model_files_full = check_model_files+'/' 

    n_files = sum(1 for file in os.listdir(check_model_files_full) if 'depths_contours' in file)

    for i in range(n_files, len(tsteps)):
        time_writing = tsteps[i]
        message2 = f"Current time {time_writing/1e6} My."
        print(message2)
        # Update all the values
        UpdatePipeline(time=tsteps[i], proxy=contour1)    
    
        # save data
        message = f"Writing depth contours files at time  {time_writing}."
        print(message)
        SetActiveSource(contour1)
        SaveData(output_model_folder +'/depths_contours_%.4d.csv' % i, proxy=contour1,
            ChooseArraysToWrite=1,
            PointDataArrays=['T', 'density', 'depth', 'viscosity'],
            AddTime=1)


# # trace generated using paraview version 5.11.2
# #import paraview
# #paraview.compatibility.major = 5
# #paraview.compatibility.minor = 11
# 
# #### import the simple module from the paraview
# from paraview.simple import *
# import os
# 
# #### disable automatic camera reset on 'Show'
# paraview.simple._DisableFirstRenderCameraReset()
# 
# path_models = ['/scratch/usr/bbkponsm/model/globalscale/sphere3d/R01e_Rodinia_2GPa_Mantle_C20MPa_f003_LR/']
# 
# for path in path_models:
#     filename2=path + 'solution.pvd'
# 
#     folder = path.split('/')[-2]
#     output_base_folder = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/'
#     check_path_from_sshfs = output_base_folder
#     
#     output_model_folder = output_base_folder + folder
#     check_model_files = check_path_from_sshfs + folder
# 
#     try:
#         os.makedirs(check_model_files)
#     except FileExistsError:
#         print("model repository for postprocess already exists")
#         # directory already exists
#         pass
# 
#     # create a new 'PVD Reader'
#     solutionpvd = PVDReader(registrationName='solution.pvd', FileName=filename2)
#     solutionpvd.PointArrays = ['velocity', 'stress', 'shear_stress', 'principal_stress_direction_1', 'principal_stress_direction_2', 'principal_stress_direction_3', 'p', 'T', 'continent', 'asthenosphere', 'strain_rate', 'nonadiabatic_temperature', 'depth', 'vertical_heat_flux', 'heat_flux_map', 'density', 'thermal_expansivity', 'specific_heat', 'viscosity', 'thermal_conductivity', 'thermal_diffusivity', 'current_cohesions', 'current_friction_angles', 'current_yield_stresses', 'plastic_yielding', 'v_r', 'v_phi', 'v_theta', 'principal_stress_1', 'principal_stress_2', 'principal_stress_3']
#     
#     # get animation scene
#     animationScene1 = GetAnimationScene()
#     
#     # get the time-keeper
#     timeKeeper1 = GetTimeKeeper()
#     
#     # update animation scene based on data timesteps
#     animationScene1.UpdateAnimationUsingDataTimeSteps()
#     
#     # Properties modified on solutionpvd
#     solutionpvd.PointArrays = ['T', 'density', 'depth', 'viscosity']
#     
#     UpdatePipeline(time=0.0, proxy=solutionpvd)
#     
#     # Properties modified on animationScene1
#     animationScene1.Stride = 20
#     
#     # create a new 'Contour'
#     contour1 = Contour(registrationName='Contour1', Input=solutionpvd)
#     contour1.ContourBy = ['POINTS', 'T']
#     contour1.Isosurfaces = [2133.2440185546875]
#     contour1.PointMergeMethod = 'Uniform Binning'
#     
#     # Properties modified on contour1
#     contour1.ContourBy = ['POINTS', 'depth']
#     contour1.Isosurfaces = [100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 1100000, 1200000, 1300000, 1400000, 1500000, 1600000, 1700000, 1800000, 1900000, 2000000, 2100000, 2200000, 2300000, 2400000, 2500000, 2600000, 2700000, 2800000, 2880000]
# #     [100000.0, 298571.4285714286, 497142.85714285716, 695714.2857142857, 894285.7142857143, 1092857.142857143, 1291428.5714285714, 1490000.0, 1688571.4285714286, 1887142.8571428573, 2085714.285714286, 2284285.7142857146, 2482857.1428571427, 2681428.5714285714, 2880000.0]
#     
#     UpdatePipeline(time=0.0, proxy=contour1)
#     
#     # save data
#     SaveData(output_model_folder +'/depths_contours.csv', proxy=contour1, WriteTimeSteps=1,
#         WriteTimeStepsSeparately=1,
#         ChooseArraysToWrite=1,
#         PointDataArrays=['T', 'density', 'depth', 'viscosity'],
#         AddTime=1)