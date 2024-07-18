#!/bin/bash


# # This script automatize the the extraction of the data and the postprocessing of the model to the making of the figures  

echo "Building sshfs connection of the remote server."

umount -f ~/Desktop/modelblogin

sshfs -o "IdentityFile=/Users/ponsm/.ssh/id_rsa_hlrn" bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm /Users/ponsm/Desktop/modelblogin -o defer_permissions

echo "Building and sending files to the server"
/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run('/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/Python_paraview/build_autopostprocess_files.m'); catch;disp('Error: The MATLAB auto-postprocess files could not be properly built.'); end; quit"

# rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/auto_generated_postprocess/Model_extract_global_3D_auto.py bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/.
# rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/Python_paraview/jobparaview_global3D.sh bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/.

# echo "Running Pvbatch on the remote server"
# ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de 'cd /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/; sbatch -p large96 -t 12:00:00 -N 1 --tasks-per-node 96 -o %j.output -e %j.ouput --mem=747000mb --account bbkponsm --job-name PPv /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/jobparaview_global3D.sh'

# # Wait for the job to finish
# while ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de 'squeue -n PPv' | grep -q PPv; do
#     echo "Job Pv batch is running on the remote server. Waiting..."
#     sleep 20  # Adjust the sleep interval
# done

# echo "Paraview Job has finished. Proceeding to the copy of the data to local."

# chmod +x /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/auto_generated_postprocess/Models_to_copy_from_remote.sh

# /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/auto_generated_postprocess/Models_to_copy_from_remote.sh

# echo "Running Matlab Postprocess and making figures"
# # /Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "run('/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/auto_generated_postprocess/Postprocess_auto_generated'); catch; end; quit"
# /Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run('/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/auto_generated_postprocess/Postprocess_auto_generated'); catch; disp('Error: The MATLAB Postprocess encountered an issue.'); end; quit"
# echo "Postprocess succeeded and figures extracted"



# # working test

# /Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "disp(['Current folder: ' pwd])"
# /Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run('/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/Python_paraview/bash_matlab_test.m'); catch; end; quit"
# matlab -nosplash -noFigureWindows -r "try; run('/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/Python_paraview/bash_matlab_test.m'); catch; end; quit"

# sshfs -o "IdentityFile=/Users/ponsm/.ssh/id_rsa_hlrn" bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm /Users/ponsm/Desktop/modelblogin -o defer_permissions

# rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/Python_paraview/test_simple_sphere.py bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/

# rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/Python_paraview/jobparaview_global3D.sh bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/

# ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de 'cd /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/; sbatch -p large96 -t 12:00:00 -N 1 --tasks-per-node 96 -o %j.output -e %j.ouput --mem=747000mb --account bbkponsm --job-name PPv /scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/jobparaview_global3D.sh'

# # Wait for the job to finish
# while ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de 'squeue -n PPv' | grep -q PPv; do
#     echo "Job is still running. Waiting..."
#     sleep 10  # Adjust the sleep interval as needed
# done



# rsync -varz -e 'ssh -i ~/.ssh/id_rsa_hlrn' --include="*statistics*" --exclude="*" bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm/model/globalscale/sphere3d/S01a_no_continents_C40MPA_TP1700_LR/ /Users/ponsm/Desktop/Model_sphere/Models_HLRN/S01a_no_continents_C40MPA_TP1700_LR/ ;


