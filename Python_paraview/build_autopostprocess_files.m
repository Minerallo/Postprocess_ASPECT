clc; clear all;
% This scripts get the model names and paths given by the user and automatically
% apply them to the Pvpython script and the mantlab postsprocess
path_to_postprocess_ASPECT = '/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/';
paraview_python_script = [path_to_postprocess_ASPECT 'Python_paraview/Model_extract_global_3D_auto.py'];
models_cluster_directory = '/scratch/usr/bbkponsm/model/globalscale/sphere3d/';
extracted_models_pvdata_directory = '/Volumes/Jerry/global_models_3d_extract/';
remote_Postprocess_scripts = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/';
remote_Postprocess_directory = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/';
local_access_to_postprocess_directory = '/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/sphere3d/';
local_access_to_postprocess_scripts_directory='/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/postprocess_scripts/';

%path of the Pvpython to modify
modified_file_path = './Model_extract_global_3D_auto.py';

cd(path_to_postprocess_ASPECT);

%% List of specific model names to give
% model_names = {'R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003', ...
%     'R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003'};
    
model_names = {'R01e_Rodinia_2GPa_Mantle_C20MPa_f003_LR'};
%     'RSRW01a_Rodinia_2GPa_Mantle_C10MPa_SW_1em16_f064_to_f003_LR',...
%     'P01a_Pangea_1GPa_Mantle_C40MPa_LR',  ...
%     'P01c_Pangea_1GPa_Mantle_C60MPa_LR',  ...
%     'P01e_Pangea_2GPa_Mantle_C20MPa_f003_LR',  ...
%     'P01f_Pangea_2GPa_Mantle_C10MPa_f005_LR',  ...
%     'P01g_Pangea_05GPa_Mantle_C40MPa_LR',  ...
%     'P01h_Pangea_2GPa_Mantle_C40MPa_LR',  ...
%     'P01i_Pangea_2GPa_Mantle_C60MPa_LR',  ...
%     'R01a_Rodinia_1GPa_Mantle_C40MPa_LR',  ...
%     'R01c_Rodinia_1GPa_Mantle_C60MPa_LR',  ...
%     'R01e_Rodinia_2GPa_Mantle_C20MPa_f003_LR', ...
%     'R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR',  ...
%     'R01g_Rodinia_05GPa_Mantle_C40MPa_LR',  ...
%     'R01h_Rodinia_2GPa_Mantle_C40MPa_LR',  ...
%     'R01i_Rodinia_2GPa_Mantle_C60MPa_LR',  ...
%     'T01a_Pangea_1GPa_Mantle_C60MPa_LR',  ...
%     'T01b_Pangea_1GPa_Mantle_C60MPa_LR',  ...
%     'T01e_Rodinia_1GPa_Mantle_C60MPa_LR'
% };


%% create a repository where all auto_generated files will go
% Create or check existence of the directory 'auto_generated_postprocess'


% Give a name or it will give for name auto_ followed by the list of models
% using their first word before "_"
output_directory_name = '';


    % Use regular expression to extract the first model names before underscore
    pattern = '([^_]+)_'; % Match one or more characters that are not underscores, followed by an underscore
    matches = regexp(model_names, pattern, 'tokens', 'once');
    % Extract the first match from each cell
    first_model_names = cellfun(@(x) x{1}, matches, 'UniformOutput', false);
    % Create a list
    model_list = strjoin(first_model_names, '_');
if strcmp(output_directory_name, '')
    output_directory = ['auto_' model_list];
else
    output_directory =  output_directory_name;
end
% disp(model_list);
if ~exist(output_directory, 'dir')
    mkdir(output_directory);
end

%% Read the content of the Python script and modify the models list
fid = fopen(paraview_python_script, 'r');
python_script_content = fread(fid, '*char')';
fclose(fid);

% Build the full paths to models by concatenating the base directory with model names
models_list_to_postprocess = cellfun(@(name) fullfile([models_cluster_directory name '/']), model_names, 'UniformOutput', false);

% Convert MATLAB cell array to Python list format
new_path_models_python = cellfun(@(x) ['''' x ''''], models_list_to_postprocess, 'UniformOutput', false);

% Replace the existing path_models in the Python script content
python_script_content = regexprep(python_script_content, 'path_models\s*=\s*\[.*?\]', ['path_models = [' strjoin(new_path_models_python, ',\n        ') ']\n']);

output_file = fullfile(output_directory, modified_file_path);

% Save the modified script
fid = fopen(output_file, 'w');
fwrite(fid, python_script_content);
fclose(fid);

disp(['Modified Pvpython script saved to: ' output_file]);

%% Copy and modify the model list for Postprocess main
% Specify the path to the original Postprocess_v1_2.m script
postprocess_script_path = [path_to_postprocess_ASPECT 'Python_paraview/Postprocess_v1_2_auto.m'];

% Specify the path to save the modified copy
postprocess_modified_script_path = [path_to_postprocess_ASPECT 'Postprocess_v1_2_auto_generated.m'];

% Define the new path_models variable for the copied script

% new_path_models = cellfun(@(name) fullfile(extracted_models_pvdata_directory, name), model_names, 'UniformOutput', false);

% Build the full paths to models by concatenating the base directory with model names
new_path_models = cellfun(@(name) [extracted_models_pvdata_directory name '/'], model_names, 'UniformOutput', false);

% Convert the models_list_to_postprocess to the desired format
formatted_path_models = cellfun(@(path) ['''' path ''''], new_path_models, 'UniformOutput', false);

path_models_variable = ['path_models = {[' strjoin(formatted_path_models, '],\n    [') ']};'];

fid_original = fopen(postprocess_script_path, 'r');
original_script_content = fread(fid_original, '*char')';
fclose(fid_original);

existing_addpath_pattern = 'addpath\(genpath\(pwd\)\);';
modified_script_content = regexprep(original_script_content, existing_addpath_pattern, ...
    ['current_file_path = pwd;', newline, 'parent_directory = fileparts(current_file_path);', newline, 'addpath(genpath(parent_directory));'], 'dotall');

% Find the existing path_models and replace it with the new one
existing_path_models_pattern = 'path_models\s*=\s*\{.*?\};';
modified_script_content = regexprep(modified_script_content, existing_path_models_pattern, path_models_variable, 'dotall');

output_file_postprocess = fullfile(output_directory, 'Postprocess_auto_generated.m');

% Save the modified copy
fid_modified = fopen(output_file_postprocess, 'w');
fwrite(fid_modified, modified_script_content);
fclose(fid_modified);

disp(['Modified matlab script saved to: ' output_file_postprocess]);



%% generate a bashscript containing the models path of the Pv data that will need to be copy to local

% Define the common variables
remote_base_directory = remote_Postprocess_directory;
local_base_directory = extracted_models_pvdata_directory;
ssh_key = '~/.ssh/id_rsa_hlrn';

% List of model names
% Combine model names into a single string separated by spaces
model_names_str = strjoin(model_names, ' ');

output_file_bash_copy = fullfile(output_directory, 'Models_to_copy_from_remote.sh');


% Save the Bash script to a file using fprintf with explicit newline characters
fid = fopen(output_file_bash_copy, 'w');
fprintf(fid, '#!/bin/bash\n\n');
% fprintf(fid, '# Set common variables\n');
% fprintf(fid, 'remote_base_directory=''%s''\n', remote_base_directory);
% fprintf(fid, 'local_base_directory=''%s''\n', local_base_directory);
% fprintf(fid, 'remote_model_directory=''%s''\n', models_cluster_directory);
% fprintf(fid, 'ssh_key=''%s''\n\n', ssh_key);
% fprintf(fid, '# Loop through the list of model names\n');
% fprintf(fid, 'for model in %s; do\n', model_names_str);
% fprintf(fid, '    # Construct the remote and local paths\n');
% fprintf(fid, '    remote_path="${remote_base_directory}${model}"\n');
% fprintf(fid, '    remote_model_path="${remote_model_directory}${model}"\n');
% fprintf(fid, '    local_path="${local_base_directory}"\n\n');
% fprintf(fid, '    # Perform rsync\n');
% fprintf(fid, '    rsync -varz -e "ssh -i $ssh_key" "bbkponsm@blogin.hlrn.de:$remote_path" "$local_path"\n');
% fprintf(fid, '    rsync -varz -e "ssh -i $ssh_key" --include="*statistics*" --exclude="*" "bbkponsm@blogin.hlrn.de:$remote_model_path" "$local_path"\n');
% fprintf(fid, '    rsync -varz -e ''ssh -i %s'' --include="*statistics*" --exclude="*" bbkponsm@blogin.hlrn.de:%s${model} %s${model}/.\n',ssh_key, models_cluster_directory, local_base_directory);

for i = 1:numel(model_names)
    fprintf(fid, 'rsync -varz -e ''ssh -i %s'' bbkponsm@blogin.hlrn.de:%s%s/ %s%s/.\n',ssh_key,remote_base_directory, model_names{i}, local_base_directory, model_names{i});
    fprintf(fid, 'rsync -varz -e ''ssh -i %s'' --include="*statistics*" --exclude="*" bbkponsm@blogin.hlrn.de:%s%s/ %s%s/.\n', ssh_key, models_cluster_directory, model_names{i}, local_base_directory, model_names{i});
    fprintf(fid, 'rsync -varz -e ''ssh -i %s'' --include="depth_average.txt" --exclude="*" bbkponsm@blogin.hlrn.de:%s%s/ %s%s/.\n', ssh_key, models_cluster_directory, model_names{i}, local_base_directory, model_names{i});
end

fclose(fid);

disp(['Bash_copy script saved to: ' output_file_bash_copy]);

%% make a jobscript for Pv batch in the autogenerated folder.sh

output_pvbatch_jobscript = fullfile(output_directory, 'jobparaview_global3D.sh');

% Save the Bash script to a file using fprintf with explicit newline characters
fid_slurm = fopen(output_pvbatch_jobscript, 'w');
fprintf(fid_slurm, '#!/bin/bash\n\n');
fprintf(fid_slurm, '#SBATCH -A bbk00014\n');
fprintf(fid_slurm, 'module load gcc/9.2.0 openmpi/gcc.9 anaconda3 llvm/9.0.0 paraview\n\n');
fprintf(fid_slurm, '# Ensure the cpus-per-task option is propagated to srun commands\n');
fprintf(fid_slurm, 'export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK\n\n');
fprintf(fid_slurm, 'mpirun --map-by socket:pe=$OMP_NUM_THREADS pvbatch %s%s/Model_extract_global_3D_auto.py\n',remote_Postprocess_scripts,output_directory);
fclose(fid_slurm);

disp(['Your Slurm script built and saved to: ' output_pvbatch_jobscript]);


%% make a copy of this file to the autogenerated folder
file_input_sh = [path_to_postprocess_ASPECT 'Python_paraview/build_autopostprocess_files.m'];
copyfile(file_input_sh, output_directory);
%% Check and create local repository corresponding to the model names
for i = 1:length(model_names)
    repository_name = model_names{i};
    repository_path = fullfile(extracted_models_pvdata_directory, repository_name);

    % Check if the repository exists
    if exist(repository_path, 'dir') == 7
        disp(['Repository "', repository_name, '" already exists.']);
    else
        % Create the repository if it doesn't exist
        mkdir(repository_path);
        disp(['Repository "', repository_name, '" created successfully.']);
    end
end

%% Check and create remote repository corresponding to the model names
for i = 1:length(model_names)
    repository_name = model_names{i};
    repository_path = fullfile(local_access_to_postprocess_directory, repository_name);

    % Check if the repository exists
    if exist(repository_path, 'dir') == 7
        disp(['Repository "', repository_name, '" already exists.']);
    else
        % Create the repository if it doesn't exist
        mkdir(repository_path);
        disp(['Repository "', repository_name, '" created successfully.']);
    end
end

%% Check and create remote repository to store the postprocess files corresponding
repository_path = fullfile(local_access_to_postprocess_scripts_directory, output_directory);

% Check if the repository exists
if exist(repository_path, 'dir') == 7
    disp(['Remote scripts repository "', output_directory, '" already exists.']);
else
    % Create the repository if it doesn't exist
    mkdir(repository_path);
    disp(['Remote scripts repository "', output_directory, '" created successfully.']);
end

%% make the sh file that will run all the scripts (could be more flexible in the future)

output_file_sh = fullfile(output_directory, 'Postprocess_auto.sh');
if numel(model_list) > 3
    last_letters = model_list(end-3:end);
else
    last_letters = model_list;
end

jobname =['Pv',last_letters];

% Save the Bash script to a file using fprintf with explicit newline characters
fid = fopen(output_file_sh, 'w');
fprintf(fid, '#!/bin/bash\n\n');
fprintf(fid, '\n');
fprintf(fid,'# This script automatize the the extraction of the data and the postprocessing of the model to the making of the figures\n');
fprintf(fid, '\n');
fprintf(fid,'echo "Building sshfs connection of the remote server."\n');
fprintf(fid, '\n');
fprintf(fid, 'umount -f ~/Desktop/modelblogin\n');
fprintf(fid, '\n');
fprintf(fid, 'sshfs -o "IdentityFile=/Users/ponsm/.ssh/id_rsa_hlrn" bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm /Users/ponsm/Desktop/modelblogin -o defer_permissions\n');
fprintf(fid, '\n');
fprintf(fid, 'echo "Building and sending files to the server"\n');
fprintf(fid, '##/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run(''/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/build_autopostprocess_files.m''); catch;disp(''Error: The MATLAB auto-postprocess files could not be properly built.''); end; quit"\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, 'rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Model_extract_global_3D_auto.py bbkponsm@blogin.hlrn.de:%s%s/.\n',output_directory,remote_Postprocess_scripts,output_directory);
fprintf(fid, 'rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/jobparaview_global3D.sh bbkponsm@blogin.hlrn.de:%s%s/.\n',output_directory,remote_Postprocess_scripts,output_directory);
fprintf(fid, '\n');
fprintf(fid, 'echo "Running Pvbatch on the remote server"\n');
fprintf(fid, ['ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de ', ...
    '''cd %s%s/; ', ...
    'sbatch -p large96 -t 12:00:00 -N 1 --tasks-per-node 96 -o %%j.output -e %%j.output --mem=747000mb --account bbkponsm --job-name %s ', ...
    '%s%s/jobparaview_global3D.sh''\n'],remote_Postprocess_scripts,output_directory,jobname,remote_Postprocess_scripts,output_directory);
fprintf(fid, '\n');
fprintf(fid, '# Wait for the job to finish\n');
fprintf(fid, 'while ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de ''squeue -n %s'' | grep -q %s; do\n',jobname,jobname);
fprintf(fid, '    echo "Job %s batch is running on the remote server. Waiting..."\n',jobname);
fprintf(fid, '    sleep 20  # Adjust the sleep interval\n');
fprintf(fid, 'done\n');
fprintf(fid, '\n');
fprintf(fid, 'echo "Paraview Job has finished. Proceeding to the copy of the data to local."\n');
fprintf(fid, '\n');
fprintf(fid, 'chmod +x /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, '/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, 'echo "Running Matlab Postprocess and making figures"\n');
fprintf(fid, '#/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "run(''/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; end; quit"\n',output_directory);
fprintf(fid, '/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run(''/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; disp(''Error: The MATLAB Postprocess encountered an issue.''); end; quit"\n',output_directory);
fprintf(fid, 'echo "Postprocess succeeded and figures extracted"\n');

fclose(fid);

disp(['Postprocess_auto.sh built and saved to: ' output_file_sh]);

%% make the sh file for restart without running pvbatch

output_file_sh = fullfile(output_directory, 'Postprocess_auto_restart.sh');

% Save the Bash script to a file using fprintf with explicit newline characters
fid = fopen(output_file_sh, 'w');
fprintf(fid, '#!/bin/bash\n\n');
fprintf(fid, '\n');
fprintf(fid,'# This script automatize the the extraction of the data and the postprocessing of the model to the making of the figures\n');
fprintf(fid, '\n');
fprintf(fid,'echo "Building sshfs connection of the remote server."\n');
fprintf(fid, '\n');
fprintf(fid, 'umount -f ~/Desktop/modelblogin\n');
fprintf(fid, '\n');
fprintf(fid, 'sshfs -o "IdentityFile=/Users/ponsm/.ssh/id_rsa_hlrn" bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm /Users/ponsm/Desktop/modelblogin -o defer_permissions\n');
fprintf(fid, '\n');
fprintf(fid, 'chmod +x /Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, '/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, 'echo "Running Matlab Postprocess and making figures"\n');
fprintf(fid, '#/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "run(''/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; end; quit"\n',output_directory);
fprintf(fid, '/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run(''/Users/ponsm/Nextcloud/Postdoc_MEET/Modelling/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; disp(''Error: The MATLAB Postprocess encountered an issue.''); end; quit"\n',output_directory);
fprintf(fid, 'echo "Postprocess succeeded and figures extracted"\n');

fclose(fid);

disp(['Postprocess_auto.sh built and saved to: ' output_file_sh]);

disp(['Go to : ', path_to_postprocess_ASPECT,output_directory, ' then run "chmod +x Postprocess_auto.sh" and "./Postprocess_auto.sh" to initiate the postprocess']);
disp(['or to update restart without running Pvbatch : run "chmod +x Postprocess_auto_restart.sh" and "./Postprocess_auto_restart.sh"']);




