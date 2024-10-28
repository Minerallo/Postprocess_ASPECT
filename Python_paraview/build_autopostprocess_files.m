clc; clear all;
% This scripts get the model names and paths given by the user and automatically
% apply them to the Pvpython script and the mantlab postsprocess
path_to_postprocess_ASPECT = '/Users/ponsm/Desktop/Postprocess_ASPECT/';
paraview_python_script = [path_to_postprocess_ASPECT 'Python_paraview/Model_extract_global_3D_auto.py'];
paraview_sphpython_script = [path_to_postprocess_ASPECT 'Python_paraview/depths_contour_auto.py'];
postprocess_sphpython_script = [path_to_postprocess_ASPECT 'Python_paraview/sph_global.py'];

models_cluster_directory = '/scratch/usr/bbkponsm/model/globalscale/sphere3d/';
extracted_models_pvdata_directory = '/Volumes/Jerry/global_models_3d_extract/';
remote_Postprocess_scripts = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/postprocess_scripts/';
remote_Postprocess_directory = '/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/';
local_access_to_postprocess_directory = '/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/sphere3d/';
local_access_to_postprocess_scripts_directory='/Users/ponsm/Desktop/modelblogin/model/globalscale/Potsprocess/postprocess_scripts/';
local_access_to_models = '/Users/ponsm/Desktop/modelblogin/model/globalscale/sphere3d/';

%path of the Pvpython to modify
modified_file_path = './Model_extract_global_3D_auto.py';
modified_file_sphpath = './depths_contour_auto.py';
modified_file_sph_postprocess_path = './sph_global.py';

postprocess_choices = ""; %Modify as needed (call for spherical harmonic "sph", "fatbox", both "sph, fatbox", or leave empty)\n\n');

cd(path_to_postprocess_ASPECT);

%% List of specific model names to give
% model_names = {'R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003', ...
%     'R01f_Rodinia_2GPa_Mantle_C10MPa_f005_LR_SB_f003'};

model_names = {'V07d_R01f_Rodinia_init_RL9_ref_minvisc1e20_ContEx10_refn_F03_f0045_SW2em16_test'};
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


%% Read the content of the Python spherical harmonic script and modify the models list
fid = fopen(paraview_sphpython_script, 'r');
python_script_content = fread(fid, '*char')';
fclose(fid);

% Build the full paths to models by concatenating the base directory with model names
models_list_to_postprocess = cellfun(@(name) fullfile([models_cluster_directory name '/']), model_names, 'UniformOutput', false);

% Convert MATLAB cell array to Python list format
new_path_models_python = cellfun(@(x) ['''' x ''''], models_list_to_postprocess, 'UniformOutput', false);

% Replace the existing path_models in the Python script content
python_script_content = regexprep(python_script_content, 'path_models\s*=\s*\[.*?\]', ['path_models = [' strjoin(new_path_models_python, ',\n        ') ']\n']);

output_file = fullfile(output_directory, modified_file_sphpath);

% Save the modified script
fid = fopen(output_file, 'w');
fwrite(fid, python_script_content);
fclose(fid);

disp(['Modified Pvpython spharmonics script saved to: ' output_file]);


%% Read the content of the Python spherical harmonic postprocess script and modify the models list
fid = fopen(postprocess_sphpython_script, 'r');
python_script_content = fread(fid, '*char')';
fclose(fid);

% Build the full paths to models by concatenating the base directory with model names
models_list_to_postprocess = cellfun(@(name) fullfile([extracted_models_pvdata_directory name '/']), model_names, 'UniformOutput', false);

% Convert MATLAB cell array to Python list format
new_path_models_python = cellfun(@(x) ['''' x ''''], models_list_to_postprocess, 'UniformOutput', false);

% Replace the existing path_models in the Python script content
% python_script_content = regexprep(python_script_content, 'data_file_path\s*=\s*\.*?\', ['data_file_path = ' strjoin(new_path_models_python, ',\n        ') ']\n']);
split_content = strsplit(python_script_content, newline);
for i = 1:length(split_content)
    if contains(split_content{i}, 'data_file_path =')
        split_content{i} = ['data_file_path = ' strjoin(new_path_models_python, ',\n        ') ''];
        break;
    end
end

python_script_content = strjoin(split_content, newline);


%% Read the content of the Python Fatbox scripts postprocess script and modify the models list
output_file = fullfile(output_directory, modified_file_sph_postprocess_path);
% Save the modified script
fid = fopen(output_file, 'w');
fwrite(fid, python_script_content);
fclose(fid);

disp(['Modified Pvpython spharmonics script saved to: ' output_file]);


% Define the base path where the models are stored
base_folder_path = extracted_models_pvdata_directory;
% if solution surfaces are saved in the external volumes then we can use 
% '/Volumes/Jerry/global_models_3d/';
% otherwise access them via sshfs, later I may do it on the cluster directly

% Define the base output directory
output_directory_fatbox_plate_boundaries = '/Users/ponsm/Nextcloud/group_monitoring_earth_evolution_through_time/Research/Michael_Pons/models/Global_model_3D/';
postprocess_fatbox_python_script = '/Users/ponsm/Desktop/software/fatbox/spherical_scripts/scripts_for_autopostprocess_ASPECT/';

% Define the list of Python scripts to modify
python_scripts = {...
    '1_Plate_boundaries_extract.py', ...
    '2_Plate_boundaries_correlation.py', ...
    '3_Plate_boundaries_slip.py', ...
    '4_Plate_boundaries_displacement.py', ...
    '5_Plate_boundaries_evolution.py'};

% Loop through each model name
for i = 1:length(model_names)
    
    % Define the full folder path for the current model
    model_folder_path = fullfile(base_folder_path, model_names{i});
    
    % Define the output directory for this model
    model_output_directory = fullfile(output_directory);
    
    % Ensure the directory exists
    if ~exist(model_output_directory, 'dir')
        mkdir(model_output_directory);
    end
    
    % Define the output directory for fatbox plate boundaries specific to this model
    specific_output_directory = fullfile(output_directory_fatbox_plate_boundaries, model_names{i}, 'plate_boundaries/');
    
    % Loop through each Python script to modify
    for j = 1:length(python_scripts)
        
        % Print which Python script is being modified
        disp(['Modifying script: ' python_scripts{j} ' for model: ' model_names{i}]);
        
        % Read the content of the current Python script
        script_path = fullfile(postprocess_fatbox_python_script, python_scripts{j});
        fid = fopen(script_path, 'r');
        python_script_content = fread(fid, '*char')';
        fclose(fid);
        
        % Replace the existing 'folder_path' in the Python script content
        split_content = strsplit(python_script_content, newline);
        for k = 1:length(split_content)
            if contains(split_content{k}, 'folder_path =')
                split_content{k} = ['folder_path = ''' model_folder_path '/'''];
            elseif contains(split_content{k}, 'output_directory =')
                split_content{k} = ['output_directory = ''' specific_output_directory ''''];
            end
        end
        
        % Reconstruct the modified Python script content
        python_script_content = strjoin(split_content, newline);
        
        % Define the output file path for the modified Python script
        output_file = fullfile(model_output_directory, python_scripts{j});
        
        % Save the modified script
        fid = fopen(output_file, 'w');
        fwrite(fid, python_script_content);
        fclose(fid);
        
        disp(['Modified Python Fatbox script for model ' model_names{i} ' saved to: ' output_file]);
    end
end

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
fprintf(fid, 'postprocess_choice="%s";  # Modify as needed (e.g., "sph", "fatbox", "sph, fatbox", or leave empty)\n\n', postprocess_choices);
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
    fprintf(fid, 'rsync -varz -e ''ssh -i %s'' --include="*statistics*" --exclude="*" bbkponsm@blogin.hlrn.de:%s%s/ %s%s/.\n', ssh_key, models_cluster_directory, model_names{i}, local_base_directory, model_names{i});
    fprintf(fid, 'rsync -varz -e ''ssh -i %s'' --include="depth_average.txt" --exclude="*" bbkponsm@blogin.hlrn.de:%s%s/ %s%s/.\n', ssh_key, models_cluster_directory, model_names{i}, local_base_directory, model_names{i});
    fprintf(fid, 'rsync -varz -e ''ssh -i %s'' --include="original.prm" --exclude="*" bbkponsm@blogin.hlrn.de:%s%s/ %s%s/.\n', ssh_key, models_cluster_directory, model_names{i}, local_base_directory, model_names{i});
    fprintf(fid, 'rsync -varz --ignore-existing -e ''ssh -i %s'' bbkponsm@blogin.hlrn.de:%s%s/ %s%s/.\n',ssh_key,remote_base_directory, model_names{i}, local_base_directory, model_names{i});
%     fprintf(fid, 'rsync -varz -e ''ssh -i %s'' --include="solution_surface-*.vtu" --exclude="*" bbkponsm@blogin.hlrn.de:%s%s/solution/ %s%s/.\n',ssh_key,models_cluster_directory, model_names{i}, local_base_directory, model_names{i});
    % Conditionally copy solution_surface files only if "fatbox" is in postprocess_choice
    fprintf(fid, 'if [[ " ${postprocess_choice[@]} " =~ " fatbox " ]]; then\n');
    fprintf(fid, '    rsync -varz -e ''ssh -i %s'' --include="solution_surface-*.vtu" --exclude="*" bbkponsm@blogin.hlrn.de:%s%s/solution/ %s%s/.\n', ssh_key, models_cluster_directory, model_names{i}, local_base_directory, model_names{i});
    fprintf(fid, 'fi\n');
end

fclose(fid);

disp(['Bash_copy script saved to: ' output_file_bash_copy]);

%% make a jobscript for Pv batch in the autogenerated folder.sh

output_pvbatch_jobscript = fullfile(output_directory, 'jobparaview_global3D.sh');

% Save the Bash script to a file using fprintf with explicit newline characters
fid_slurm = fopen(output_pvbatch_jobscript, 'w');
%%%%RL9
fprintf(fid_slurm, '#!/bin/bash\n\n');
fprintf(fid_slurm, '#SBATCH -A bbk00014\n');
fprintf(fid_slurm, 'module load HRZIBenv sw.clx.el9 slurm anaconda3/2023.09\n\n');
fprintf(fid_slurm, 'module unload openmpi/gcc/5.0.3\n\n');
fprintf(fid_slurm, 'module load impi\n\n');
fprintf(fid_slurm, '\n');
fprintf(fid_slurm, '\n');
fprintf(fid_slurm, '/sw/comm/impi/mpi/2021.13/bin/mpirun -np 96 /sw/viz/paraview/x86_64.el9/ParaView-headless-5.13.1-egl-MPI-Linux-Python3.10-x86_64/bin/pvbatch %s%s/Model_extract_global_3D_auto.py\n',remote_Postprocess_scripts,output_directory);
fclose(fid_slurm);
%%%CentOS7
% fprintf(fid_slurm, '#!/bin/bash\n\n');
% fprintf(fid_slurm, '#SBATCH -A bbk00014\n');
% fprintf(fid_slurm, 'module unload gcc/9.3.0\n\n');
% fprintf(fid_slurm, 'module load gcc/9.2.0 openmpi/gcc.9 anaconda3 llvm/9.0.0 paraview\n\n');
% fprintf(fid_slurm, '# Ensure the cpus-per-task option is propagated to srun commands\n');
% fprintf(fid_slurm, 'export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK\n\n');
% fprintf(fid_slurm, 'mpirun --map-by socket:pe=$OMP_NUM_THREADS pvbatch %s%s/Model_extract_global_3D_auto.py\n',remote_Postprocess_scripts,output_directory);
% fclose(fid_slurm);

disp(['Your Slurm script built and saved to: ' output_pvbatch_jobscript]);


%% make a jobscript for Pv batch in the autogenerated folder.sh for contour_dephs

output_pvbatch_jobscript = fullfile(output_directory, 'jobparaview_global3D_depths_contour.sh');

% Save the Bash script to a file using fprintf with explicit newline characters
fid_slurm = fopen(output_pvbatch_jobscript, 'w');
%%%%RL9
fprintf(fid_slurm, '#!/bin/bash\n\n');
fprintf(fid_slurm, '#SBATCH -A bbk00014\n');
fprintf(fid_slurm, 'module load HRZIBenv sw.clx.el9 slurm anaconda3/2023.09\n\n');
fprintf(fid_slurm, 'module unload openmpi/gcc/5.0.3\n\n');
fprintf(fid_slurm, 'module load impi\n\n');
fprintf(fid_slurm, '\n');
fprintf(fid_slurm, '\n');
fprintf(fid_slurm, '/sw/comm/impi/mpi/2021.13/bin/mpirun -np 96 /sw/viz/paraview/x86_64.el9/ParaView-headless-5.13.1-egl-MPI-Linux-Python3.10-x86_64/bin/pvbatch %s%s/depths_contour_auto.py\n',remote_Postprocess_scripts,output_directory);
fclose(fid_slurm);

%%%%old 2023 CentOS7
% fprintf(fid_slurm, '#!/bin/bash\n\n');
% fprintf(fid_slurm, '#SBATCH -A bbk00014\n');
% fprintf(fid_slurm, 'module unload gcc/9.3.0\n\n');
% fprintf(fid_slurm, 'module load gcc/9.2.0 openmpi/gcc.9 anaconda3 llvm/9.0.0 paraview\n\n');
% fprintf(fid_slurm, '# Ensure the cpus-per-task option is propagated to srun commands\n');
% fprintf(fid_slurm, 'export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK\n\n');
% fprintf(fid_slurm, 'mpirun --map-by socket:pe=$OMP_NUM_THREADS pvbatch %s%s/depths_contour_auto.py\n',remote_Postprocess_scripts,output_directory);
% fclose(fid_slurm);

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


%% make the sh file for spharmonics
output_file_sh = fullfile(output_directory, 'spharmonics.sh');

% Save the Bash script to a file using fprintf with explicit newline characters
fid = fopen(output_file_sh, 'w');

fprintf(fid, '#!/bin/bash\n\n');
fprintf(fid, '# Name of the Conda environment\n');
fprintf(fid, 'ENV_NAME="extra_postprocess"\n\n');
fprintf(fid, '# Path to your Python script\n');
fprintf(fid, 'PYTHON_SCRIPT="sph_global.py"\n\n');
fprintf(fid, '# Activate the Conda environment\n');
fprintf(fid, 'source "$(conda info --base)/etc/profile.d/conda.sh"\n');
fprintf(fid, 'conda activate $ENV_NAME\n\n');
fprintf(fid, '# Run the Python script\n');
fprintf(fid, 'python $PYTHON_SCRIPT\n\n');
fprintf(fid, '# Deactivate the Conda environment\n');
fprintf(fid, 'conda deactivate\n');

fclose(fid);

disp(['spharmonics.sh built and saved to: ' output_file_sh]);

%% make the sh file for Fatbox
% Specify the output file for the Bash script
output_file_sh = fullfile(output_directory, 'Fatbox_plate_boundaries.sh');

% Save the Bash script to a file using fprintf with explicit newline characters
fid = fopen(output_file_sh, 'w');

fprintf(fid, '#!/bin/bash\n\n');
fprintf(fid, '# Name of the Conda environment\n');
fprintf(fid, 'ENV_NAME="extra_postprocess"\n\n');
fprintf(fid, '# List of Python scripts to execute\n');
fprintf(fid, 'PYTHON_SCRIPTS=(\n');
fprintf(fid, '    "1_Plate_boundaries_extract.py" \\\n');
fprintf(fid, '    "2_Plate_boundaries_correlation.py" \\\n');
fprintf(fid, '    "3_Plate_boundaries_slip.py" \\\n');
fprintf(fid, '    "4_Plate_boundaries_displacement.py" \\\n');
fprintf(fid, '    "5_Plate_boundaries_evolution.py"\n');
fprintf(fid, ')\n\n');
fprintf(fid, '# Activate the Conda environment\n');
fprintf(fid, 'source "$(conda info --base)/etc/profile.d/conda.sh"\n');
fprintf(fid, 'conda activate $ENV_NAME\n\n');
fprintf(fid, '# Loop through and run each Python script\n');
fprintf(fid, 'for SCRIPT in "${PYTHON_SCRIPTS[@]}"; do\n');
fprintf(fid, '    python "$SCRIPT"\n');
fprintf(fid, 'done\n\n');
fprintf(fid, '# Deactivate the Conda environment\n');
fprintf(fid, 'conda deactivate\n');

fclose(fid);

disp(['Fatbox_plate_boundaries.sh built and saved to: ' output_file_sh]);

%% make the sh file that will run all the scripts (could be more flexible in the future)

output_file_sh = fullfile(output_directory, 'Postprocess_auto.sh');
if numel(model_list) > 3
    last_letters = model_list(end-3:end);
else
    last_letters = model_list;
end

jobname =['Pv',last_letters];
jobname_sph =['Pvs',last_letters];

% Save the Bash script to a file using fprintf with explicit newline characters
fid = fopen(output_file_sh, 'w');
fprintf(fid, '#!/bin/bash\n\n');
fprintf(fid, '\n');
fprintf(fid,'# This script automatize the the extraction of the data and the postprocessing of the model to the making of the figures\n');
fprintf(fid, '\n');
% Define postprocess_choice as a comma-separated list of options or leave empty for no post-process
fprintf(fid, 'postprocess_choice="%s";  # Modify as needed (e.g., "sph", "fatbox", "sph, fatbox", or leave empty)\n\n', postprocess_choices);
fprintf(fid, '\n');
fprintf(fid,'echo "Building sshfs connection of the remote server."\n');
fprintf(fid, '\n');
fprintf(fid, '#umount -f ~/Desktop/modelblogin\n');
fprintf(fid, '\n');
fprintf(fid, 'sshfs -o "IdentityFile=/Users/ponsm/.ssh/id_rsa_hlrn" bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm /Users/ponsm/Desktop/modelblogin -o defer_permissions\n');
fprintf(fid, '\n');
fprintf(fid, 'echo "Building and sending files to the server"\n');
fprintf(fid, '##/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run(''/Users/ponsm/Desktop/Postprocess_ASPECT/%s/build_autopostprocess_files.m''); catch;disp(''Error: The MATLAB auto-postprocess files could not be properly built.''); end; quit"\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, 'rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Model_extract_global_3D_auto.py bbkponsm@blogin.hlrn.de:%s%s/.\n',output_directory,remote_Postprocess_scripts,output_directory);
fprintf(fid, 'rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Desktop/Postprocess_ASPECT/%s/jobparaview_global3D.sh bbkponsm@blogin.hlrn.de:%s%s/.\n',output_directory,remote_Postprocess_scripts,output_directory);
fprintf(fid, 'rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Desktop/Postprocess_ASPECT/%s/depths_contour_auto.py bbkponsm@blogin.hlrn.de:%s%s/.\n',output_directory,remote_Postprocess_scripts,output_directory);
fprintf(fid, 'rsync -varz -e "ssh -i /Users/ponsm/.ssh/id_rsa_hlrn" /Users/ponsm/Desktop/Postprocess_ASPECT/%s/jobparaview_global3D_depths_contour.sh bbkponsm@blogin.hlrn.de:%s%s/.\n',output_directory,remote_Postprocess_scripts,output_directory);
fprintf(fid, '\n');
fprintf(fid, 'echo "Running Pvbatch %s on the remote server"\n',jobname);
fprintf(fid, ['ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de ', ...
    '''cd %s%s/; ', ...
    '/usr/local/slurm/bin/sbatch -p cpu-clx:large -t 12:00:00 -N 1 --tasks-per-node 96 -o %%j.output -e %%j.error --mem=747000mb --account bbkponsm --job-name %s ', ...
    '%s%s/jobparaview_global3D.sh''\n'],remote_Postprocess_scripts,output_directory,jobname,remote_Postprocess_scripts,output_directory);
% Add conditional check for "sph" in postprocess_choice
fprintf(fid, 'if [[ "$postprocess_choice" == *"sph"* ]]; then\n');
fprintf(fid, '    echo "Running Pvbatch %s on the remote server for depth contour processing"\n', jobname_sph);
fprintf(fid, ['    ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de ', ...
    '''cd %s%s/; ', ...
    '/usr/local/slurm/bin/sbatch -p cpu-clx:large -t 12:00:00 -N 1 --tasks-per-node 96 -o %%j.output -e %%j.error --mem=747000mb --account bbkponsm --job-name %s ', ...
    '%s%s/jobparaview_global3D_depths_contour.sh''\n'], remote_Postprocess_scripts, output_directory, jobname_sph, remote_Postprocess_scripts, output_directory);
fprintf(fid, 'fi\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '# Wait for the job to finish\n');
fprintf(fid, 'while ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de ''/usr/local/slurm/bin/squeue -n %s'' | grep -q %s; do\n',jobname,jobname);
fprintf(fid, '    echo "Job %s batch is running on the remote server. Waiting..."\n',jobname);
fprintf(fid, '    sleep 20  # Adjust the sleep interval\n');
fprintf(fid, 'done\n');
fprintf(fid, '\n');
fprintf(fid, 'echo "Paraview Job has finished. Proceeding to the copy of the data to local."\n');
fprintf(fid, '\n');
fprintf(fid, 'chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, 'echo "Running Matlab Postprocess and making figures"\n');
fprintf(fid, '#/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "run(''/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; end; quit"\n',output_directory);
fprintf(fid, '/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run(''/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; disp(''Error: The MATLAB Postprocess encountered an issue.''); end; quit"\n',output_directory);

fprintf(fid, '\n');

% Check if postprocess_choice is empty
fprintf(fid, 'if [ -z "$postprocess_choice" ]; then\n');
fprintf(fid, '    echo "No postprocess selected. Skipping postprocess steps."\n');
fprintf(fid, 'else\n');

% Split postprocess_choice into an array
fprintf(fid, '    IFS=", " read -r -a choices <<< "$postprocess_choice"\n');

% Check if "sph" is one of the choices and, if so, wait for the job and copy models
fprintf(fid, '    if [[ " ${choices[@]} " =~ " sph " ]]; then\n');
fprintf(fid, '        while ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de ''/usr/local/slurm/bin/squeue -n %s'' | grep -q %s; do\n', jobname_sph, jobname_sph);
fprintf(fid, '            echo "Job %s batch is running on the remote server. Waiting..."\n', jobname_sph);
fprintf(fid, '            sleep 20  # Adjust the sleep interval\n');
fprintf(fid, '        done\n\n');

% Copy models from remote if "sph" is selected
fprintf(fid, '        /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n', output_directory);
fprintf(fid, '    fi\n\n');

% Loop through each choice in postprocess_choice and execute relevant scripts
fprintf(fid, '    for choice in "${choices[@]}"; do\n');

% Run spherical harmonics post-process if "sph" is selected
fprintf(fid, '        if [ "$choice" = "sph" ]; then\n');
fprintf(fid, '            echo "Running spherical harmonics postprocess. Ensure the required Conda environment is loaded."\n');
fprintf(fid, '            echo "If the environment \\"extra_postprocess\\" does not exist, create it using the following commands:"\n');
fprintf(fid, '            echo "  conda create --name extra_postprocess python=3.10.9"\n');
fprintf(fid, '            echo "  conda activate extra_postprocess"\n');
fprintf(fid, '            echo "  conda install -c conda-forge cartopy multiprocess cryptography cryptography-vectors pyvista pyshtools meshio dask geopy"\n\n');
fprintf(fid, '            chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n', output_directory);
fprintf(fid, '            /Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n', output_directory);
fprintf(fid, '        fi\n\n');

% Run fatbox post-process if "fatbox" is selected
fprintf(fid, '        if [ "$choice" = "fatbox" ]; then\n');
fprintf(fid, '            echo "Running Fatbox_plate_boundaries.sh"\n');
fprintf(fid, '            chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
fprintf(fid, '            /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
fprintf(fid, '        fi\n\n');

% End of loop
fprintf(fid, '    done\n\n');
fprintf(fid, '    echo "Postprocess succeeded and figures extracted"\n');
fprintf(fid, 'fi\n');  



%%CentOS7
% fprintf(fid, 'while ssh -i /Users/ponsm/.ssh/id_rsa_hlrn bbkponsm@blogin.hlrn.de ''squeue -n %s'' | grep -q %s; do\n',jobname_sph,jobname_sph);
% fprintf(fid, '    echo "Job %s batch is running on the remote server. Waiting..."\n',jobname_sph);
% fprintf(fid, '    sleep 20  # Adjust the sleep interval\n');
% fprintf(fid, 'done\n');
% 
% fprintf(fid, '\n');
% fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
% fprintf(fid, '\n');
% 
% fprintf(fid, ['echo "Running spherical harmonics postprocess. Ensure the required Conda environment is loaded. ' ...
%     'If the environment \\"extra_postprocess\\" does not exist, create it using the following commands:"\n']);
% fprintf(fid, 'echo "  conda create --name extra_postprocess python=3.10.9"\n');
% fprintf(fid, 'echo "  conda activate extra_postprocess"\n');
% fprintf(fid, 'echo "  conda install -c conda-forge cartopy multiprocess cryptography cryptography-vectors pyvista pyshtools meshio dask geopy"\n');
% fprintf(fid, '\n');
% fprintf(fid, 'chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n',output_directory);
% fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n',output_directory);
% fprintf(fid, '\n');
% fprintf(fid, 'echo "Running Fatbox_plate_boundaries.sh"\n');
% fprintf(fid, '\n');
% fprintf(fid, 'chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
% fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
% fprintf(fid, '\n');
% fprintf(fid, 'echo "Postprocess succeeded and figures extracted"\n');


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
fprintf(fid, '#umount -f ~/Desktop/modelblogin\n');
fprintf(fid, '\n');
fprintf(fid, 'sshfs -o "IdentityFile=/Users/ponsm/.ssh/id_rsa_hlrn" bbkponsm@blogin.hlrn.de:/scratch/usr/bbkponsm /Users/ponsm/Desktop/modelblogin -o defer_permissions\n');
fprintf(fid, '\n');
fprintf(fid, 'chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
fprintf(fid, '\n');
fprintf(fid, 'echo "Running Matlab Postprocess and making figures"\n');
fprintf(fid, '#/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "run(''/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; end; quit"\n',output_directory);
fprintf(fid, '/Applications/MATLAB_R2022b.app/bin/matlab -nodisplay -r "try; run(''/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Postprocess_auto_generated''); catch; disp(''Error: The MATLAB Postprocess encountered an issue.''); end; quit"\n',output_directory);
fprintf(fid, '\n');

% Split postprocess_choice into an array
fprintf(fid, 'postprocess_choice="%s";  # Modify as needed (e.g., "sph", "fatbox", "sph, fatbox", or leave empty)\n\n', postprocess_choices);

fprintf(fid, 'IFS=", " read -r -a choices <<< "$postprocess_choice"\n\n');

% Loop through each choice in postprocess_choice and execute relevant scripts
fprintf(fid, 'for choice in "${choices[@]}"; do\n');

% Run spherical harmonics post-process if "sph" is selected
fprintf(fid, '    if [ "$choice" = "sph" ]; then\n');
fprintf(fid, '        echo "Running spherical harmonics postprocess. Ensure the required Conda environment is loaded."\n');
fprintf(fid, '        echo "If the environment \\"extra_postprocess\\" does not exist, create it using the following commands:"\n');
fprintf(fid, '        echo "  conda create --name extra_postprocess python=3.10.9"\n');
fprintf(fid, '        echo "  conda activate extra_postprocess"\n');
fprintf(fid, '        echo "  conda install -c conda-forge cartopy multiprocess cryptography cryptography-vectors pyvista pyshtools meshio dask geopy"\n\n');
fprintf(fid, '        chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n', output_directory);
fprintf(fid, '        /Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n', output_directory);
fprintf(fid, '    fi\n\n');

% Run fatbox post-process if "fatbox" is selected
fprintf(fid, '    if [ "$choice" = "fatbox" ]; then\n');
fprintf(fid, '        echo "Running Fatbox_plate_boundaries.sh"\n');
fprintf(fid, '        chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
fprintf(fid, '        /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
fprintf(fid, '    fi\n\n');

% End of loop
fprintf(fid, 'done\n\n');
fprintf(fid, 'echo "Postprocess succeeded and figures extracted"\n');

%%%CentOS7
% fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Models_to_copy_from_remote.sh\n',output_directory);
% fprintf(fid, '\n');
% fprintf(fid, ['echo "Running spherical harmonics postprocess. Ensure the required Conda environment is loaded. ' ...
%     'If the environment \\"extra_postprocess\\" does not exist, create it using the following commands:"\n']);
% fprintf(fid, 'echo "  conda create --name extra_postprocess python=3.10.9"\n');
% fprintf(fid, 'echo "  conda activate extra_postprocess"\n');
% fprintf(fid, 'echo "  conda install -c conda-forge cartopy multiprocess cryptography cryptography-vectors pyvista pyshtools meshio dask"\n');
% fprintf(fid, '\n');
% fprintf(fid, 'chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n',output_directory);
% fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/spharmonics.sh\n',output_directory);
% fprintf(fid, '\n');
% fprintf(fid, 'echo "Postprocess succeeded and figures extracted"\n');
% fprintf(fid, 'echo "Running Fatbox_plate_boundaries.sh"\n');
% fprintf(fid, '\n');
% fprintf(fid, 'chmod +x /Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
% fprintf(fid, '/Users/ponsm/Desktop/Postprocess_ASPECT/%s/Fatbox_plate_boundaries.sh\n', output_directory);
% fprintf(fid, '\n');
% fprintf(fid, 'echo "Postprocess succeeded and figures extracted"\n');


fclose(fid);

disp(['Postprocess_auto.sh built and saved to: ' output_file_sh]);

disp(['Go to : ', path_to_postprocess_ASPECT,output_directory, ' then run "chmod +x Postprocess_auto.sh" and "./Postprocess_auto.sh" to initiate the postprocess']);
disp(['or to update restart without running Pvbatch : run "chmod +x Postprocess_auto_restart.sh" and "./Postprocess_auto_restart.sh"']);




