# Set common variables
remote_base_directory="/scratch/usr/bbkponsm/model/globalscale/Potsprocess/sphere3d/"
local_base_directory="/Volumes/Jerry/global_models_3d_extract/"
ssh_key="~/.ssh/id_rsa_hlrn"

# Loop through the list of model names
for model in "${@}"; do
    # Construct the remote and local paths
    remote_path="${remote_base_directory}${model}"
    local_path="${local_base_directory}"

    # Perform rsync
    rsync -varz -e "ssh -i $ssh_key" "bbkponsm@blogin.hlrn.de:$remote_path" "$local_path"
done
