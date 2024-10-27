# %%
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.ndimage import gaussian_filter
import pyshtools as pysh
import pygmt
from scipy.interpolate import griddata
import warnings
from cryptography.utils import CryptographyDeprecationWarning
from shapely.errors import ShapelyDeprecationWarning

with warnings.catch_warnings():
    warnings.filterwarnings('ignore', category=CryptographyDeprecationWarning)
    warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

    
#folder_path = '/Volumes/Jerry/global_models_3d/R01e_Rodinia_2GPa_Mantle_C20MPa_f003_LR/'
data_file_path = 
# '/Volumes/Jerry/global_models_3d_extract/R01e_Rodinia_2GPa_Mantle_C20MPa_f003_LR/'


#model_start_time = 1000  # [Ma]
lmax = 24
spacing = 180 / (2 * lmax + 2)
file_interval = 10

# Function to extract the numeric part of the filename
def extract_number(filename):
    file_number = int(filename.split('_')[-1].split('.')[0])
    return file_number

# List comprehension to get all filenames starting with 'depths_contours'
file_list = [filename for filename in os.listdir(data_file_path) if filename.startswith('depths_contours')]
sorted_file_list = sorted(file_list, key=extract_number)

# Ensure there are enough files in the directory
if len(sorted_file_list) < file_interval:
    raise ValueError(f"There are not enough files in the directory to select files at an interval of {file_interval}.")

# Select files at the specified interval
selected_files = sorted_file_list[::file_interval]

# Extract model name from folder_path
model_name = os.path.basename(os.path.normpath(data_file_path))

# Define the output directory
output_dir = f'/Users/ponsm/Nextcloud/group_monitoring_earth_evolution_through_time/Research/Michael_Pons/models/Global_model_3D/{model_name}/mantle_flow_spharmonics'
os.makedirs(output_dir, exist_ok=True)

# Initialize results list
all_results = []

# Open CSV file for writing
csv_filename = os.path.join(output_dir, 'depth_degree_power_data.csv')
with open(csv_filename, 'w') as csv_file:
    csv_file.write('File,Time (Myr),Depth (km),Degree,Power\n')

    # Loop over selected files and process each
    for file_idx, selected_file in enumerate(selected_files):
        print(f"Processing file {file_idx + 1}/{len(selected_files)}: {selected_file}")
        
        file_path = os.path.join(data_file_path, selected_file)

        data = pd.read_csv(file_path, delimiter=',', skiprows=1, 
                           names=['Time', 'Points:0', 'Points:1', 'Points:2', 'T', 'density', 'depth', 'viscosity'])
        data['r'] = np.sqrt(data['Points:0']**2 + data['Points:1']**2 + data['Points:2']**2)
        data['latitude'] = 90 - np.arccos(data['Points:2'] / data['r']) * 180 / np.pi
        data['longitude'] = 180 + np.arctan2(data['Points:1'], data['Points:0']) * 180 / np.pi

        unique_depths = data['depth'].unique()

        for depth_idx, test_depth in enumerate(unique_depths):
            print(f"  Processing depth {depth_idx + 1}/{len(unique_depths)}: {test_depth}")

            depth_data = data[data['depth'] == test_depth]
            latitudes = depth_data['latitude'].values
            longitudes = depth_data['longitude'].values
            viscosities = depth_data['viscosity'].values

            lon_grid, lat_grid = np.meshgrid(np.linspace(0.5, 359.5, 360), np.linspace(-89.5, 89.5, 180))
            viscosities_interp = griddata((longitudes, latitudes), viscosities, (lon_grid, lat_grid), method='linear')

            valid_mask = ~np.isnan(viscosities_interp)
            viscosities_interp_filled = griddata(
                (lon_grid[valid_mask], lat_grid[valid_mask]),
                viscosities_interp[valid_mask],
                (lon_grid, lat_grid),
                method='nearest'
            )

            viscosities_interp_filled_log = np.log10(viscosities_interp_filled)
            viscosities_vector = viscosities_interp_filled_log.flatten()
            lon_vector = lon_grid.flatten()
            lat_vector = lat_grid.flatten()

            data_in_columns = pd.DataFrame({
                'longitude': lon_vector,
                'latitude': lat_vector,
                'viscosities': viscosities_vector
            }).drop_duplicates(subset=['longitude', 'latitude']).dropna()

            data_grid = pygmt.sphinterpolate(data=data_in_columns, spacing=spacing, region="g")

            tlons = data_grid.coords['lon'].values
            tlats = data_grid.coords['lat'].values
            Lon, Lat = np.meshgrid(tlons, tlats)
            tdata = data_grid.values

            tgrid = pysh.SHGrid.from_array(np.flipud(tdata))
            cilm = tgrid.expand(lmax_calc=lmax, csphase=-1)
            power_per_l = cilm.spectrum()
            degrees = np.arange(lmax + 1, dtype=float)
            valid_indices = degrees > 0

            filtered_degrees = degrees[valid_indices]
            filtered_power_per_l = power_per_l[valid_indices]

            time = data.iat[0, data.columns.get_loc('Time')]
            time_str = f'{time / 1e6:.2f}'

            for deg, pwr in zip(filtered_degrees, filtered_power_per_l):
                all_results.append((time, test_depth, deg, pwr))
                csv_file.write(f'{selected_file},{time_str},{test_depth / 1000:.2f},{deg},{pwr}\n')

        # Plotting for each file
        df_plot = pd.DataFrame(all_results, columns=['time', 'depth', 'degree', 'power'])
        df_plot['depth'] = df_plot['depth'] / 1000  # Convert depth to kilometers
        df_plot['power2'] = np.sqrt(df_plot['power'])
        df_plot['power3'] = np.log10(df_plot['power'])

        # Create pivot table
        pivot_table = df_plot.pivot_table(index="depth", columns="degree", values="power3")

        # Apply Gaussian filter for smoothing
        smooth_pivot_table = gaussian_filter(pivot_table, sigma=1)  # Adjust sigma for desired smoothness

        # Set the colormap limits
        vmin, vmax = -4, -2  # Example limits; adjust as needed

        # Plotting
        plt.figure(figsize=(12, 8))
        sns.heatmap(smooth_pivot_table, cmap="Reds", cbar_kws={'label': 'log10(Power)'}, 
                    xticklabels=pivot_table.columns, yticklabels=pivot_table.index, vmin=vmin, vmax=vmax)

        # Set the title with the time included
        plt.title(f'SPH Power Spectrum by Depth and Degree (Time: {time_str} Myr)')
        plt.xlabel('Spherical Harmonic Degree')
        plt.ylabel('Depth [km]')

        # Customize x-axis to show integer labels
        ax = plt.gca()
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{int(x)}'))

        # Manually set y-axis ticks to match the depth values
        plt.yticks(ticks=np.arange(len(pivot_table.index)), labels=[f'{d:.1f}' for d in pivot_table.index], rotation=0)

        plt.gca()
        # .invert_yaxis()  # Optionally invert the y-axis to have depth increasing downward
        
        # Save the plot with a filename using leading zeros
        plot_filename = f'{model_name}_file_{file_idx + 1:03d}.png'
        plot_filepath = os.path.join(output_dir, plot_filename)
        plt.savefig(plot_filepath)
        plt.close()

        print(f"Saved plot to {plot_filepath}")


# %%

# Read the data from CSV
df_all_results = pd.read_csv(csv_filename)
print(df_all_results)

# df_all_results['Log_Power'] = np.log10(df_all_results['Power'])


# Get unique depths
unique_depths = df_all_results['Depth (km)'].unique()

# Create the output directory if it doesn't exist
output_subdir = os.path.join(output_dir, 'sph_per_depths')
os.makedirs(output_subdir, exist_ok=True)

# Function to plot the heatmap
def plot_heatmap(data, title, filename):
    plt.figure(figsize=(12, 8))
    
    # Create pivot table for the heatmap
    # pivot_table = data.pivot_table(index='Degree', columns='Time (Myr)', values='Log_Power', aggfunc='mean')
    pivot_table = data.pivot_table(index='Degree', columns='Time (Myr)', values='Power', aggfunc='mean')

    # Set the colormap limits
    # vmin, vmax = -4, 0  # Example limits; adjust as needed
    # Plot heatmap
    # sns.heatmap(pivot_table, cmap="Reds", cbar_kws={'label': 'log10(Power)'}, vmin=vmin, vmax=vmax)
    sns.heatmap(pivot_table, cmap="Reds", cbar_kws={'label': 'log10(Power)'})
    plt.title(title)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Harmonic Degree')
    plt.tight_layout()

    # Save the plot
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot to {filename}")

# # Plot for the surface (assuming surface is the shallowest depth)
# surface_depth = unique_depths.min()
# df_surface = df_all_results[df_all_results['Depth (km)'] == surface_depth]
# plot_heatmap(df_surface, f'Spherical Harmonic Degree vs. Time at Surface (Depth: {surface_depth:.2f} km)', 
#              os.path.join(output_dir, 'harmonic_degree_vs_time_surface.png'))

# Plot for each depth
for depth in unique_depths:
    df_depth = df_all_results[df_all_results['Depth (km)'] == depth]
    plot_heatmap(df_depth, f'Spherical Harmonic Degree vs. Time at Depth: {depth:.2f} km', 
                 os.path.join(output_subdir, f'harmonic_degree_vs_time_depth_{depth:.2f}_km.png'))

# %%



