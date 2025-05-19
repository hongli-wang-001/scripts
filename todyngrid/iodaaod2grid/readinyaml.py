import yaml

def create_config_yaml(filename='config.yaml'):
    config = {
        'input_file': 'pm25.nc',
        'obs_variable': 'ObsValue/particulatematter2p5Surface',
        'output_image': 'v10_interpolated_pm25.png',
        'output_netcdf': 'v10_interpolated_pm25.nc',
        'grid_file': 'dyn_grid_spec.nc'
    }

    with open(filename, 'w') as f:
        yaml.dump(config, f, sort_keys=False)

    print(f"Configuration file '{filename}' has been created.")

# Call the function
create_config_yaml()

