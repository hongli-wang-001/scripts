import pandas as pd

# Load the original CSV file
infile = 'your_input_file.csv'  # Replace with your actual input file name
outfile = 'your_output_file.csv'  # Replace with your desired output file name

# Read the CSV file
df = pd.read_csv(infile)

# Update the header
new_columns = [
    "time_stamp", "sensor_index", "humidity_a", "temperature_a", "pressure_a",
    "pm2.5_atm_a", "pm2.5_atm_b", "pm2.5_cf_1_a", "pm2.5_cf_1_b",
    "pm10.0_atm_a", "pm10.0_atm_b", "pm10.0_cf_1_a", "pm10.0_cf_1_b",
    "id", "location_type", "latitude", "longitude", "datetime_local",
    "datetime_utc", "PM25", "OZONE"
]

# Assign the new column names
df.columns = new_columns

# Write the updated DataFrame back to a CSV file
df.to_csv(outfile, index=False)

