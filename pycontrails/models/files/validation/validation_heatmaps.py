import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob

# Define the path to the CSV files
csv_files = glob.glob("/home/ktait98/pycontrails_kt/pycontrails/models/files/validation/*.csv")

# Initialize empty DataFrames for r_squared and nrmse
r_squared_df = pd.DataFrame()
nrmse_df = pd.DataFrame()

# Read each CSV file and append the data to the DataFrames
for file in csv_files:
    print(file)
    df = pd.read_csv(file, index_col=0)
    location = file.split('/')[-1].split('.')[0]  # Extract the location name from the file name
    r_squared_df[location] = df['r_squared']
    # print min and max r_squared
    print("min r_squared: ", df['r_squared'].min())
    print("max r_squared: ", df['r_squared'].max())

    nrmse_df[location] = df['nrmse']
    # print min and max nrmse but exclude inf values
    print("min nrmse: ", df['nrmse'].min())
    print("max nrmse: ", df['nrmse'].replace([float('inf')], 0).max())

# Plot the heatmap for r_squared
plt.figure(figsize=(10, 8))
sns.heatmap(r_squared_df, annot=False, linewidths=0, cmap="viridis")
plt.title("R² Heatmap")
plt.xlabel("Location")
plt.ylabel("Species")
plt.show()

# Plot the heatmap for nrmse
plt.figure(figsize=(10, 8))
sns.heatmap(nrmse_df, annot=False, linewidths=0, cmap="viridis")
plt.title("NRMSE Heatmap")
plt.xlabel("Location")
plt.ylabel("Species")
plt.show()