import pandas as pd
import matplotlib.pyplot as plt

# Replace 'your_file.csv' with the path to your CSV file
csv_file = 'mu_excess.csv'

# Load the CSV data into a DataFrame
df = pd.read_csv(csv_file)

# Extract Acceptance Ratio and Density (rho) columns
acceptance_ratio = df['Acc_ratio']
density = df['rho']

# Create the plot
plt.figure(figsize=(10, 6))

# Plot Acceptance Ratio vs. Density
plt.plot(density, acceptance_ratio, marker='o', color='b', label='Acceptance Ratio vs. rho')

# Add labels and title
plt.xlabel('Density (rho)')
plt.ylabel('Acceptance Ratio')
plt.title('Acceptance Ratio vs. Density')
plt.legend()

# Show the plot
plt.grid(True)
plt.show()
