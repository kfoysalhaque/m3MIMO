import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load the CSV file
file_path = 'BER_8x8_beam_tracking3.csv'
data = pd.read_csv(file_path)

font_size = 55

# Clean up the data by renaming columns and dropping unnecessary rows
data.columns = ['Transmission', '-60', '-50', '-40', '-30', '-20', '-10', '0', '10', '20', '30', '40', '50', '60']
data = data.drop(0).reset_index(drop=True)

# Convert data to numeric where applicable
data.iloc[:, 1:] = data.iloc[:, 1:].apply(pd.to_numeric)

# Set up the plot
plt.figure(figsize=(12, 8))

# Define different line styles and markers for each technology
line_styles = ['-', '--', '-.', ':', '-', '--', '-.', ':']
markers = ['o', 's', 'D', '^', 'v', 'p', '*', 'x']
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange']

# Plot each line
for i, row in data.iterrows():
    plt.plot(data.columns[1:].astype(float), row[1:]*100, label=row['Transmission'],
             linestyle=line_styles[i % len(line_styles)],
             marker=markers[i % len(markers)], markersize=8, linewidth=3,
             color=colors[i % len(colors)])

# Beautify the plot
# plt.title('Multi-modal vs Single Modal Communication', fontsize=25)
plt.xlabel('RX location (degree)', fontsize=font_size)
plt.ylabel('BER (%)', fontsize=font_size)
plt.minorticks_on()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Set specific x-ticks
x_ticks = [-60, -40, -20, 0, 20, 40, 60]
plt.xticks(x_ticks, fontsize=font_size-5)

# Set specific y-ticks
y_ticks = [0, 10, 20, 30]
plt.yticks(y_ticks, fontsize=font_size-5)

plt.ylim(0, 30)
plt.legend(fontsize=font_size-10, loc='upper center', bbox_to_anchor=(0.5, 0.98), ncol=1, framealpha=0.8, columnspacing=0.3, handlelength=0.7)
plt.tight_layout()

# Save and show the plot
plt.savefig('BER_8x8_tracking.png', dpi=300, format='png')
plt.savefig('BER_8x8_tracking.eps', dpi=300, format='eps')
plt.show()