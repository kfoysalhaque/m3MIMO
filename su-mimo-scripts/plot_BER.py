import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# Load the CSV file
file_path = 'BER_8x8.csv'
data = pd.read_csv(file_path)

font_size = 45

# Clean up the data by renaming columns and dropping unnecessary rows
data.columns = ['Transmission', '1', '1.5', '2', '2.5', '3', '3.5', '4', '4.5', '5', '5.5', '6', '6.5', '7', '7.5', '8', '8.5', '9', '9.5', '10', '10.5', '11', '11.5', '12']
data = data.drop(0).reset_index(drop=True)

# Convert data to numeric where applicable
data.iloc[:, 1:] = data.iloc[:, 1:].apply(pd.to_numeric)

# Set up the plot
plt.figure(figsize=(10, 6))

# Define different line styles and markers for each technology
line_styles = ['-', '--', '-.', ':', '-', '--', '-.', ':']
markers = ['o', 's', 'D', '^', 'v', 'p', '*', 'x']
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'orange']

# Plot each line
for i, row in data.iterrows():
    plt.plot(data.columns[1:].astype(float), row[1:]*100, label=row['Transmission'],
             linestyle=line_styles[i % len(line_styles)],
             marker=markers[i % len(markers)], markersize=8, linewidth=2,
             color=colors[i % len(colors)])

# Beautify the plot
# plt.title('Multi-modal vs Single Modal Communication', fontsize=25)
plt.xlabel('Transmission distance (m)', fontsize=font_size)
plt.ylabel('BER (%)', fontsize=font_size)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Set specific x-ticks
x_ticks = [1, 3, 5, 7, 9, 11]
plt.xticks(x_ticks, fontsize=font_size-5)

# Set specific y-ticks
y_ticks = [0, 4, 8, 12, 16]
plt.yticks(y_ticks, fontsize=font_size-5)

plt.ylim(0, 16)
plt.legend(fontsize=font_size-10, loc='upper center', bbox_to_anchor=(0.4, 0.95), ncol=1, framealpha=0.0)
plt.tight_layout()

# Save and show the plot
plt.savefig('BER_8x8_su-mimo.png', dpi=300, format='png')
plt.savefig('BER_8x8_su-mimo.eps', dpi=300, format='eps')
plt.show()
