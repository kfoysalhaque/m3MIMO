# Imports

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors
import matplotlib as mpl

# mpl.rcParams['font.serif'] = 'Palatino'
# mpl.rcParams['text.usetex'] = 'true'
# mpl.rcParams['text.latex.preamble'] = r'\usepackage{newtxmath}'
# mpl.rcParams['font.size'] = 22

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def add_labels(rects):
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{:.0f}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2 , height),
                    xytext=(1, 1),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=font_size-10)

# Change Params Here:

font_size = 55

file_string = 'gesture_accuracy_cross_study'
columns_as_group = ["Conference\nroom", "Lab"]

data_file = file_string + '.txt'
fig_pdf_file = file_string + '.pdf'
fig_eps_file = file_string + '.eps'

data = np.loadtxt(data_file)

#row_1 = data[:, :]
#row_2 = data[1, :]

x = np.arange(len(columns_as_group))
width = 0.15
fig, ax = plt.subplots()
fig.set_figheight(7)
fig.set_figwidth(10)

# Adjust the x values for tick positions
group_width = len(columns_as_group)

rects1 = ax.bar(x, data, width, color='skyblue', hatch='/', edgecolor="black", linewidth=2)
#rects2 = ax.bar(x + width, row_2, width, color='lightcyan', hatch=' \ ', edgecolor="black", linewidth=2)

add_labels(rects1)
#add_labels(rects2)

ax.set_ylabel('Accuracy (%)', fontsize=font_size)
ax.set_xticks(x)
ax.set_xticklabels(columns_as_group, fontsize=font_size)
# ax.set_xlabel('No. of sub-channels', fontsize=35)

ax.set_ylim(60, 105)
ax.set_yticks([70, 85, 100], [70, 85, 100], fontsize=font_size)

#ax.legend([rects1], rows_as_bar, loc='upper center', ncol=2, fontsize=30, framealpha=0.0)

plt.grid(axis='y', linestyle='--', linewidth=0.5)
# plt.title('Redundant channel occupation with one CSI monitor', fontsize= 35)

plt.tight_layout()

plt.savefig(fig_pdf_file, dpi=300, format='pdf')
plt.savefig(fig_eps_file, dpi=300, format='eps')

plt.show()
print('1')
