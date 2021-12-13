import numpy as np
import matplotlib.pyplot as plt
import re
import grapher_lib

color_dict = {
    1 : 'b',
    2 : 'g',
    4 : 'r',
    8 : 'm',
    16 : 'brown',
}
marker_dict = {
    1 : 'o',
    2 : 'v',
    4 : 's',
    8 : '*',
    16 : '+',
}

num_ranks_dict = {}
times_dict = {}

infile_template = "HPC_Runs_1153/Dahlquist_1D_Testing_Dev18_4/result_t%s_r1.log"
file_numbers = (1, 2, 4, 8, 16)
x_dict, y_dict, serial_times_dict = grapher_lib.get_graph_info(
    infile_template,
    file_numbers,
    "num ranks",
    "iteration",
    dim=1,
    tolerance=1e-14
)

fig, ax = plt.subplots(nrows=3, ncols=2)

fig.suptitle("Iterations By ParaReal Implementation")
fig.supxlabel("Number of Processes")
fig.supylabel("Avg # of Iterations to Converge")
fig.set_size_inches(8, 8)

# First subplot; all-in-one
for i in file_numbers:
    color = color_dict[i]
    marker = marker_dict[i]
    ax[0, 0].scatter(x_dict[i], y_dict[i], color=color, marker=marker, label="tf = %s"%i)
    ax[0, 0].set_title("All Timespans")

# Other subplots; one-set-each
for n, i in enumerate(file_numbers):
    color = color_dict[i]
    marker = marker_dict[i]
    ax[(1+n)//2, (1+n)%2].scatter(x_dict[i], y_dict[i], color=color, marker=marker)
    ax[(1+n)//2, (1+n)%2].set_title("$t_f$ = %s"%(i))

fig.legend(loc=2, prop={'size': 6}, ncol=1)
plt.tight_layout()
plt.show()
