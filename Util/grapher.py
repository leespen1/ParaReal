import numpy as np
import matplotlib.pyplot as plt

def get_durations(data_str):
    lines = data_str.strip(" \n").split("\n")
    return np.array([float(line.split(",")[0]) for line in lines])

def get_nd_pts(data_str):
    lines = data_str.strip(" \n").split("\n")
    return np.array([convert_nd_point_str(line.split(",")[1]) for line in lines])


#def process_data_str(data_str):
#    return np.array([process_nd_point_line(line) for line in data_str.split("\n")])
#
#def process_nd_point_line(line_str):
#    return np.array([convert_nd_point_str(entry) for entry in line_str.split(",")])
#
def convert_nd_point_str(pt_str):
    return np.array(pt_str.strip(" \n[]").split(";")).astype(float)
#
#def convert_results_section():
#    pass

data_list = ["", ""]
with open(input("Input File")) as infile:
    output_i = None
    for line in infile:
        if line.strip() == "ParaReal":
            output_i = 0
            continue
        elif line.strip() == "Serial":
            output_i = 1
            continue
        elif line.strip() == "Duration,Solution":
            continue
        elif output_i is not None:
            data_list[output_i] += line

parareal_data_str = data_list[0].strip(" \n")
serial_data_str = data_list[1].strip(" \n")

parareal_pts  = get_nd_pts(parareal_data_str)[1:]
parareal_durs = get_durations(parareal_data_str)[1:]

serial_pts  = get_nd_pts(serial_data_str)
serial_durs = get_durations(serial_data_str)

serial_norm = np.linalg.norm(serial_pts[0,:])
error_ary = np.array(
    [np.linalg.norm(serial_pts[0,:]-pt) for pt in parareal_pts]
)
frac_error_ary = np.array(
    [err/serial_norm for err in error_ary]
)


fig, ax = plt.subplots()
#ax.scatter(parareal_durs, error_ary)
ax.scatter(parareal_durs, frac_error_ary)
ax.set_title("Dev from Serial Sol vs Time Taken")
ax.set_xlabel("Time Taken (s)")
ax.set_ylabel("Fractional Deviation from Serial Solution")
ax.axvline(x=serial_durs[0], color='r', label="Serial Time")

machine_epsilon = np.finfo(float).eps
lower_perf_lim = None
for i, error in enumerate(error_ary):
    if error <= machine_epsilon:
        #ax.axvline(x=parareal_durs[i], color='g', 
        #           label="ParaReal Time - Perfect Accuracy Reached")
        lower_perf_lim = parareal_durs[i]
#ax.hlines(y=0.2, xmin=4, xmax=20, linewidth=2, color='r')
#ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_yscale("symlog") # Uses a linear scale near zero to allow showing points very close to zero. But this is giving me problems

xlim0, xlim1 = ax.get_xlim()
if lower_perf_lim is not None:
    ax.axvspan(lower_perf_lim, xlim1, alpha=0.2, color="g",
               label = "Machine-$\epsilon$ Accuracy")
ax.set_xlim(xlim0, max(lower_perf_lim, xlim1))

plt.legend()
plt.show()


## Example program output
#ParaReal
#Duration,Solution
#5.331e-06,[2.44140625;4.8828125;7.32421875]
#0.00126738,[2.706916392;5.413832784;8.120749177]
#0.002452182,[2.717744522;5.435489044;8.153233566]
#Serial
#Duration,Solution
#0.002422362,[2.717942121;5.435884242;8.153826363]

