import numpy as np
import matplotlib.pyplot as plt
import re

def get_durations(data_str):
    lines = data_str.strip(" \n").split("\n")
    return np.array([float(line.split(",")[0]) for line in lines])


def get_nd_pts(data_str):
    lines = data_str.strip(" \n").split("\n")
    return np.array([convert_nd_point_str(line.split(",")[1]) for line in lines])


def convert_nd_point_str(pt_str):
    return np.array(pt_str.strip(" \n[]").split(";")).astype(float)


"""
Given a section, return the parareal and serial data (if any). Results WILL
include coarse solve, which the user is likely uninterested in.
"""
def get_results_str_from_section(section_lines):

    data_list = ["", ""]
    output_i = None
    for line in section_lines:
        if line.strip() == "ParaReal":
            output_i = 0
            continue
        elif line.strip() == "Serial":
            output_i = 1
            continue
        elif line.strip() == "Duration,Solution":
            continue
        elif line.strip() == "End Run":
            continue
        elif output_i is not None:
            data_list[output_i] += line.strip("\n") + "\n" # Make sure has 1 newline

    return data_list[0].strip(" \n"), data_list[1].strip(" \n") # Make sure no lingering whitespace


def get_section_strs(file_name):
    section_str = ""
    with open(file_name) as infile:
        for line in infile:
            if line.strip() == "Start Run":
                section_str = ""
            elif line.strip() == "End Run":
                yield section_str
                section_str = ""
            else:
                section_str += line

def get_serial_time_and_sol_from_results_str(results_str):
    return get_durations(results_str)[0], get_nd_pts(results_str)[0] # Serial solution should only be one line


def get_parareal_data_from_results_str(results_str, serial_sol, include_coarse=False):
    results_str = results_str.strip(" \n")

    start_i = 0 if include_coarse else 1

    parareal_pts  = get_nd_pts(results_str)[start_i:]
    parareal_durs = get_durations(results_str)[start_i:]

    serial_norm = np.linalg.norm(serial_sol)
    error_ary = np.array(
        [np.linalg.norm(serial_sol - pt) for pt in parareal_pts]
    )
    frac_error_ary = np.array(
        [err/serial_norm for err in error_ary]
    )
    return np.column_stack((parareal_durs, frac_error_ary))


def get_time_to_meet_tolerance(parareal_data, tolerance=None):
    if tolerance is None:
        tolerance = np.finfo(float).eps # Use machine epsilon if nothing else available
    for dur, frac_err in parareal_data:
        if frac_err <= tolerance:
            return dur
    return -1


def get_iter_to_meet_tolerance(parareal_data, tolerance=None, coarse_included=False):
    if tolerance is None:
        tolerance = np.finfo(float).eps # Use machine epsilon if nothing else available
    for i, (dur, frac_err) in enumerate(parareal_data):
        if frac_err <= tolerance:
            return i+1 # Account for zero-indexing
    return -1


def get_section_info_from_str_1d(section_str):
    section_info = {}
    section_info["lambda"] = re.search(r"(?<=u' = )[-\d.]+", section_str).group(0)
    section_info["num ranks"] = re.search(r"(?<=Num Ranks: )[\d]+", section_str).group(0)
    section_info["u0"] = re.search(r"(?<=u0 = )[-\d.]+", section_str).group(0)
    section_info["t0"] = re.search(r"(?<=t0 = )[-\d.]+", section_str).group(0)
    section_info["tf"] = re.search(r"(?<=tf = )[-\d.]+", section_str).group(0)
    section_info["solvers per rank"] = re.search(r"(?<=Fine Solvers Per Rank: )[\d]+", section_str).group(0)

    return section_info


def get_section_info_from_str_2d_polar(section_str):
    # Needs to be tested
    section_info = {}
    #section_info["lambda"] = re.search(r"(?<=u' = )[-\d.]+", section_str).group(0)
    section_info["num ranks"] = re.search(r"(?<=Num Ranks: )[\d]+", section_str).group(0)
    section_info["u0"] = re.search(r"(?<=u0 = )[-\d.()]+", section_str).group(0)
    section_info["t0"] = re.search(r"(?<=t0 = )[-\d.]+", section_str).group(0)
    section_info["tf"] = re.search(r"(?<=tf = )[-\d.]+", section_str).group(0)
    section_info["r"] = re.search(r"(?<=r = )[-\d.]+", section_str).group(0)
    section_info["theta"] = re.search(r"(?<=theta = )[-\d.]+", section_str).group(0)
    section_info["solvers per rank"] = re.search(r"(?<=Fine Solvers Per Rank: )[\d]+", section_str).group(0)

    return section_info


def get_graph_info(infile_template, file_numbers, x_type, y_type, dim=1, tolerance=1e-14):
    if y_type not in ("time", "iteration"):
        raise ValueError("Unexpected Graph y_type")
    if x_type not in ("num ranks", "theta"):
        raise ValueError("Unexpected Graph x_type")
    if dim not in (1, 2):
        raise ValueError("Unexpected Dimension")


    x_dict = {}
    y_dict = {}
    serial_times_dict = {}

    serial_time = None
    serial_sol = None
    section_info = None
    for i in file_numbers:
        infile_name = infile_template%(i)
        x_dict[i] = []
        y_dict[i] = []
        for section_str in get_section_strs(infile_name):
            parareal_str, serial_str = get_results_str_from_section(section_str.split("\n"))
            # Should really only do this once, from the first run
            if serial_str:
                serial_time, serial_sol = get_serial_time_and_sol_from_results_str(serial_str)
                serial_times_dict[i] = serial_time


            parareal_data = get_parareal_data_from_results_str(parareal_str, serial_sol)

            if dim == 1:
                section_info = get_section_info_from_str_1d(section_str)
            else:
                section_info = get_section_info_from_str_2d_polar(section_str)

            if x_type == "num ranks":
                x = section_info["num ranks"]
            else:
                x = section_info["theta"]

            if y_type == "time":
                y = get_time_to_meet_tolerance(parareal_data, tolerance=tolerance)
            else:
                y = get_iter_to_meet_tolerance(parareal_data, tolerance=tolerance)

            if y != -1:
                x_dict[i].append(x)
                y_dict[i].append(y)

    return x_dict, y_dict, serial_times_dict



