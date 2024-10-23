#!/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from src.models.frames import cartesian2geodetic
from src.models.frames.frames import ecef2ned


#script to compare PVAT files from different sources (for example, compare INS integration from audelpro and NEXA)

def read_file(file, is_csv, ignore_header, token_map, time_millis=False):
    fh = open(file)
    lines = fh.readlines()
    count = 0

    data = {}

    for line in lines:
        if ignore_header and count == 0:
            count = 1
            continue

        if is_csv:
            tokens = line.split(",")
        else:
            tokens = line.split()
        # print(tokens)
        if time_millis:
            sow = float(tokens[token_map[0]]) / 1000.0
        else:
            sow = float(tokens[token_map[0]])

        pos_x = float(tokens[token_map[1]])
        pos_y = float(tokens[token_map[2]])
        pos_z = float(tokens[token_map[3]])
        pos = np.array([pos_x, pos_y, pos_z])

        vel_x = float(tokens[token_map[4]])
        vel_y = float(tokens[token_map[5]])
        vel_z = float(tokens[token_map[6]])
        vel = np.array([vel_x, vel_y, vel_z])

        att_x = float(tokens[token_map[7]])
        att_y = float(tokens[token_map[8]])
        att_z = float(tokens[token_map[9]])
        att = np.array([att_x, att_y, att_z])

        data[sow] = [pos, vel, att]

    return data

def plot_error(error, title, y_label):
    # compute rms
    rms = {}
    for k,v in error.items():
        rms[k] = np.linalg.norm(v)

    x = [i[0] for i in error.values()]
    y = [i[1] for i in error.values()]
    z = [i[2] for i in error.values()]

    fig, ax = plt.subplots()
    ax.set_title(f"RMS {title} error")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(y_label)
    ax.plot(rms.keys(), rms.values())

    fig, ax = plt.subplots()
    ax.set_title(f"{title} error")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(y_label)
    ax.plot(rms.keys(), x, label="north")
    ax.plot(rms.keys(), y, label="east")
    ax.plot(rms.keys(), z, label="down")
    ax.legend()


def main():

    ######### TRAJETORIA 8 #############
    # NEXA TRAJ
    input_file1 = "C:\\NEXA_delivery\\tests\\audelpro\\test6_new_yaw\\input_ref\\reference_traj_smooth_yaw.csv"
    data1 = read_file(input_file1, True, False, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], time_millis=True)
    input_file2 = "C:\\NEXA_delivery\\tests\\audelpro\\test6_new_yaw\\nexa\\04_out_NEXAto\\IMU_trajectory___fw__.ecef.txt"
    data2 = read_file(input_file2, False, True, [0, 1, 2, 3, 4, 5, 6, 7,8,9])


    error_pos = {}
    error_vel = {}

    for epoch in data1.keys():
        if epoch in data1 and epoch in data2:
            d1 = data1[epoch]  # ref data
            d2 = data2[epoch]  # est data

            ref_pos = d1[0]  # ref pos
            ref_vel = d1[1]  #

            est_pos = d2[0]
            est_vel = d2[1]

            llh = cartesian2geodetic(*ref_pos)
            error_pos[epoch] = ecef2ned(est_pos, ref_pos, llh)
            error_vel[epoch] = ecef2ned(est_vel, ref_vel, llh)

    plot_error(error_pos, "position", "error [m]")
    plot_error(error_vel, "velocity", "error [m/s]")

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    # ax.plot(points[0], points[1], points[2], marker='x')
    # ax.scatter(*points.T[0], color='red')

    #traj = data1.values()
    #x = [i[0][0] for i in traj]
    #y = [i[0][1] for i in traj]
    #z = [i[0][2] for i in traj]

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot(x, y, z)
    plt.show()

main()