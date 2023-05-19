import numpy as np
from src.models.attitude.attitude import euler2exp_att, exp_att2euler, euler2quat
from src.models.attitude.attitude import *

def read_file(file, token_map):
    fh = open(file)
    lines = fh.readlines()

    data = []

    for line in lines:

        tokens = line.split(",")
        sow = float(tokens[token_map[0]])

        pos_x = float(tokens[token_map[1]])
        pos_y = float(tokens[token_map[2]])
        pos_z = float(tokens[token_map[3]])
        pos = np.array([pos_x, pos_y, pos_z])

        vel_x = float(tokens[token_map[4]])
        vel_y = float(tokens[token_map[5]])
        vel_z = float(tokens[token_map[6]])
        vel = np.array([vel_x, vel_y, vel_z])

        att_x = float(tokens[token_map[9]])  # roll
        att_y = float(tokens[token_map[8]])  # pitch
        att_z = float(tokens[token_map[7]])  # yaw
        euler = np.array([att_x, att_y, att_z])
        quat = euler2exp_att(pos, euler)
        print(euler)

        data.append(np.array([sow, *pos, *vel, *quat]))

    return np.array(data)
def test_exp_att():


    pvat_file = "C:\\NEXA_delivery\\tests\\audelpro\\test6_new_yaw\\input_ref\\reference_traj_smooth_yaw.csv"
    pvat_exp_att_file="C:\\NEXA_delivery\\tests\\audelpro\\test6_new_yaw\\input_ref\\reference_traj_smooth_yaw_expatt.csv"
    pvat_data = read_file(pvat_file, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    np.savetxt(pvat_exp_att_file, pvat_data, delimiter=",")

def test_get_quat():


    pvat_file = "C:\\NEXA_delivery\\tests\\audelpro\\test6_new_yaw\\input_ref\\reference_traj_smooth_yaw.csv"
    pvat_quat_file="C:\\NEXA_delivery\\tests\\audelpro\\test6_new_yaw\\input_ref\\reference_traj_smooth_yaw.csv_quat"
    pvat_data = read_file(pvat_file, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    np.savetxt(pvat_quat_file, pvat_data, delimiter=",")


def test_atts():
    #dcm = np.array([[1,2,3],[4,5,6],[7,8,9]])
    #quat = dcm2quat(dcm)

    #print(dcm)
    #print(quat)

    pos = np.array([4837583.0, -310732.0, 4132333.0])
    euler = np.array([-1.555,0.22,0.011]) # roll pitch yaw
    lat,lon,h = cartesian2geodetic(*pos)
    print(lat, lon, h)
    dcm_e_n = latlon2dcm_e_n(lat, lon)
    print("dcm_e_n")
    print(dcm_e_n)

    dcm_n_b = euler2dcm(euler)
    print("dcm_n_b")
    print(dcm_n_b)

    dcm_e_b = dcm_n_b @ dcm_e_n
    print("dcm_e_b")
    print(dcm_e_b)

    q = dcm2quat(dcm_e_b)
    print("final q")
    print(q)

    print("\n\n")
    att = euler2exp_att(pos, euler)
    euler2 = exp_att2euler(pos, att)
    print(euler)
    print(euler2)

def main():
    #test_exp_att()
    #test_get_quat()

    #test_to_ned()
    test_atts()
    #pos = np.array([4.837583308680000e+06 , -3.107328416690000e+05  , 4.132333731440000e+06])
    #euler = np.array([0,0,2.777584889938198480e-04])
    #exp_att = euler2exp_att(pos,euler)
    #print(exp_att)

print("#------------------------------------------------#")
print("#           Welcome to PyGNSSFix TESTS           #")
print("#------------------------------------------------#\n")

if __name__ == "__main__":
    main()
