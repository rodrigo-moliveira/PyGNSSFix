from src.models.frames.attitude import *


def quat_testing():
    x_A = np.array([1023.23, 904.23, -50.323]) # x vector in A
    euler_angles = np.array([0.34, 2.3, 1.7]) # from A to B
    dcm_A_B = euler2dcm(euler_angles) # DCM from A to B
    x_B_dcm = dcm_A_B @ x_A # x_B vector rotated with DCMs
    quat = dcm2quat(dcm_A_B)
    x_B_quat = quat2dcm(quat) @ x_A
    print("-----------------------------------")
    print("Test 1: DCM and quat2dcm conversion")
    print("dcm_A_B = ", dcm_A_B)
    print("quat = ", quat)
    print("x_A = ", x_A)
    print("x_B_dcm = ", x_B_dcm)
    print("x_B_quat = ", x_B_quat)

    print("\n\n")
    print("-----------------------------------------")
    print("Test 2: Successive DCM and quat rotations")
    euler_A = np.array([0.34, 2.3, 1.7])  # from A to B
    euler_B = np.array([-0.83, 1.2948, -0.48])  # from B to C
    dcm_A_B = euler2dcm(euler_A)
    dcm_B_C = euler2dcm(euler_B)
    dcm_A_C = dcm_B_C @ dcm_A_B # DCM_A^C = DCM_B^C * DCM_A^B
    quat_A_B = dcm2quat(dcm_A_B)
    quat_B_C = dcm2quat(dcm_B_C)
    quat_A_C = quat_multiply(quat_B_C, quat_A_B) # quat_A^C = quat_B^C * quat_A^B
    dcm_A_C_from_quat = quat2dcm(quat_A_C)
    x_C_dcm = dcm_A_C @ x_A
    x_C_quat = dcm_A_C_from_quat @ x_A
    print("dcm_A_C = ", dcm_A_C)
    print("dcm_A_C_from_quat = ", dcm_A_C_from_quat)
    print("x_C_dcm = ", x_C_dcm)
    print("x_C_quat = ", x_C_quat)

    print("\n\n")
    print("-------------------------")
    print("Test 3: Rotating a vector")
    euler_A = np.array([-1.2323, 0.234, 1.7])  # from A to B
    dcm_A = euler2dcm(euler_A)
    quat_A = dcm2quat(dcm_A)
    x_B_dcm = dcm_A @ x_A
    x_B_quat = quat2dcm(quat_A) @ x_A
    x_A_quat = np.array([x_A[0], x_A[1], x_A[2], 0])
    x_B_quat_2 = quat_multiply(quat_multiply(quat_A, x_A_quat), quat_conjugate(quat_A)) # x_B = q * x_A * q_conj
    x_B_quat_3 = quat_rotate(quat_A, x_A)
    print("x_B_dcm = ", x_B_dcm)
    print("x_B_quat = ", x_B_quat)
    print("x_B_quat_2 = ", x_B_quat_2[0:3])
    print("x_B_quat_3 = ", x_B_quat_3)


if __name__ == "__main__":
    quat_testing()
