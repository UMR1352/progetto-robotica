import numpy as np
from numpy import sin, cos, pi, arctan2, arccos, arcsin, hypot
from scipy.linalg import inv

A = [0, -0.425, -0.3922, 0, 0, 0]
D = [0.1625, 0, 0, 0.1333, 0.0997, 0.0996]


def ur5Inverse(p60, R60):
    T60 = np.asmatrix(np.vstack(
        (np.column_stack((R60, p60)), np.array([0, 0, 0, 1]))))
    # Finding th1
    p50 = (T60 * np.matrix([0, 0, -D[5], 1]).T).A1
    th1_1 = arctan2(p50[1], p50[0]) + arccos(D[3] /
                                             hypot(p50[1], p50[0])) + pi / 2
    th1_2 = arctan2(p50[1], p50[0]) - arccos(D[3] /
                                             hypot(p50[1], p50[0])) + pi / 2

    # Finding th5
    th5_1 = arccos((p60[0] * sin(th1_1) - p60[1] * cos(th1_1) - D[3]) / D[5])
    th5_2 = -arccos((p60[0] * sin(th1_1) - p60[1] * cos(th1_1) - D[3]) / D[5])
    th5_3 = arccos((p60[0] * sin(th1_2) - p60[1] * cos(th1_2) - D[3]) / D[5])
    th5_4 = -arccos((p60[0] * sin(th1_2) - p60[1] * cos(th1_2) - D[3]) / D[5])

    # Just a bit more of magic stuff
    T06 = inv(T60)
    Xhat = T06[0:3, 0]
    Yhat = T06[0:3, 1]

    th6_1 = arctan2((-Xhat[1] * sin(th1_1) + Yhat[1] * cos(th1_1)) / sin(th5_1),
                    (Xhat[0] * sin(th1_1) - Yhat[0] * cos(th1_1)) / sin(th5_1))
    th6_2 = arctan2((-Xhat[1] * sin(th1_1) + Yhat[1] * cos(th1_1)) / sin(th5_2),
                    (Xhat[0] * sin(th1_1) - Yhat[0] * cos(th1_1)) / sin(th5_2))
    th6_3 = arctan2((-Xhat[1] * sin(th1_1) + Yhat[1] * cos(th1_2)) / sin(th5_3),
                    (Xhat[0] * sin(th1_2) - Yhat[0] * cos(th1_2)) / sin(th5_3))
    th6_4 = arctan2((-Xhat[1] * sin(th1_1) + Yhat[1] * cos(th1_2)) / sin(th5_4),
                    (Xhat[0] * sin(th1_2) - Yhat[0] * cos(th1_2)) / sin(th5_4))

    T41m = inv(T10f(th1_1)) * T60 * inv(T65f(th6_1)) * inv(T54f(th5_1))
    p41_1 = np.asarray(T41m[0:3, 3]).flatten()
    p41xz_1 = hypot(p41_1[0], p41_1[2])
    print(p41xz_1)

    T41m = inv(T10f(th1_1)) * T60 * inv(T65f(th6_2)) * inv(T54f(th5_2))
    p41_2 = np.asarray(T41m[0:3, 3]).flatten()
    p41xz_2 = hypot(p41_2[0], p41_2[2])

    T41m = inv(T10f(th1_2)) * T60 * inv(T65f(th6_3)) * inv(T54f(th5_3))
    p41_3 = np.asarray(T41m[0:3, 3]).flatten()
    p41xz_3 = hypot(p41_3[0], p41_3[2])

    T41m = inv(T10f(th1_2)) * T60 * inv(T65f(th6_4)) * inv(T54f(th5_4))
    p41_4 = np.asarray(T41m[0:3, 3]).flatten()
    p41xz_4 = hypot(p41_4[0], p41_4[2])

    # There are for real 8 fucking possible values for th3
    th3_1 = arccos((p41xz_1 ** 2 - A[1] ** 2 - A[2] ** 2) / (2 * A[1] * A[2]))
    th3_2 = arccos((p41xz_2 ** 2 - A[1] ** 2 - A[2] ** 2) / (2 * A[1] * A[2]))
    th3_3 = arccos((p41xz_3 ** 2 - A[1] ** 2 - A[2] ** 2) / (2 * A[1] * A[2]))
    th3_4 = arccos((p41xz_4 ** 2 - A[1] ** 2 - A[2] ** 2) / (2 * A[1] * A[2]))
    th3_5 = -th3_1
    th3_6 = -th3_2
    th3_7 = -th3_3
    th3_8 = -th3_4

    # Same shit for th2
    th2_1 = arctan2(-p41_1[2], -p41_1[0]) - \
        arcsin((-A[2] * sin(th3_1)) / p41xz_1)
    th2_2 = arctan2(-p41_2[2], -p41_2[0]) - \
        arcsin((-A[2] * sin(th3_2)) / p41xz_2)
    th2_3 = arctan2(-p41_3[2], -p41_3[0]) - \
        arcsin((-A[2] * sin(th3_3)) / p41xz_3)
    th2_4 = arctan2(-p41_4[2], -p41_4[0]) - \
        arcsin((-A[2] * sin(th3_4)) / p41xz_4)
    th2_5 = arctan2(-p41_1[2], -p41_1[0]) - \
        arcsin((A[2] * sin(th3_1)) / p41xz_1)
    th2_6 = arctan2(-p41_2[2], -p41_2[0]) - \
        arcsin((A[2] * sin(th3_2)) / p41xz_2)
    th2_7 = arctan2(-p41_3[2], -p41_3[0]) - \
        arcsin((A[2] * sin(th3_3)) / p41xz_3)
    th2_8 = arctan2(-p41_4[2], -p41_4[0]) - \
        arcsin((A[2] * sin(th3_4)) / p41xz_4)

    T43m = T32f(th3_1).I * T21f(th2_1).I * T10f(th1_1).I * \
        T60 * T65f(th6_1).I * T54f(th5_1).I
    Xhat43 = T43m[0:2, 0].A1
    th4_1 = arctan2(Xhat43[1], Xhat43[0])

    T43m = T32f(th3_2).I * T21f(th2_2).I * T10f(th1_1).I * \
        T60 * T65f(th6_2).I * T54f(th5_2).I
    Xhat43 = T43m[0:2, 0].A1
    th4_2 = arctan2(Xhat43[1], Xhat43[0])

    T43m = T32f(th3_3).I * T21f(th2_3).I * T10f(th1_2).I * \
        T60 * T65f(th6_3).I * T54f(th5_3).I
    Xhat43 = T43m[0:2, 0].A1
    th4_3 = arctan2(Xhat43[1], Xhat43[0])

    T43m = T32f(th3_4).I * T21f(th2_4).I * T10f(th1_2).I * \
        T60 * T65f(th6_4).I * T54f(th5_4).I
    Xhat43 = T43m[0:2, 0].A1
    th4_4 = arctan2(Xhat43[1], Xhat43[0])

    T43m = T32f(th3_5).I * T21f(th2_5).I * T10f(th1_1).I * \
        T60 * T65f(th6_1).I * T54f(th5_1).I
    Xhat43 = T43m[0:2, 0].A1
    th4_5 = arctan2(Xhat43[1], Xhat43[0])

    T43m = T32f(th3_6).I * T21f(th2_6).I * T10f(th1_1).I * \
        T60 * T65f(th6_2).I * T54f(th5_2).I
    Xhat43 = T43m[0:2, 0].A1
    th4_6 = arctan2(Xhat43[1], Xhat43[0])

    T43m = T32f(th3_7).I * T21f(th2_7).I * T10f(th1_2).I * \
        T60 * T65f(th6_3).I * T54f(th5_3).I
    Xhat43 = T43m[0:2, 0].A1
    th4_7 = arctan2(Xhat43[1], Xhat43[0])

    T43m = T32f(th3_8).I * T21f(th2_8).I * T10f(th1_2).I * \
        T60 * T65f(th6_4).I * T54f(th5_4).I
    Xhat43 = T43m[0:2, 0].A1
    th4_8 = arctan2(Xhat43[1], Xhat43[0])

    return np.matrix([[th1_1, th2_1, th3_1, th4_1, th5_1, th6_1],
                      [th1_1, th2_2, th3_2, th4_2, th5_2, th6_2],
                      [th1_2, th2_3, th3_3, th4_3, th5_3, th6_3],
                      [th1_2, th2_4, th3_4, th4_4, th5_4, th6_4],
                      [th1_1, th2_5, th3_5, th4_5, th5_1, th6_1],
                      [th1_1, th2_6, th3_6, th4_6, th5_2, th6_2],
                      [th1_2, th2_7, th3_7, th4_7, th5_4, th6_4],
                      [th1_2, th2_8, th3_8, th4_8, th5_4, th6_4]])


def T10f(th):
    return np.matrix([[cos(th), -sin(th), 0, 0],
                      [sin(th), cos(th), 0, 0],
                      [0, 0, 1, D[0]],
                      [0, 0, 0, 1]])


def T21f(th):
    return np.matrix([[cos(th), -sin(th), 0, 0],
                      [0, 0, -1, 0],
                      [sin(th), cos(th), 0, 0],
                      [0, 0, 0, 1]])


def T32f(th):
    return np.matrix([[cos(th), -sin(th), 0, A[1]],
                      [sin(th), cos(th), 0, 0],
                      [0, 0, 1, D[2]],
                      [0, 0, 0, 1]])


def T43f(th):
    return np.matrix([[cos(th), -sin(th), 0, A[2]],
                      [sin(th), cos(th), 0, 0],
                      [0, 0, 1, D[3]],
                      [0, 0, 0, 1]])


def T54f(th):
    return np.matrix([[cos(th), -sin(th), 0, 0],
                      [0, 0, -1, -D[4]],
                      [sin(th), cos(th), 0, 0],
                      [0, 0, 0, 1]])


def T65f(th):
    return np.matrix([[cos(th), -sin(th), 0, 0],
                      [0, 0, 1, D[5]],
                      [-sin(th), -cos(th), 0, 0],
                      [0, 0, 0, 1]])
