import numpy as np
from numpy import cos, sin

A = [0, -0.425, -0.3922, 0, 0, 0]
D = [0.1625, 0, 0, 0.1333, 0.0997, 0.0996]


def ur5Direct(Th):
    T10f = np.matrix([[cos(Th[0]), -sin(Th[0]), 0, 0],
                      [sin(Th[0]), cos(Th[0]), 0, 0],
                      [0, 0, 1, D[0]],
                      [0, 0, 0, 1]])
    T21f = np.matrix([[cos(Th[1]), -sin(Th[1]), 0, 0],
                      [0, 0, -1, 0],
                      [sin(Th[1]), cos(Th[1]), 0, 0],
                      [0, 0, 0, 1]])
    T32f = np.matrix([[cos(Th[2]), -sin(Th[2]), 0, A[1]],
                      [sin(Th[2]), cos(Th[2]), 0, 0],
                      [0, 0, 1, D[2]],
                      [0, 0, 0, 1]])
    T43f = np.matrix([[cos(Th[3]), -sin(Th[3]), 0, A[2]],
                      [sin(Th[3]), cos(Th[3]), 0, 0],
                      [0, 0, 1, D[3]],
                      [0, 0, 0, 1]])
    T54f = np.matrix([[cos(Th[4]), -sin(Th[4]), 0, 0],
                      [0, 0, -1, -D[4]],
                      [sin(Th[4]), cos(Th[4]), 0, 0],
                      [0, 0, 0, 1]])
    T65f = np.matrix([[cos(Th[5]), -sin(Th[5]), 0, 0],
                      [0, 0, 1, D[5]],
                      [-sin(Th[5]), -cos(Th[5]), 0, 0],
                      [0, 0, 0, 1]])
    T06 = T10f * T21f * T32f * T43f * T54f * T65f

    return T06[0:3, 3], T06[0:3, 0:3]
