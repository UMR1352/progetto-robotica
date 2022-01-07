import numpy as np
from scipy.spatial.transform.rotation import Rotation as R

def eul2rotm(e):
    return R.from_euler('zyx', e)

def rotm2eul(r):
    return R.as_euler('zyx')

def p2pMotionPlan(dk, ik, xEs, phiEs, xEf, phiEf, minT, maxT, dt):
    qEs = ik(xEs, eul2rotm(phiEs.transpose()))
    qEf = ik(xEf, eul2rotm(phiEf.transpose()))
    qEs = qEs[0, :]
    qEf = qEf[0, :]

    M = np.matrix([[1, minT, minT ** 2, minT ** 3],
                   [0, 1, 2 * minT, 3 * minT ** 2],
                   [1, maxT, maxT ** 2, maxT ** 3],
                   [0, 1, 2 * maxT, 3 * maxT ** 2]])

    A = np.matrix()
    for i in range(qEs[0, :].length()):
        b = np.array([qEs[i], 0, qEf(i), 0])
        a = M.inverse() * b
        A = np.concatenate((A, a.transpose()), axis=0) # Add a' as row
    
    Th = np.matrix()
    xE = np.matrix()
    phiE = np.matrix()
    for t in range(minT, maxT, dt):
        th = np.array(t)
        for i in range(qEs.length()):
            q = A[i, 0] + A[i, 1] * t + A[i, 2] * t ** 2 + A[i, 3] * t ** 3
            th = np.concatenate((th, q), axis=1) # Add q as col
        Th = np.concatenate((Th, th), axis=0) # Add th as row
        [mx, mR] = dk(th[1:6]) # This won't work since python can't pattern match or destructure... fuck this shit
        t_mx_row = np.concatenate((np.array([t]), mx.transpose()), axis=1)
        xE = np.concatenate((xE, t_mx_row), axis=0)
        t_rot_row = np.concatenate((np.array([t]), rotm2eul(mR)), axis=1)
        phiE = np.concatenate((phiE, t_rot_row), axis=0)
    
    return Th, xE, phiE
    