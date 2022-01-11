import numpy as np
from scipy.spatial.transform.rotation import Rotation as R


def eul2rotm(e):
    return R.from_euler('zyx', e).as_matrix()


def rotm2eul(r):
    return np.matrix(R.from_matrix(r.A).as_euler('zyx'))


def p2pMotionPlan(dk, ik, xEs, phiEs, xEf, phiEf, minT, maxT, dt):
    qEs = ik(xEs, eul2rotm(phiEs.transpose()))
    qEf = ik(xEf, eul2rotm(phiEf.transpose()))
    qEs = qEs[0, :].A1
    qEf = qEf[0, :].A1

    M = np.matrix([[1, minT, minT ** 2, minT ** 3],
                   [0, 1, 2 * minT, 3 * minT ** 2],
                   [1, maxT, maxT ** 2, maxT ** 3],
                   [0, 1, 2 * maxT, 3 * maxT ** 2]])

    A = []
    for i in range(qEs.size):
        b = np.matrix([[qEs[i], 0, qEf[i], 0]]).T
        a = (M.I * b).T.A1
        A.append(a)

    A = np.matrix(A)
    Th = []
    xE = []
    phiE = []
    for t in np.arange(minT, maxT, dt):
        th = [t]
        for i in range(qEs.size):
            q = A[i, 0] + A[i, 1] * t + A[i, 2] * t ** 2 + A[i, 3] * t ** 3
            th.append(q)
        Th.append(th)

        mx, mR = dk(th[1:7])
        t_mx_row = np.concatenate((np.matrix([[t]]), mx.T), axis=1)
        xE.append(t_mx_row.A1.tolist())
        t_rot_row = np.concatenate((np.matrix([[t]]), rotm2eul(mR)), axis=1)
        phiE.append(t_rot_row.A1.tolist())

    return Th, xE, phiE
