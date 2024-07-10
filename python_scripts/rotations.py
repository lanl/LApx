import numpy as np

def rotationX(alpha):
    Rx =np.zeros([3,3])
    c = np.cos(np.deg2rad(alpha))
    s = np.sin(np.deg2rad(alpha))
    Rx[0,0] = 1.
    Rx[1,1] = c
    Rx[1,2] = -s
    Rx[2,1] = s
    Rx[2,2] = c
    return Rx

def rotationY(alpha):
    Ry =np.zeros([3,3])
    c = np.cos(np.deg2rad(alpha))
    s = np.sin(np.deg2rad(alpha))

    Ry[0,0] = c
    Ry[0,2] = s
    Ry[1,1] = 1.
    Ry[2,0] = -s
    Ry[2,2] = c
    return Ry

def rotationZ(alpha):
    Rz =np.zeros([3,3])
    c = np.cos(np.deg2rad(alpha))
    s = np.sin(np.deg2rad(alpha))

    Rz[0,0] = c
    Rz[0,1] = -s
    Rz[1,0] = s
    Rz[1,1] = c
    Rz[2,2] = 1.
    return Rz

def bungeRotationMatrixSample2Crystal(ph, th, tm):
    r_mtx = np.zeros([3,3])
    cph = np.cos(ph)
    sph = np.sin(ph)
    cth = np.cos(th)
    sth = np.sin(th)
    ctm = np.cos(tm)
    stm = np.sin(tm)

    r_mtx[0,0]=ctm*cph-sph*stm*cth
    r_mtx[1,0]=-stm*cph-sph*ctm*cth
    r_mtx[2,0]=sph*sth
    r_mtx[0,1]=ctm*sph+cph*stm*cth
    r_mtx[1,1]=-sph*stm+cph*ctm*cth
    r_mtx[2,1]=-sth*cph
    r_mtx[0,2]=sth*stm
    r_mtx[1,2]=ctm*sth
    r_mtx[2,2]=cth

    return r_mtx

def bungeRotationMatrixCrystal2Sample(ph, th, tm):
    return np.transpose(bungeRotationMatrixSample2Crystal(ph, th, tm))

def bungeAnlgesFromRotationMatrixSample2Crystal(r_mtx) :
    tol = 1e-6

    th=np.arccos(r_mtx[2,2])
    if (abs(r_mtx[2,2]) > (1.-tol) ):
       tm=0.
       ph=np.atan2(r_mtx[0,1],r_mtx[0,0])
    else :
       sth=np.sin(th)
       tm=np.arctan2(r_mtx[0,2]/sth,r_mtx[1,2]/sth)
       ph=np.arctan2(r_mtx[2,0]/sth,-r_mtx[0,1]/sth)

    return np.asarray([ph, th, tm])

def bungeAnlgesFromRotationMatrixCrystal2Sample(r_mtx) :
    return bungeAnlgesFromRotationMatrixSample2Crystal(r_mtx.T)

def rotateVectorToZ(vec):
    R=np.zeros([3,3])
    v1 = np.zeros([3])
    # normalize vector
    vec /= np.linalg.norm(vec)
    w = np.abs(vec)
    if ((w[2] >= w[1] and w[1] >= w[0]) or (w[1] >= w[2] and w[2] >= w[0])):
        # vec(0) is the smallest component
        v1[0] = 1
    elif ((w[2] >= w[0] and w[0] >= w[1]) and ( w[0] >= w[2] and w[2] >= w[1]) ):
        # vec(1) is the smallest component
        v1[1] = 1
    else :
        # vec(2) is the smallest component
        v1[2] = 1;

    v1 -= np.dot(v1, vec)*vec
    v1 /= np.linalg.norm(v1)

    v0 = np.cross(v1,vec)

    for i in range(3):
        R[0,i] = v0[i]
        R[1,i] = v1[i]
        R[2,i] = vec[i]
    return R

def rotVec1ToVec2(vec1, vec2):
  rot1_to_z = rotateVectorToZ(vec1)
  rot2_to_z = rotateVectorToZ(vec2)
  return np.dot(rot2_to_z.T, rot1_to_z)


if __name__ == "__main__":
    b = np.zeros([3])
    n = np.zeros([3])
    b[0] = 0
    b[1] = 1
    b[2] = -1
    n[0] = 1
    n[1] = 1
    n[2] = 1

    print("b", b)
    print("n", n)
    R1_c2S = rotVec1ToVec2(b, np.asarray([0.,1.,0.]))
    print("R1_c2S", R1_c2S)
    b1 = np.dot(R1_c2S, b)
    n1 = np.dot(R1_c2S, n)
    print("b1", b1)
    print("n1", n1)
    R2_c2S = rotVec1ToVec2(n1, np.asarray([0.,0.,1.]))
    b2 = np.dot(R2_c2S, b1)
    n2 = np.dot(R2_c2S, n1)
    print("b2", b2)
    print("n2", n2)

    R_c2s = np.dot(R2_c2S,R1_c2S)
    print("R*b",np.dot(R_c2s,b/np.linalg.norm(b)))
    print("R*n",np.dot(R_c2s,n/np.linalg.norm(n)))
    angles = bungeAnlgesFromRotationMatrixCrystal2Sample(R_c2s)
    print(np.rad2deg(angles))
    print( np.deg2rad((np.rad2deg(angles)-90.))/(np.pi/2.))


    # rotatate 5 degrees around the y Axis
    R3 = rotationZ(1.)
    R_c2s = np.dot(R3,R_c2s)
    print("R*b",np.dot(R_c2s,b/np.linalg.norm(b)))
    print("R*n",np.dot(R_c2s,n/np.linalg.norm(n)))
    angles = bungeAnlgesFromRotationMatrixCrystal2Sample(R_c2s)
    print(np.rad2deg(angles))
    print( np.deg2rad((np.rad2deg(angles)-90.))/(np.pi/2.))

    # print(bungeRotationMatrixCrystal2Sample(angles[0],angles[1], angles[2]))
