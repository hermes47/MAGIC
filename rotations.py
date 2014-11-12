import numpy as np

''' Returns the cross product matrix of a vector, ie the matrix:
     0  -z   y
     z   0  -x
    -y   x   0   '''
def crossProductMatrix(w):
    x = w[0]
    y = w[1]
    z = w[2]
    xProdMat = np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return xProdMat

''' Returns the normalised cross product of two vectors. '''
def normCrossProduct(a, b):
    w = np.cross(a, b)
    w = w/np.linalg.norm(w)
    return w

''' Returns the angle between two vectors. If the optional delta argument is supplied,
delta is added to each component of the first vector, and the angle returned is the angle
between that modified vector and the second one.'''
def vectAngle(vecA, vecB):
    theta = np.arccos(np.dot(vecA,vecB)/(np.linalg.norm(vecA)*np.linalg.norm(vecB)))
    return theta

''' Returns the rotation matrix that will rotate about a given axis, w, by a given
angle, theta.'''
def rotationMatrix(w, theta):
    ''' Given a unit vector, w, the rotation matrix is given as R = cos(t)I + sin(t)U + (1-cos(t))wXw
    where I is the identity matrix and U is the cross product matrix.'''
    compA = np.cos(theta)*np.identity(3, dtype=np.float64)
    compB = np.sin(theta)*crossProductMatrix(w)
    compC = (1-np.cos(theta))*np.outer(w,w)
    R = compA + compB + compC
    return R
        
''' Determines the longest distance between two atoms in a given molecule.
Returns the normalised vector representing the longest axis.'''       
def longAxis(mol): #may also need to return the atm indicies that define the vector, or just the indicies
    i = 0
    maxLength = 0.
    maxLengthVec = np.array([[0.,0.,0.]])
    while i < mol.getNumAtms()-1:
        j = i+1
        while j < mol.getNumAtms():
            vecIJ = mol.getAtms()[i].getXYZ(vec=True) - mol.getAtms()[j].getXYZ(vec=True)
            if np.linalg.norm(vecIJ) > maxLength:
                maxLength = np.linalg.norm(vecIJ)
                maxLengthVec = vecIJ
            j += 1
        i += 1
    return maxLengthVec/maxLength

if __name__ == "__main__":
    print("words")