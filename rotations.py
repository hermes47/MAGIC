import numpy as np


def crossProductMatrix(vectorInput):
    '''
    Returns the cross product matrix of a vector, ie the matrix:
     0  -z   y
     z   0  -x
    -y   x   0 
    '''
    if vectorInput.shape == (1,3):
        vectorInput = vectorInput.T
        x = vectorInput[0,0]
        y = vectorInput[1,0]
        z = vectorInput[2,0]
    elif vectorInput.shape == (3,):
        x = vectorInput[0]
        y = vectorInput[1]
        z = vectorInput[2]
    xProdMat = np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    return xProdMat

def rotationVectorDetermine(rotMolVec, refMolVec):
    # determine the rotation vector
    #w = np.cross(rotMolVec,refMolVec)
    w = np.cross(refMolVec, rotMolVec)
    # normalise the rotation vector
    w = w/np.linalg.norm(w)
    return w #np.array([1,1,1])/np.linalg.norm(np.array([1,1,1]))

def rotationAngleDetermine(rotMolVec, refMolVec, delta=False):
    if delta:
        # add same amount to each of x,y,z, such that dx^2 + dy^2 + dz^2 = delta
        for i in (0,1,2):
            if rotMolVec[i] > refMolVec[i]: rotMolVec[i] += np.sqrt(abs(delta/3.))
            elif rotMolVec[i] < refMolVec[i]: rotMolVec[i] -= np.sqrt(abs(delta/3.))
            
    # determine the new angle between the two
    newTheta = np.arccos(np.dot(rotMolVec,refMolVec)/(np.linalg.norm(rotMolVec)*np.linalg.norm(refMolVec)))
    # return the target theta
    return newTheta #1*np.pi/180

def rotateMolecule(rotMol, refMol, rotMolAtmIndex=-1, refMolAtmIndex=-1, displacement=0, vect=np.array([0,0,0]), ang=None, applyTo=None):
    rotMolVec = vect
    refMolVec = vect
    if rotMolAtmIndex is not -1:
        rotMolVec = rotMol.getAtmWithIndex(rotMolAtmIndex).getXYZ(vec=True)
    if refMolAtmIndex is not -1:
        refMolVec = refMol.getAtmWithIndex(refMolAtmIndex).getXYZ(vec=True)
    if np.linalg.norm(vect) == 0.:
        vect = rotationVectorDetermine(rotMolVec,refMolVec)
    if ang is None:
        ang = rotationAngleDetermine(rotMolVec, refMolVec, displacement)
    '''
    Given a unit vector, w, the rotation matrix is given as R = cos(t)I + sin(t)U + (1-cos(t))wXw
    where I is the identity matrix and U is the cross product matrix
    '''
    rotMat = np.cos(ang)*np.identity(3, dtype=np.float64) + np.sin(ang)*crossProductMatrix(vect) + (1-np.cos(ang))*np.outer(vect,vect)
    if applyTo is None:
        applyTo = []
        for num in range(rotMol.getNumAtms()): applyTo.append(rotMol.getAtms()[num].getAtmIndex())
    for atm in rotMol.getAtms():
        if atm.getAtmIndex() in applyTo:
            newCoordVec = np.dot(rotMat,atm.getXYZ(vec=True))
            atm.setXYZ(newCoordVec, vec=True)
        
def longAxis(mol):
    # determines the longest distance between two atoms in the molecule, and returns the vector between the two
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
    return maxLengthVec

def dihedralRotate(mol, dihedralIndex, bondIndex, atmList, numSteps, stepSize, printOut=False):
    print("stuff")

if __name__ == "__main__":
    print("words")