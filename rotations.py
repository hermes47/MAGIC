''' Contains a number of methods that are of use in rotational operations.
Methods do not act on the arguments they are supplied, and all return the
desired output.'''
import numpy as np

''' Returns the cross product matrix of a vector, ie the matrix:
     0  -z   y
     z   0  -x
    -y   x   0   '''
def crossProductMatrix(w):
    if w.shape == (1,3): # one of these shapes needs to go as it should be unused
        w = w.T
        x = w[0,0]
        y = w[1,0]
        z = w[2,0]
    elif w.shape == (3,):
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
def vectAngle(vecA, vecB, delta=None):
    if delta is not None:
        for i in (0,1,2):
            if vecA[i] > vecB[i]: vecA[i] += np.sqrt(abs(delta/3.))
            elif vecA[i] < vecB[i]: vecA[i] -= np.sqrt(abs(delta/3.))
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

def rotateMolecule(self, refMol, rotMolAtmIndex=-1, refMolAtmIndex=-1, displacement=0, vect=np.array([0,0,0]), ang=None, applyTo=None):
    rotMolVec = vect
    refMolVec = vect
    if rotMolAtmIndex is not -1:
        rotMolVec = self.getAtmWithIndex(rotMolAtmIndex).getXYZ(vec=True)
    if refMolAtmIndex is not -1:
        refMolVec = refMol.getAtmWithIndex(refMolAtmIndex).getXYZ(vec=True)
    if np.linalg.norm(vect) == 0.:
        vect = normCrossProduct(rotMolVec,refMolVec)
    if ang is None:
        ang = vectAngle(rotMolVec, refMolVec, displacement)
    rotMat = np.cos(ang)*np.identity(3, dtype=np.float64) + np.sin(ang)*crossProductMatrix(vect) + (1-np.cos(ang))*np.outer(vect,vect)
    if applyTo is None:
        applyTo = []
        for num in range(self.getNumAtms()): applyTo.append(self.getAtms()[num].getAtmIndex())
    for atm in self.getAtms():
        if atm.getAtmIndex() in applyTo:
            newCoordVec = np.dot(rotMat,atm.getXYZ(vec=True))
            atm.setXYZ(newCoordVec, vec=True)
            
        
''' Determines the longest distance between two atoms in a given molecule.
Returns the normalised vector representing the longest axis.'''       
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
    return maxLengthVec/np.linalg.norm(maxLengthVec)

def dihedralRotate(mol, dihedralIndex, bondIndex, atmList, numSteps, stepSize, printOut=False):
    print("stuff")

if __name__ == "__main__":
    print("words")