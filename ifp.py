
class IFP ():
    # parse the given IFP file and store it in memory
    def __init__(self, ifpData):
        self.title = []
        self.forcefield = ""
        self.massAtomType = []
        self.bondStretchType = []
        self.bondAngleBendType = []
        self.impDihedralType = []
        self.torsDihedralType = []
        self.singleAtomLJPair = []
        self.mixedAtomLJPair = []
        self.specAtomLJPair = []
        print("IFP data loading not yet implemented")
    def __repr__(self, *args, **kwargs):
        return ""
    ''' Methods to return the data stored in the class.'''
    def getTitle(self):
        return self.title
    def getForcefield(self):
        return self.forcefield
    def getMassAtomTypes(self):
        return self.massAtomType
    def getBondStretchTypes(self):
        return self.bondStretchType
    def getBondAngleBendTypes(self):
        return self.bondAngleBendType
    def getImpDihedralTypes(self):
        return self.impDihedralType
    def getTorsDihedralTypes(self):
        return self.torsDihedralType
    def getSingleAtomLJPair(self):
        return self.singleAtomLJPair
    def getMixedAtomLJPair(self):
        return self.mixedAtomLJPair
    def getSpecAtomLJPair(self):
        return self.specAtomLJPair
    ''' Methods to set the data in the class that are not lists.'''
    def setForcefield(self, ff):
        self.forcefield = ff
    ''' Methods to add to the various lists. '''
    def addMassAtomType(self, mat):
        # mat is a dictionary with the keys: N ATMAS and ATMASN
        # N is the atom mass code
        # ATMAS is the mass of the atom mass code
        # ATMASN is the name of the atom mass code
        self.massAtomType.append(mat)
    def addBondStretchType(self, bst):
        # bst is a dictionary with the keys: N CB HB B0
        # N is the bond stretch type code
        # CB is the quartic bond-stretch force constant
        # HB is the harmonic bond stretch force constant
        # B0 is the ideal bond length
        self.bondStretchType.append(bst)
    def addBondAngleBendType(self, bab):
        # bab is a dictionary with the keys: N CT CHT and T0
        # N is the angle bend type code
        # CT is the non-harmonic force constant
        # CHT is the harmonic force constant
        # T0 is the ideal bond angle
        self.bondAngleBendType.append(bab)
    def addImpDihedralType(self, idt):
        # idt is a dictionary with the keys: N CQ Q0
        # N is the improper dihedral-angle type code
        # CQ is the force constant
        # Q0 is the ideal improper angle
        self.impDihedralType.append(idt)
    def addTorsDihedralType(self, tdt):
        # tdt is a dictionary with the keys: N CP PD and MP
        # N is the dihedral-angle type code
        # CP is the force constant
        # PD is the phase shift
        # MP is the multiplicity
        self.torsDihedralType.append(tdt)
    def addSingleAtomLJPair(self, saljp):
        # saljp is a dictionary with the keys:  N TYPE C6 C12_1 C12_2 C12_3 CS6 CS12 INT
        # N is the atom type code
        # TYPE is the atom type name
        # C6 is the c6 parameter
        # C12_1-3 are the possible C12 parameters
        # CS6 is the c6 parameter for LJ14Pair
        # CS12 is the c12 parameter for LJ14Pair
        # Int is a list of length (num atm types) assigning each paired interaction to a C12 value
        self.singleAtomLJPair.append(saljp)
    #def addMixedAtomLJPair(self, maljp):
        # maljp is a dictionary with the keys: 
if __name__ == "__main__":
    ifp=IFP("stuff")