import io

class IFP ():
    # parse the given IFP file and store it in memory
    def __init__(self):
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
        self.maxMassType = 0
        self.maxStretchType = 0
        self.maxBendType = 0
        self.maxImpType = 0
        self.maxTorsType = 0
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
    def addTitleLine(self, line):
        self.title.append(line)
    def addMassAtomType(self, mat):
        # mat is a dictionary with the keys: N ATMAS and ATMASN
        # N is the atom mass code
        # ATMAS is the mass of the atom mass code
        # ATMASN is the name of the atom mass code
        self.massAtomType.append(mat)
        if mat['N'] > self.MaxMassType(): self.maxMassType = mat['N']
    def addBondStretchType(self, bst):
        # bst is a dictionary with the keys: N CB HB B0
        # N is the bond stretch type code
        # CB is the quartic bond-stretch force constant
        # HB is the harmonic bond stretch force constant
        # B0 is the ideal bond length
        self.bondStretchType.append(bst)
        if bst['N'] > self.MaxStretchType(): self.maxStretchType = bst['N']
    def addBondAngleBendType(self, bab):
        # bab is a dictionary with the keys: N CT CHT and T0
        # N is the angle bend type code
        # CT is the non-harmonic force constant
        # CHT is the harmonic force constant
        # T0 is the ideal bond angle
        self.bondAngleBendType.append(bab)
        if bab['N'] > self.MaxBendType(): self.maxBendType = bab['N']
    def addImpDihedralType(self, idt):
        # idt is a dictionary with the keys: N CQ Q0
        # N is the improper dihedral-angle type code
        # CQ is the force constant
        # Q0 is the ideal improper angle
        self.impDihedralType.append(idt)
        if idt['N'] > self.MaxImpType(): self.maxImpType = idt['N']
    def addTorsDihedralType(self, tdt):
        # tdt is a dictionary with the keys: N CP PD and MP
        # N is the dihedral-angle type code
        # CP is the force constant
        # PD is the phase shift
        # MP is the multiplicity
        self.torsDihedralType.append(tdt)
        if tdt['N'] > self.MaxTorsType(): self.maxTorsType = tdt['N']
    def addSingleAtomLJPair(self, dt):
        self.singleAtomLJPair.append(dt)
    def addMixedAtomLJPair(self, dt):
        self.mixedAtomLJPair.append(dt)
    def addSpecAtomLJPair(self, dt):
        self.specAtomLJPair.append(dt)
    def NumMassTypes(self):
        return len(self.getMassAtomTypes())
    def NumStretchTypes(self):
        return len(self.getBondStretchTypes())
    def NumBendTypes(self):
        return len(self.getBondAngleBendTypes())
    def NumImpTypes(self):
        return len(self.getImpDihedralTypes())
    def NumTorsTypes(self):
        return len(self.getTorsDihedralTypes())
    def MaxMassType(self):
        return self.maxMassType
    def MaxStretchType(self):
        return self.maxStretchType
    def MaxBendType(self):
        return self.maxBendType
    def MaxImpType(self):
        return self.maxImpType
    def MaxTorsType(self):
        return self.maxTorsType
    def getTorsWithMultiplicity(self, n):
        torreturn = []
        for dihed in self.getTorsDihedralTypes():
            if dihed['MP'] == n: torreturn.append(dihed)
        return torreturn
    def printIFP(self):
        data = {}
        data['titleLines'] = '\n'.join(x for x in self.getTitle())
        data['ff'] = self.getForcefield()
        data['masscount'] = self.NumMassTypes()
        data['maxmass'] =self.MaxMassType()
        data['stretchcount'] = self.NumStretchTypes()
        data['maxstretch'] =self.MaxStretchType()
        data['bendcount'] = self.NumBendTypes()
        data['maxbend'] =self.MaxBendType()
        data['impcount'] = self.NumImpTypes()
        data['maximp'] =self.MaxImpType()
        data['torscount'] = self.NumTorsTypes()
        data['maxtors'] =self.MaxTorsType()
        data['massatomtypelines'] = '\n'.join('{N}     {ATMAS}   {ATMASN}'.format(**x) for x in self.getMassAtomTypes())
        data['bondstretchtypelines'] = '\n'.join('{N}  {CB}  {HB}  {B0}'.format(**x) for x in self.getBondStretchTypes())
        data['bondanglebendtypelines'] = '\n'.join('{N}  {CT}  {CHT}  {T0}'.format(**x) for x in self.getBondAngleBendTypes())
        data['impdihedraltypelines'] = '\n'.join('{N}  {CQ}  {Q0}'.format(**x) for x in self.getImpDihedralTypes())
        data['torsdihedraltypelines'] = '\n'.join('{N}  {CP}  {PD}  {MP}'.format(**x) for x in self.getTorsDihedralTypes())
        data['singleatomLJ'] = '\n'.join(x for x in self.getSingleAtomLJPair())
        data['mixedatomLJ'] = '\n'.join(x for x in self.getMixedAtomLJPair())
        data['specatomLJ'] = '\n'.join(c for c in self.getSpecAtomLJPair())
        f = io.StringIO()
        print('''TITLE
{titleLines}
END
FORCEFIELD
{ff}
END
MASSATOMTYPECODE
  {masscount}    {maxmass}
{massatomtypelines}
END
BONDSTRETCHTYPECODE
{stretchcount}    {maxstretch}
{bondstretchtypelines}
END
BONDANGLEBENDTYPECODE
{bendcount}    {maxbend}
{bondanglebendtypelines}
END
IMPDIHEDRALTYPECODE
{impcount}    {maximp}
{impdihedraltypelines}
END
TORSDIHEDRALTYPECODE
{torscount}    {maxtors}
{torsdihedraltypelines}
END
SINGLEATOMLJPAIR
{singleatomLJ}
END
MIXEDATOMLJPAIR
{mixedatomLJ}
END
SPECATOMLJPAIR
{specatomLJ}
END
'''.format(**data), file=f)
        return f.getvalue()
if __name__ == "__main__":
    from genmethods import parseIFP
    ffData = parseIFP('54A7', './')
    print(ffData.printIFP())