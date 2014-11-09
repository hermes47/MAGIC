import random
import io
import sys
import numpy
''' Molecule class.
Collection of atoms, bonds, angles, impropers and dihedrals
'''
class Molecule(object):
    def __init__(self):
        self.compndName = "UNKNOWN_" + str(random.randint(1000000,9999999))
        self.title = ""
        self.authors = []
        self.atms = []
        self.bonds = []
        self.angles = []
        self.impropers = []
        self.dihedrals = []
        self.masterPDB = "MASTER"
        self.mtbTitle = []
        self.forcefield = ""
        self.physicalConstants = []
        self.linkExclusions = 0
        self.precedingExclusions = 0
        self.numJoiningPoints = 0
        self.joiningPoints = [] # stored as index of atom
        self.transverseBonds = [] # stored as index of bond
    ''' Method so that the Molecule can be printed vaguely nicely '''
    def __repr__(self, *args, **kwargs):
        dataString = (self.getCompndName(),
                      self.getNumAtms(),
                      self.getNumBonds(),
                      self.getNumAngles(),
                      self.getNumImpropers(),
                      self.getNumDihedrals())
        return ": ".join(str(v) for v in dataString)
    ''' Methods for getting individual data points from the molecule '''
    def getCompndName(self):
        return self.compndName
    def getTitle(self):
        return self.title
    def getAuthors(self):
        return self.authors
    def getNumAtms(self):
        return len(self.atms)
    def getNumBonds(self):
        return len(self.bonds)
    def getNumAngles(self):
        return len(self.angles)
    def getNumImpropers(self):
        return len(self.impropers)
    def getNumDihedrals(self):
        return len(self.dihedrals)
    def getAtm_Index(self, aIndex):
        return self.atms[aIndex]
    def getAtm_Pos(self, aPos):
        return self.atms[aPos - 1]
    def getAtms(self):
        return self.atms
    def getAtmWithIndex(self, index):
        for atm in self.getAtms():
            if atm.getAtmIndex() == index: return atm
    def getBond(self, bIndex):
        return self.bonds[bIndex]
    def getBonds(self):
        return self.bonds
    def getAngles(self):
        return self.angles
    def getImpropers(self):
        return self.impropers
    def getDihedrals(self):
        return self.dihedrals
    def getChargeGroups(self):
        chargeGroups = []
        tempChargeGroup = []
        for atm in self.getAtms():
            tempChargeGroup.append(atm.getAtmIndex())
            if atm.getChargeGroupCode() == 1:
                chargeGroups.append(tempChargeGroup)
                tempChargeGroup = []
        return chargeGroups
    def getMasterPDB(self):
        return self.masterPDB
    def getBondIDs(self):
        bondIDs = []
        for b in self.getBonds():
            bondIDs.append(b.getBondID())
        return bondIDs
    def getMTBTitle(self):
        return self.mtbTitle
    def getForceField(self):
        return self.forcefield
    def getPhysicalConstants(self):
        return self.physicalConstants
    def getLinkExclusions(self):
        return self.linkExclusions
    def getPrecedingExclusions(self):
        return self.precedingExclusions
    def getNumJoiningPoints(self):
        return self.numJoiningPoints
    def getJoiningPoints(self):
        return self.joiningPoints
    def getTransBonds(self):
        return self.transverseBonds
    def getTransBondsTuple(self, index):
        broken = tuple(map(int,self.transverseBonds[index].split(',')))
        return broken
    ''' Methods for getting collections of data points from the molecule '''
    def getPDBStart(self):
        return "COMPND    " + self.getCompndName()
    def getPDBEnd(self):
        #data = self.masterPDB.split()
        #return "%6s%9s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s\nEND" %(data[0],data[1],data[2],data[3],
        #                                                        data[4],data[5],data[6],data[7],
        #                                                        data[8],data[9],data[10],data[11],
        #                                                        data[12],)
        return "END"
    def getDipole(self):
        a = numpy.array([0.,0.,0.])
        for atm in self.getAtms():
            a += atm.getDipoleVector()
        return a
    def getDipoleMagnitude(self):
        return numpy.linalg.norm(self.getDipole())
    def getMolecularCoordinates(self, dipole=False):
        if dipole:
            return numpy.vstack((atm.getDipoleVector() for atm in self.getAtms()))
        else:
            return numpy.vstack((atm.getXYZ(vec=True) for atm in self.getAtms()))
    ''' Methods for setting individual data points in the molecule '''
    def setCompndName(self, sName):
        if sName.startswith("COMPND"):
            self.compndName = sName[10:]
        else: self.compndName = sName
    def setTitle(self, sTitle):
        if sTitle.startswith("TITLE"):
            self.title = sTitle[10:]
        else: self.title = sTitle
    def setMasterPDB(self, sMaster):
        self.masterPDB = sMaster
    def setForceField(self, ff):
        self.forcefield = ff
    def setLinkExclusions(self, le):
        self.linkExclusions = le
    def setPrecedingExclusions(self, pe):
        self.precedingExclusions = pe
    def setNumJoiningPoints(self, jp):
        self.numJoiningPoints = jp
    ''' Methods for adding items to the lists '''
    def addAtm(self, atmData, copyAtm=False):
        if not copyAtm:
            at = Atom(atmData)
            self.atms.append(at)
        elif copyAtm:
            self.atms.append(atmData)
    def addBond(self, bondData, copyBond=False):
        if not copyBond:
            bd = Bond(bondData)
        elif copyBond:
            bd = bondData
        if bd.getBondID() not in self.getBondIDs(): # only add the bond if it does not yet exist in the molecule
            self.bonds.append(bd)
        else: 
            print("BondID: "+bd.getBondID()+" rejected as already exists", file=sys.stderr)
    def addAngle(self, angleData, copyAngle=False):
        if not copyAngle:
            an = Angle(angleData)
            self.angles.append(an)
        elif copyAngle:
            self.angles.append(angleData)
    def addImproper(self, improperData, copyImp=False):
        if not copyImp:
            im = Improper(improperData)
            self.impropers.append(im)
        elif copyImp:
            self.impropers.append(improperData)
    def addDihedral(self, dihedralData, copyDihed=False):
        if not copyDihed:
            di = Dihedral(dihedralData)
            self.dihedrals.append(di)
        elif copyDihed:
            self.dihedrals.append(dihedralData)
    def addAuthor(self, sAuthor):
        if not sAuthor.startswith("Author"):
            self.authors.append(sAuthor)
        else: self.authors.append(sAuthor[10:])
    def addMTBTitleLine(self, line):
        self.mtbTitle.append(line)
    def addPhysicalConstant(self, const):
        self.physicalConstants.append(const)
    def addJoiningPoint(self, jp, copyJP=False):
        if not copyJP:
            self.joiningPoints.append(jp - 1) # comes in as a position, stored as an Index
        elif copyJP:
            self.joiningPoints.append(jp)
        self.numJoiningPoints += 1
    def setJoiningPoint(self, l):
        self.joiningPoints = l
        self.numJoiningPoints = len(l)
    def addTransBond(self, bd):
        # if the argument type is a string, break it at the comma, cast to ints, and rejoin as a tuple
        self.transverseBonds.append(bd)
    def setTransBonds(self, l):
        self.transverseBonds = l 
    ''' Methods for higher level operations '''
    
    ''' method to centre the molecule, either on a given atmIndex or a given set of XYZ coordinates'''
    def centreOn(self, data):
        if type(data) is int: # get the required data
            data = self.getAtm_Index(data).getXYZ()
        for atom in self.atms:
            atom.setXYZ([a-b for a, b in zip(atom.getXYZ(),data)])
    ''' method for deleting an atom from the molecule. Also deletes any bonds that the molecule is in 
    and shifts the bond indexs to account for the deleted bonds'''
    ##### DEPRECITED #####
    def delAtm(self, atmIndex):
        tempAtms = []
        # create a new list of atoms which is missing the one that needs to be deleted
        for atom in self.atms:
            if atom.getAtmIndex() != atmIndex:
                tempAtms.append(atom)
        # create a new list of bonds which is missing those bonds that involve the atom to be deleted
        tempBonds = []
        for bond in self.bonds:
            if not bond.containsAtmIndex(atmIndex):
                tempBonds.append(bond)
        # create a change list detailing what each bond index needs to change to
        cmpList = []
        shift = 0
        for i in range(len(tempAtms)):
            if self.atms[i+shift].getAtmIndex() != tempAtms[i].getAtmIndex():
                cmpList.append((i+shift,-1))
                shift += 1
            cmpList.append((i+shift,i)) # (x,y) such that in an existing bond, x becomes y
        # change the bond indexs for each bond in the tempBonds
        for bond in tempBonds:
            bond.setAtmAIndex(cmpList[bond.getAtmAIndex()][1])
            bond.setAtmBIndex(cmpList[bond.getAtmBIndex()][1])
        # fix the atmIndexs within each atom
        for i in range(len(tempAtms)):
            tempAtms[i].setAtmIndex(i)
        # fix the MasterPDB string to account for the lost atoms
        newMaster = list(map(int, self.getMasterPDB().split()[1:]))
        for i in range(len(newMaster)):
            if newMaster[i] > 0:
                newMaster[i] -= 1
        newMasterString = "MASTER    " + " ".join(i for i in map(str, newMaster)) + "\n"
        self.setMasterPDB(newMasterString)
        # change the atms and bonds lists to their new ones, and fix the counts
        self.atms = tempAtms
        self.bonds = tempBonds
        self.numAtms = len(tempAtms)
        self.numBonds = len(tempBonds)
    ''' method for shifting the indices of each atm, and the bonded info, by a given amount '''
        # DEPRECIATED
    def shiftIndex(self, shiftCount):
        # shift the atoms
        for atom in self.atms:
            newIndex = atom.getAtmIndex() + shiftCount
            atom.setAtmIndex(newIndex)
        # shift the bond indexing
        for bond in self.bonds:
            newAIndex = bond.getAtmAIndex() + shiftCount
            newBIndex = bond.getAtmBIndex() + shiftCount
            bond.setAtmAIndex(newAIndex)
            bond.setAtmBIndex(newBIndex)
    ''' method for printing out a PDB of the molecule '''
    def printPDB(self):
        # open string stream for write
        f = io.StringIO()
        print(self.getPDBStart(), file=f)
        for atom in self.getAtms():
            print(atom.getPDBString(), file=f)
        for bond in self.getBonds():
            if bond.getAtmAIndex() > bond.getAtmBIndex():
                temp = bond.getAtmBIndex()
                bond.setAtmBIndex(bond.getAtmAIndex())
                bond.setAtmAIndex(temp)
            print(bond.getPDBString(), file=f)
        #for i in range(1,self.getNumAtms()+1):
        #    bondedAtoms = []
        #    for bond in self.getBonds():
        #        if bond.getAtmAPos() == i:
        #            bondedAtoms.append(bond.getAtmBPos())
        #        elif bond.getAtmBPos() == i:
        #            bondedAtoms.append(bond.getAtmAPos())
        #    bondedAtoms = list(set(bondedAtoms))
        #    bondedAtoms.sort()
        #    words = "CONECT%5d" %i +" "+" ".join('%5d' %a for a in bondedAtoms)
        #    print(words)
        print(self.getPDBEnd(), file=f)
        return f.getvalue()
        
    ''' method for printing out the MTB of the molecule '''
    def printMTB(self):
        # open string stream for write
        f = io.StringIO()
        print("MTBUILDBLSOLUTE", file=f)
        print(self.getCompndName(), file=f)
        print("# atoms", file=f)
        print("%5d %5d" %(self.getNumAtms(), 0), file=f)
        for a in self.getAtms():
            atomString = "%5d %5s %5d %5d %10.6f %5d %4d " + " ".join( str(b) for b in a.getNonBondExcludeAtms())
            print(atomString %(a.getAtmPos(),a.getAtmName(),a.getIACM(),a.getMassNum(),a.getCharge(),a.getChargeGroupCode(),len(a.getNonBondExcludeAtms())), file=f)
        print("# bonds", file=f)
        print("%5d" % self.getNumBonds(), file=f)
        for b in self.getBonds():
            print("%5d %5d %5d" %(b.getAtmAPos(), b.getAtmBPos(), b.getBondTypeCode()), file=f)
        print("# angles", file=f)
        print("%5d" % self.getNumAngles(), file=f)
        for a in self.getAngles():
            print("%5d %5d %5d %5d" %(a.getAtmAPos(), a.getAtmBPos(), a.getAtmCPos(), a.getAngleTypeCode()), file=f)
        print("# impropers", file=f)
        print("%5d" % self.getNumImpropers(), file=f)
        for a in self.getImpropers():
            print("%5d %5d %5d %5d %5d" %(a.getAtmAPos(), a.getAtmBPos(), a.getAtmCPos(), a.getAtmDPos(), a.getImpTypeCode()), file=f)
        print("# dihedrals", file=f)
        print("%5d" % self.getNumDihedrals(), file=f)
        for a in self.getDihedrals():
            print("%5d %5d %5d %5d %5d" %(a.getAtmAPos(), a.getAtmBPos(), a.getAtmCPos(), a.getAtmDPos(), a.getDihedTypeCode()), file=f)
        print("# LJ Exceptions", file=f)
        print(0, file=f)
        print("END", file=f)
        return f.getvalue()
    def printChargeData(self):
        f = io.StringIO()
        for a in self.getAtms():
            print("%s    %6.3f %6.3f %6.3f " %(a.getAtmName()[0], a.getX(), a.getY(), a.getZ()), file=f)
        return f.getvalue()

''' Atom class.
Description:
describes an atom as taken from a PDB file. Has no information of connectivity, but has positional information as well as all the information
required to regenerate the PDB file line that was used to generate the atom. Atom positions, ie those that are relevant to a molecule, are stored
as they come from the PDB file. Thus and atom with position x, is the x'th atom in the list, accessed from index x-1.

Variables:
lineType:    what type of line in PDB file the atom comes from. (HETATM or ATOM). Does nothing except regenerate the PDB line. Type = str.
atmIndex:    the position that the atom takes in the molecule's atom list. Starts counting from 1 as is taken directly from PDB file. Type = int.
atmName :    the name of the atom. Does nothing except regenerate the PDB file. Type = str.
resName :    the name of the residue. Does nothing except regenerate the PDB file. Type = str.
resIndex:    the index of the residue, as taken directly from PDB file. Does nothing except regenerate the PDB file. Type = int.
xPos    :    the position of the atom along the X-axis in space. Type = float.
yPos    :    the position of the atom along the Y-axis in space. Type = float.
zPos    :    the position of the atom along the Z-axis in space. Type = float.
crys1   :    the first crystalographic data point taken from the PDB. Does nothing except regenerate the PDB file. Type = float.
crys2   :    the second crystalographic data point taken from the PDB. Does nothing except regenerate the PDB file. Type = float.
atmType :    the elemental symbol for the atom. Is used to help determine whether a given atom will be used as a joining point or not. Type = str
'''
class Atom(object):
    ''' initialisation method'''
    def __init__(self, pdbString):
        # initialise all the data points from the pdb string
        # NEEDS reworking to not work on split but to work on column numbers
        self.lineType, self.atmIndex, self.atmName, self.resName, self.resIndex, self.xPos, self.yPos, self.zPos, self.crys1, self.crys2, self.atmType = pdbString.split()
        #self.lineType, self.atmIndex, self.atmName, self.resName, self.resIndex, self.xPos, self.yPos, self.zPos, self.crys1, self.crys2 = pdbString.split()
        # cast to int and float those that are needed to be
        self.atmIndex = int(self.atmIndex)
        self.resIndex = int(self.resIndex)
        self.xPos = float(self.xPos)
        self.yPos = float(self.yPos)
        self.zPos = float(self.zPos)
        self.crys1 = float(self.crys1)
        self.crys2 = float(self.crys2)
        # data unique to the mtb file
        self.IACM = 0
        self.massNum = 0
        self.charge = 0.
        self.chargeGroupCode = 0
        self.numNonBondExclude = 0
        self.nonBondExcludeAtms = []
        self.tier = 0
    def __repr__(self):
        return repr((self.atmIndex, self.charge, self.chargeGroupCode))
    ''' Methods for getting the indivdual data points from with the atom class '''
    def getLineType(self):
        return self.lineType
    def getAtmIndex(self): # returns the index for list access of the atom, starts at 0
        return self.atmIndex - 1
    def getAtmPos(self): # returns the position of the atom in the list, starts at 1
        return self.atmIndex
    def getAtmName(self):
        return self.atmName
    def getResName(self):
        return self.resName
    def getResIndex(self):
        return self.resIndex
    def getX(self):
        return round(self.xPos,3)
    def getY(self):
        return round(self.yPos,3)
    def getZ(self):
        return round(self.zPos,3)
    def getCrys1(self):
        return round(self.crys1,3)
    def getCrys2(self):
        return round(self.crys2,3)
    def getAtmType(self):
        return self.atmType
    def getIACM(self):
        return self.IACM
    def getMassNum(self):
        return self.massNum
    def getCharge(self):
        return self.charge
    def getChargeGroupCode(self):
        return self.chargeGroupCode
    def getNumNonBondExclude(self):
        return len(self.nonBondExcludeAtms)
    def getNonBondExcludeAtms(self):
        return self.nonBondExcludeAtms
    def getTier(self):
        return self.tier
    ''' Methods for getting combinations of data points from the atom class '''
    def getXYZ(self, vec=False):
        if vec:
            return numpy.array([self.xPos,self.yPos,self.zPos],dtype=numpy.float64)
        else:
            return [self.xPos,self.yPos,self.zPos]
    def getDipoleVector(self):
        return self.getCharge()*numpy.array([self.getX(),self.getY(),self.getZ()])
    ''' Methods for setting/changing data points within the atom class'''
    def setAtmIndex(self, i):
        self.atmIndex = i + 1 # sets the index of the atom, based on the integer index passed. as atmIndex stores the position of the atom, is +1
    def setAtmPos(self, i):
        self.atmIndex = i # sets the position of the atom, based on the integer position passed.
    def setAtmName(self, s):
        self.atmName = s
    def setResName(self, s):
        self.resName = s
    def setResIndex(self, i):
        self.resIndex = i
    def setX(self, f):
        self.xPos = f
    def setY(self, f):
        self.yPos = f
    def setZ(self, f):
        self.zPos = f
    def setCrys1(self, f):
        self.crys1 = f
    def setCrys2(self, f):
        self.crys2 = f
    def setXYZ(self, l, vec=False):
        if vec:
            if l.shape == (1,3):
                l = l.T
            self.setX(float(l[0]))
            self.setY(float(l[1]))
            self.setZ(float(l[2]))
        else:
            self.setX(l[0])
            self.setY(l[1])
            self.setZ(l[2])
    def setIACM(self, iacm):
        self.IACM = iacm
    def setMassNum(self, mass):
        self.massNum = mass
    def setCharge(self, ch):
        self.charge = ch
    def setChargeGroupCode(self, cgm):
        self.chargeGroupCode = cgm
    def setNumNonBondExclude(self, num):
        self.numNonBondExclude = num
    def setNonBondExcludeAtms(self, l):
        self.nonBondExcludeAtms = l
    def addNonBondExlude(self, a):
        self.nonBondExcludeAtms.append(a)
    def addMTBData(self, dataList):
        self.setIACM(int(dataList[0]))
        self.setMassNum(int(dataList[1]))
        self.setCharge(float(dataList[2]))
        self.setChargeGroupCode(int(dataList[3]))
        if int(dataList[4]) != 0:
            self.setNonBondExcludeAtms(list(map(int,dataList[5:])))
    def setTier(self, t):
        self.tier = t
    ''' Methods for higher order operations '''
    def getPDBString(self):
        #123456|78901|2|3456|7890|1|23456|7890|12345678|90123456|78901234|567890|123456|789012|3456|78|90
        #HETATM|   29| |H10 |MCW5| |    0|    |   2.552|   1.304|  -1.322|  1.00|  0.00|      |    | H|  
        formattedString = "%-6s%5d %-4s%5s%5d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  " %(self.getLineType(), 
                                                                               self.getAtmPos(), 
                                                                               self.getAtmName(), 
                                                                               self.getResName(),
                                                                               self.getResIndex(),
                                                                               self.getX(),
                                                                               self.getY(),
                                                                               self.getZ(),
                                                                               self.getCrys1(),
                                                                               self.getCrys2(), 
                                                                               self.getAtmType())
        return formattedString
        
''' Bond class.
'''
class Bond(object):
    ''' initisalisation method '''
    def __init__(self, bondData):
        # initialise all the data from the bond data string or list
        if type(bondData) is str:
            self.atmA, self.atmB, self.bondLength = bondData.split()
        elif type(bondData) is list or type(bondData) is tuple:
            self.atmA = bondData[0]
            self.atmB = bondData[1]
            self.bondLength = bondData[2]
        # cast to ints and floats those that need to be, to make sure
        self.atmA = int(self.atmA)
        self.atmB = int(self.atmB)
        self.bondLength = float(self.bondLength)
        # make sure atmA is the lower index one
        if self.atmA > self.atmB:
            temp = self.atmA
            self.atmA = self.atmB
            self.atmB = temp
        self.bondTypeCode = 0
    ''' Methods for getting the individual data points from within the bond class '''
    def getAtmAIndex(self):
        return self.atmA - 1
    def getAtmAPos(self):
        return self.atmA
    def getAtmBIndex(self):
        return self.atmB - 1
    def getAtmBPos(self):
        return self.atmB
    def getBondLength(self):
        return self.bondLength
    def getBondTypeCode(self):
        return self.bondTypeCode
    ''' Methods for getting combinations of data points from the bond class '''
    def getAtms(self):
        return [self.atmA,self.atmB]
    def getPDBString(self):
        formattedString = "%6s%5d%5d" % ("CONECT",self.getAtmAPos(),self.getAtmBPos())
        functions = [self.getAtmAPos(),
                     self.getAtmBPos(),
                     ]
        s = "CONECT "
        for func in functions:
            try:
                s += func + " "
            except TypeError:
                s += str(func) + " "
        return formattedString
    def getBondID(self):
        return str(self.getAtmAIndex())+","+str(self.getAtmBIndex())
    ''' Methods for setting/changing data point within the bond class '''
    def setAtmAIndex(self, i):
        self.atmA = i + 1
    def setAtmAPos(self, i):
        self.atmA = i
    def setAtmBIndex(self, i):
        self.atmB = i + 1
    def setAtmBPos(self, i):
        self.atmB = i
    def setBondLength(self, f):
        self.bondLength = f
    def setBondTypeCode(self, mcb):
        self.bondTypeCode = mcb
    ''' Methods for higher level operations '''
    def containsAtmIndex(self, aIndex):
        if self.getAtmAIndex() == aIndex or self.getAtmBIndex() == aIndex: return True
        else: return False
    def containsAtmPos(self, aPos):
        if self.getAtmBPos() == aPos or self.getAtmAPos() == aPos: return True
        else: return False

''' Angle class.
'''
class Angle(object):
    ''' initialisation method '''
    def __init__(self, angleData):
        # everything need by the Angle is given in the angleData
        self.atmA, self.atmB, self.atmC, self.angleTypeCode = angleData
    ''' Methods for getting datapoints from the angle class '''
    def getAtmAIndex(self):
        return self.atmA - 1
    def getAtmAPos(self):
        return self.atmA
    def getAtmBIndex(self):
        return self.atmB - 1
    def getAtmBPos(self):
        return self.atmB
    def getAtmCIndex(self):
        return self.atmC - 1
    def getAtmCPos(self):
        return self.atmC
    def getAngleTypeCode(self):
        return self.angleTypeCode
    def getAtms(self):
        return [self.atmA,self.atmB,self.atmC]
    ''' Methods for setting/changind data points within the angle class '''
    def setAtmAIndex(self, i):
        self.atmA = i + 1
    def setAtmAPos(self, i):
        self.atmA = i
    def setAtmBIndex(self, i):
        self.atmB = i + 1
    def setAtmBPos(self, i):
        self.atmB = i
    def setAtmCIndex(self, i):
        self.atmC = i + 1
    def setAtmCPos(self, i):
        self.atmC = i
    def setAngleTypeCode(self, mcb):
        self.angleTypeCode = mcb
    
''' Improper class.
'''
class Improper(object):
    ''' initialisation method '''
    def __init__(self, impData):
        # everything need by the improper is given in the impData
        self.atmA, self.atmB, self.atmC, self.atmD, self.impTypeCode = impData
    ''' Methods for getting datapoints from the improper class '''
    def getAtmAIndex(self):
        return self.atmA - 1
    def getAtmAPos(self):
        return self.atmA
    def getAtmBIndex(self):
        return self.atmB - 1
    def getAtmBPos(self):
        return self.atmB
    def getAtmCIndex(self):
        return self.atmC - 1
    def getAtmCPos(self):
        return self.atmC
    def getAtmDIndex(self):
        return self.atmD - 1
    def getAtmDPos(self):
        return self.atmD
    def getImpTypeCode(self):
        return self.impTypeCode
    def getAtms(self):
        return [self.atmA,self.atmB,self.atmC,self.atmD]
    ''' Methods for setting/changind data points within the improper class '''
    def setAtmAIndex(self, i):
        self.atmA = i + 1
    def setAtmAPos(self, i):
        self.atmA = i
    def setAtmBIndex(self, i):
        self.atmB = i + 1
    def setAtmBPos(self, i):
        self.atmB = i
    def setAtmCIndex(self, i):
        self.atmC = i + 1
    def setAtmCPos(self, i):
        self.atmC = i
    def setAtmDIndex(self, i):
        self.atmD = i + 1
    def setAtmDPos(self, i):
        self.atmD = i
    def setImpTypeCode(self, mcb):
        self.impTypeCode = mcb

''' Dihedral class.
'''
class Dihedral(object):
    ''' initialisation method '''
    def __init__(self, dihedData):
        # everything need by the improper is given in the dihedData
        self.atmA, self.atmB, self.atmC, self.atmD, self.dihedTypeCode = dihedData
    ''' Methods for getting datapoints from the dihedral class '''
    def getAtmAIndex(self):
        return self.atmA - 1
    def getAtmAPos(self):
        return self.atmA
    def getAtmBIndex(self):
        return self.atmB - 1
    def getAtmBPos(self):
        return self.atmB
    def getAtmCIndex(self):
        return self.atmC - 1
    def getAtmCPos(self):
        return self.atmC
    def getAtmDIndex(self):
        return self.atmD - 1
    def getAtmDPos(self):
        return self.atmD
    def getDihedTypeCode(self):
        return self.dihedTypeCode
    def getAtms(self):
        return [self.atmA,self.atmB,self.atmC,self.atmD]
    ''' Methods for setting/changind data points within the dihedral class '''
    def setAtmAIndex(self, i):
        self.atmA = i + 1
    def setAtmAPos(self, i):
        self.atmA = i
    def setAtmBIndex(self, i):
        self.atmB = i + 1
    def setAtmBPos(self, i):
        self.atmB = i
    def setAtmCIndex(self, i):
        self.atmC = i + 1
    def setAtmCPos(self, i):
        self.atmC = i
    def setAtmDIndex(self, i):
        self.atmD = i + 1
    def setAtmDPos(self, i):
        self.atmD = i
    def setDihedTypeCode(self, mcb):
        self.dihedTypeCode = mcb