import random
import io
import sys
import numpy, genmethods, openbabel
from config import ATOMIC_CHARGE
import rotations
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
        num = 0
        for a in self.atms:
            if a.getAtmPos() > 0: num += 1
        return num
    def getNumBonds(self):
        num = 0
        for b in self.bonds:
            if b.getAtmAPos() > 0 and b.getAtmBPos() > 0: num +=1
        return num
    def getNumAngles(self):
        num = 0
        for b in self.angles:
            if b.getAtmAPos() > 0 and b.getAtmBPos() > 0 and b.getAtmCPos() > 0: num +=1
        return num
    def getNumImpropers(self):
        num = 0
        for b in self.impropers:
            if b.getAtmAPos() > 0 and b.getAtmBPos() > 0 and b.getAtmCPos() > 0 and b.getAtmDPos() > 0: num +=1
        return num
    def getNumDihedrals(self):
        num = 0
        for b in self.dihedrals:
            if b.getAtmAPos() > 0 and b.getAtmBPos() > 0 and b.getAtmCPos() > 0 and b.getAtmDPos() > 0: num +=1
        return num
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
    def getMoleculeCharge(self):
        charge = 0.
        for atm in self.getAtms():
            charge += atm.getCharge()
        return charge
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
        if self.getNumAtms() > 0:
            for atm in self.getAtms(): atm.setResName(self.compndName)
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
    
    ''' Method to center the molecule, either on a given atmIndex or a given set of XYZ coordinates'''
    def centreOn(self, data):
        if type(data) is int: # get the required data
            data = self.getAtm_Index(data).getXYZ()
        for atom in self.atms:
            atom.setXYZ([a-b for a, b in zip(atom.getXYZ(),data)])
    ''' Method to rotate the molecule, or defined part there of, using the given rotation matrix. '''
    def rotate(self, rotMat, atms=None):
        if atms is None:
            atms = []
            for num in range(self.getNumAtms()): atms.append(self.getAtms()[num].getAtmIndex())
        for atm in self.getAtms():
            if atm.getAtmIndex() in atms:
                atm.setXYZ(numpy.dot(rotMat,atm.getXYZ(vec=True)),vec=True)
            
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
            if atom.getAtmPos() > 0:
                print(atom.getPDBString(), file=f)
        for bond in self.getBonds():
            if bond.getAtmAIndex() > bond.getAtmBIndex():
                temp = bond.getAtmBIndex()
                bond.setAtmBIndex(bond.getAtmAIndex())
                bond.setAtmAIndex(temp)
            if bond.getAtmAPos() > 0 and bond.getAtmBPos() > 0:
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
        print('''TITLE
File          : 54a7.mtb
Force field   : 54A7 (condensed-phase simulations)
Reference     : Schmid et al. Eur. Biophys. J. 2011, 40, 843-856
File content  : Molecular topology building blocks (alpha amino acids, nucleic acids, lipids)
Format        : GROMOS11
Initial file  : WFVG, AM, CO, JD, Zuerich, September 2010
Time stamp    : PHH, Fri Aug  3 17:23:08 CEST 2012
Remark        : The main general changes from the 53A(B)6 force-field
                version to the 54A(B)7 force-field version involve
                - changed phi/psi torsional parameters of peptides
                - new choline CH3 atom type
                - new VdW parameters for Na+ and Cl-
                - addition of the AIB building block
Remark        : PHH, May 2011
                - finalized GROMOS11 file distribution
                - enforced sequential ordering by (central) atom numbers in covalent terms
                  for all files (no effect on GROMOS11 make_top; make_top enforced it anyway,
                  but now, we avoid the big list of warnings)
Modifications : [list below changes after May 2011 - with initials, location and date - and update time stamp]
PHH, 15.09.2011: Corrected the phi/psi dihedral potentials for residue AIB in the G96 and G11
                 54A7 and 54B7 mtb files (they had not been updated from 53A(B)6; found by Alpesh Malde)
PHH, 15.09.2011: Corrected the atom charges for residue DPPC in the G96 and G11 54A7 mtb files
                 (they had not been updated from 53A6; found by Alpesh Malde). Note that the
                 DPPC charges in 54B7 are kept the same as in 53B6, i.e. not updated (they differ
                 between 53A6 and 54A7, and the definition of a scheme for designing new 54B7 charges 
                 is neither obvious nor really urgent).
PHH, 09.11.2011: Reintroduced a FORCEFIELD block in all GROMOS11 files.
PHH, 09.11.2011: Changed atom name H3 by H2 in NH2 patch.
PHH, 09.11.2011: Introduced a copy of patch D5OH named 5OH, the latter meant for RNA instead of DNA.
PHH, 22.11.2011: Finalized the #@BLOCKTYPE comments in mtb files listing file name, residue
                 code (BLK=...), function (SOL,INI,TER,SVT), type (TYPE=APEP,BPEP,DNUC,RNUC,
                 HEXP,HEXU,MOLE), and full name (NAME=...): intended for later use by make_top.
PHH, 26.06.2012: Introduced MAKETOPVERSION blocks in all GROMOS11 files
PHH, 12.07.2012: Removed all LINKEXCLUSIONS and PHYSICALCONSTANTS blocks from GROMOS11 
                 auxiliary mtb files (now only included in the main mtb file). As a result 
                 (and also since MTBUILBLSOLVENT is only to be found there), make_top must
                 always be called with inclusion of the main mtb file (in addition to 
                 the possible auxiliary ones).
END
FORCEFIELD
54A7
END
MAKETOPVERSION
1.0
END
PHYSICALCONSTANTS
# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)
  0.1389354E+03
# HBAR: Planck's constant HBAR = H/(2* PI)
  0.6350780E-01
# SPDL: Speed of light (in nm/ps)
  2.9979245800E05
# BOLTZ: Boltzmann's constant
  8.31441E-03
END
LINKEXCLUSIONS
#nearest neighbour exclusions when linking
#NRNE
    2
END''', file=f)
        print("MTBUILDBLSOLUTE", file=f)
        print(self.getCompndName(), file=f)
        print("# atoms", file=f)
        print("%5d %5d" %(self.getNumAtms(), 0), file=f)
        for a in self.getAtms():
            atomString = "%5d %5s %5d %5d %10.6f %5d %4d " + " ".join( str(b) for b in a.getNonBondExcludeAtms())
            if a.getAtmPos() > 0:
                print(atomString %(a.getAtmPos(),a.getAtmName(),a.getIACM(),a.getMassNum(),a.getCharge(),a.getChargeGroupCode(),len(a.getNonBondExcludeAtms())), file=f)
        print("# bonds", file=f)
        print("%5d" % self.getNumBonds(), file=f)
        for b in self.getBonds():
            if b.getAtmAPos() > 0 and b.getAtmBPos() > 0:
                print("%5d %5d %5d" %(b.getAtmAPos(), b.getAtmBPos(), b.getBondTypeCode()), file=f)
        print("# angles", file=f)
        print("%5d" % self.getNumAngles(), file=f)
        for a in self.getAngles():
            if a.getAtmAPos() > 0 and a.getAtmBPos() > 0 and a.getAtmCPos() > 0:
                print("%5d %5d %5d %5d" %(a.getAtmAPos(), a.getAtmBPos(), a.getAtmCPos(), a.getAngleTypeCode()), file=f)
        print("# impropers", file=f)
        print("%5d" % self.getNumImpropers(), file=f)
        for a in self.getImpropers():
            if a.getAtmAPos() > 0 and a.getAtmBPos() > 0 and a.getAtmCPos() > 0 and a.getAtmDPos() > 0:
                print("%5d %5d %5d %5d %5d" %(a.getAtmAPos(), a.getAtmBPos(), a.getAtmCPos(), a.getAtmDPos(), a.getImpTypeCode()), file=f)
        print("# dihedrals", file=f)
        print("%5d" % self.getNumDihedrals(), file=f)
        for a in self.getDihedrals():
            if a.getAtmAPos() > 0 and a.getAtmBPos() > 0 and a.getAtmCPos() > 0 and a.getAtmDPos() > 0:
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
    ''' Types are dihedral, rigid, bonds '''
    def printGAMESS(self,dihed, angle, type="dihedral", solvent=True):
        atms = {'A':dihed.getAtmAPos(),'B':dihed.getAtmBPos(),
                'C':dihed.getAtmCPos(),'D':dihed.getAtmDPos()}
        info = {'runtyp':'OPTIMIZE','charge':int(self.getMoleculeCharge()),'solvent':'\n $PCM SOLVNT=WATER $END',
                'ifzmat':'3,{A},{B},{C},{D}'.format(**atms),'name':self.getCompndName(),'dihedralangle':angle,
                'freezevals':str(round(angle,3)),'zmat':'','nzvar':3*self.getNumAtms()-6}
        # use openbabel to generate a Z-matrix that can be manipulated
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdb","mopin")
        obMol = openbabel.OBMol()
        obConversion.ReadString(obMol, self.printPDB())
        zMat = obConversion.WriteString(obMol).split('\n')[3:-1]
        zmatLines = []
        for i in range(len(zMat)):
            zMat[i] = zMat[i].split()
        for atm in range(len(zMat)):
            line = zMat[atm][0]
            if int(zMat[atm][9]) != 0: line += '  ' + '  '.join(str(x) for x in zMat[atm][1:])
            elif int(zMat[atm][8]) != 0: line += '  ' + '  '.join(str(x) for x in (zMat[atm][1],zMat[atm][2],zMat[atm][3],zMat[atm][4],zMat[atm][7],zMat[atm][8]))
            elif int(zMat[atm][7]) != 0: line += '  ' + '  '.join(str(x) for x in zMat[atm][1:2])
            zmatLines.append(line)
        info['zmat'] = '\n'.join(x for x in zmatLines)
        # check that the dihedral is in the z-matrix, and change it if not
        if int(zMat[atms['D']-1][7]) != atms['C'] or int(zMat[atms['D']-1][8] != atms['B']) or int(zMat[atms['D']-1][9] != atms['A']):
            zMat[atms['D']-1][7] = str(atms['C'])
            zMat[atms['D']-1][8] = str(atms['B'])
            zMat[atms['D']-1][9] = str(atms['A'])
            zMat[atms['D']-1][5] = round(angle,3)
            zMat[atms['D']-1][3] = round(rotations.vectAngle(self.getAtmWithIndex(atms['D']-1).getXYZ(vec=True)-self.getAtmWithIndex(atms['C']-1).getXYZ(vec=True), self.getAtmWithIndex(atms['B']-1).getXYZ(vec=True)-self.getAtmWithIndex(atms['C']-1).getXYZ(vec=True))*180/numpy.pi,6)
            zMat[atms['D']-1][1] = round(numpy.linalg.norm(self.getAtmWithIndex(atms['D']-1).getXYZ(vec=True)-self.getAtmWithIndex(atms['C']-1).getXYZ(vec=True)),6)
        if type == 'rigid': info['runtyp'] = 'ENERGY'
        if not solvent: info['solvent'] = ''
        if type == 'bonds':
            for bond in self.getBonds():
                rigid = ',\n1,%d,%d' %(bond.getAtmAPos(),bond.getAtmBPos())
                freeze = round(numpy.linalg.norm(self.getAtmWithIndex(bond.getAtmAIndex()).getXYZ(vec=True)-self.getAtmWithIndex(bond.getAtmBIndex()).getXYZ(vec=True)),6)
                info['ifzmat'] += rigid
                info['freezevals'] += ',\n'+str(freeze)
        izmatLine = []
        for atm in range(1,len(zMat)):
            #print(zMat[atm])
            line = ''
            if int(zMat[atm][7]) != 0: line += '1,%s,%d' % (zMat[atm][7],atm+1)
            if int(zMat[atm][8]) != 0: line += ',2,%s,%s,%d' % (zMat[atm][8],zMat[atm][7],atm+1)
            if int(zMat[atm][9]) != 0: line += ',3,%s,%s,%s,%d' % (zMat[atm][9],zMat[atm][8],zMat[atm][7],atm+1)
            izmatLine.append(line)
        info['izmat'] =''+',\n'.join(x for x in izmatLine)
        f = io.StringIO()
        #PRINT'12345678901234567890123456789012345678901234567890123456789012345678901234567890
        print(''' $CONTRL SCFTYP=RHF RUNTYP={runtyp} COORD=ZMTMPC UNITS=ANGS ICHARG={charge} 
DFTTYP=B3LYP MAXIT=200 NZVAR={nzvar} $END{solvent}
 $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 NPFUNC=0 $END
 $SCF MAXDII=40 $END
 $STATPT NSTEP=500 DXMAX=0.05 $END
 $ZMAT 
IZMAT(1)={izmat}
IFZMAT(1)={ifzmat} 
FVALUE(1)={freezevals} $END
 $DATA
{name} calculation
C1
{zmat}

   $END'''.format(**info), file=f)
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
        # tracking data
        self.addedCharge = 0.
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
    def addCharge(self, charge):
        self.addedCharge += charge
        self.setCharge(self.getCharge()+charge)
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
        if self.getAtmPos() > 0:
            return formattedString
        else:
            return
        
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
        return [self.atmA-1,self.atmB-1,self.atmC-1,self.atmD-1]
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