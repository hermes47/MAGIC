from math import sqrt
from random import randrange
from config import NAME_LENGTH, RANDOM_NAME_CHARS, MTB_PDB_ROOT_NAMES, NUM_TIERS, CHECKED_ATOMS, LOADED_MOLECULES, ROOT_DIRECTORY, STRING_BREAK_CHARACTERS, ATOMIC_RADII, CHARGE_RULES
from copy import copy
from builtins import range
import numpy
import os
''' General purpose methods that aren't specific  to any given class '''
  
''' Method to calculate the distance between 2 sets of XYZ data.
    Returns a double.'''
def calcLength(xyzA, xyzB):
    length = 0.
    for i in range(len(xyzA)):
        length += (xyzA[i] - xyzB[i])**2
    return sqrt(length)

''' Generate a NAME_LENGTH character random string. Characters taken from the config
    character list. Return the string, if it isn't already in use. Otherwise
    generate another.'''
def genRandomString():
    while True:
        randStr = ""
        savedFiles = []
        for _, dirs, _ in os.walk(ROOT_DIRECTORY+'/OutputFiles/'):
            for direct in dirs:
                savedFiles.append(direct)
        for _ in range(NAME_LENGTH):
            randStr += RANDOM_NAME_CHARS[randrange(0,len(RANDOM_NAME_CHARS))]
        if (randStr not in MTB_PDB_ROOT_NAMES) and (randStr not in savedFiles):
            MTB_PDB_ROOT_NAMES.append(randStr)
            return randStr

''' Break string s into length sized pieces, and return a list of the pieces '''
def fragmentString(s, length=NAME_LENGTH):
    l = []
    for c in range(0,len(s),length):
        l.append(s[c:c+length])
    return l

''' Load a given file type with the given rootName '''
def loadFile(rootName, folder, fileExtension, removeComments='#'):
    with open(folder+rootName+fileExtension) as fh:
        file = fh.readlines()
    shortFile = []
    if removeComments:
        for line in file:
            if not line.startswith(removeComments): # ignore all the full line comments immediately
                shortFile.append(line[:-1]) # and get rid of the trailing new line, just cos
    return shortFile

''' Method for parsing an IFP file. '''
def parseIFP(rootName, folder):
    ifpData = loadFile(rootName, "Input/", ".ifp", True)

''' Method for parsing a PDB file.
Takes a rootname and a molecule as arguments.
Parses the PDB file and builds it into the molecule.
Returns the molecule. '''
def parsePDB(rootName, mol, hasMTB=False, passedData=None):
    # load the PDB file into memory
    if not passedData:
        pdbData = loadFile(rootName, "Input/PDBFiles/", ".pdb")
    else: pdbData = passedData
    if not hasMTB:
        # read the PDB file line by line
        for line in pdbData:
            # if the line is the compound name, set the name correctly
            if line.startswith("COMPND"):
                parsedLine = line[:-1]
                mol.setCompndName(parsedLine)
            # if the line is the title, set the title correctly
            elif line.startswith("TITLE"):
                parsedLine = line[:-1]
                mol.setTitle(parsedLine)
            # if the line contains author information, add it to the author list
            elif line.startswith("AUTHOR"):
                parsedLine = line[:-1]
                mol.addAuthor(parsedLine)
            # if the line is an atom line, add the atom to the pdbAtmData
            elif line.startswith("HETATM") or line.startswith("ATOM"):
                #lineType, atmIndex, atmName, resName, resIndex, x, y, z, crys1, crys2, atmType
                mol.addAtm(line)
            # if the line is connect data, add each bond to the pdbBondData, making sure that it is a new, unique bond
            elif line.startswith("CONECT"):
                cols = list(map(int, line.split()[1:]))
                rootAtm = cols[0]
                connAtms = cols[1:]
                rootCoords = mol.getAtm_Pos(rootAtm).getXYZ()
                for ca in connAtms:
                    if rootAtm > ca: # arrange so that the lowest index comes first
                        connectString = str(ca)+","+str(rootAtm)
                    else: connectString = str(rootAtm)+","+str(ca)
                    if connectString not in mol.getBondIDs(): # is this check redundant as the addBond method already checks?
                        caCoords = mol.getAtm_Pos(ca).getXYZ()
                        length = calcLength(rootCoords, caCoords)
                        mol.addBond((rootAtm, ca, length))
            # finally if the line is the MASTER line, set that in the molecule
            elif line.startswith("MASTER"):
                mol.setMasterPDB(line)
    else: print("MTB parsing prior to PDB parsing not yet implemented")
    
            
''' Method for parsing an MTB file.
Takes a rootName and a molecule as arguments. Optional argument hasPDB=True for if done PDB parsing or not.
Parse the MTB file and adds the info to the required atoms, bonds, etc.
Returns the molecule.'''
def parseMTB(rootName, mol, hasPDB=True):
    # load the mtb file into memory
    mtbData = loadFile(rootName, "Input/MTBFiles/", ".mtb")
    # now read in the mtb file line by line
    boolBlock = {"TITLE":False,"FORCEFIELD":False,"PHYSICALCONSTANTS":False,"LINKEXCLUSIONS":False,"MTBUILDBLSOLUTE":False}
    counts = {"atoms":0,"bonds":0,"dihedrals":0,"impropers":0,"angles":0,"LJExceptions":0,"locating":0}
    maximums = {"atoms":0,"bonds":0,"dihedrals":0,"impropers":0,"angles":0,"LJExceptions":0}
    for line in mtbData:
        # close out the block if line is END
        if line == "END":
            for i in boolBlock: boolBlock[i] = False
        else: # otherwise do all the things
            # Load things in depending on what block we're in.
            if boolBlock["TITLE"]:
                # append each line of the title block to the end of the molecules mtbTitle list
                mol.addMTBTitleLine(line)
            elif boolBlock["FORCEFIELD"]:
                # set the forcefield of the molecule
                mol.setForceField(line)
            elif boolBlock["PHYSICALCONSTANTS"]:
                # add the physical constant
                mol.addPhysicalConstant(line)
            elif boolBlock["LINKEXCLUSIONS"]:
                # add the link exclusions
                mol.setLinkExclusions(int(line))
            elif boolBlock["MTBUILDBLSOLUTE"]:
                # first line in this block gives name of the molecule
                if counts["locating"] == 0:
                    mol.setCompndName(line)
                    counts["locating"] += 1
                # second line gives the number of atoms, and number of preceding exclusions
                elif counts["locating"] == 1:
                    numAtms, numPreExclude = map(int, line.split())
                    maximums["atoms"] = numAtms # note: this assumes that the number of atoms in both the pdb
                                                # and mtb are the same. Is a possible point of breaking that could
                                                # be handled at some point
                    mol.setPrecedingExclusions(numPreExclude)
                    counts["locating"] += 1
                # next numAtms lines give the atom data
                elif counts["atoms"] < maximums["atoms"]:
                    atmData = line.split()  # split the line. Fragments are: ATOM, ATOMNAME, IACM, MASSCODE, CHARGE, 
                                            # CHARGEGROUP, numEXCLUDEDNEIGHBOURS, excludedneighbour1...
                    # delete any comments from the line
                    if "#" in atmData: del atmData[atmData.index("#"):]
                    # assume atm, and atmName are the same as the already loaded pdb file and add the MTB data to the atom
                    mol.getAtm_Pos(int(atmData[0])).addMTBData(atmData[2:])
                    counts["atoms"] += 1
                # next line gives the number of bonds
                elif counts["locating"] == 2:
                    maximums["bonds"] = int(line) 
                    # assume numBonds is same between PDB and MTB, as only print a warning if it's not true
                    if maximums["bonds"] != mol.getNumBonds(): print("MTB does nat have same number of bonds as PDB")
                    counts["locating"] += 1
                # next numBonds lines give the bond Data
                elif counts["bonds"] < maximums["bonds"]:
                    bondData = line.split()
                    # delete any comments from the line
                    if "#" in bondData: del bondData[bondData.index("#"):]
                    # everything is now ints, so make it so
                    bondData = tuple(map(int, bondData))
                    bondCode = str(bondData[0]-1)+","+str(bondData[1]-1) # bondIDs stored based on Index
                    if bondData[1] < bondData[0]: bondCode = str(bondData[1])+","+str(bondData[0]) # shouldn't ever happen
                    # as bonds are added to the end of the bond list, the index of the bondID will be 
                    # the same as the index of the bond it relates to
                    if bondCode not in mol.getBondIDs(): mol.addBond((min(bondData[0],bondData[1]),max(bondData[0],bondData[1]),0))
                    mol.getBond(mol.getBondIDs().index(bondCode)).setBondTypeCode(bondData[2])
                    counts["bonds"] += 1
                # next line gives the number of angles
                elif counts["locating"] == 3:
                    maximums["angles"] = int(line)
                    counts["locating"] += 1
                # next numAngles lines give the angle data
                elif counts["angles"] < maximums["angles"]:
                    angleData = line.split()
                    # delete any comments from the line
                    if "#" in angleData: del angleData[angleData.index("#"):]
                    # cast to ints
                    angleData = tuple(map(int, angleData))
                    # add the angle
                    mol.addAngle(angleData)
                    counts["angles"] +=1
                # next line give the number of impropers
                elif counts["locating"] == 4:
                    maximums["impropers"] = int(line)
                    counts["locating"] += 1
                # next numImproper lines give the improper data
                elif counts["impropers"] < maximums["impropers"]:
                    impData = line.split()
                    if "#" in impData: del impData[impData.index("#"):]
                    impData = tuple(map(int,impData))
                    mol.addImproper(impData)
                    counts["impropers"] += 1
                # next line gives the number of dihedrals
                elif counts["locating"] == 5:
                    maximums["dihedrals"] = int(line)
                    counts["locating"] += 1
                # next numDihedral lines give the dihedral data
                elif counts["dihedrals"] < maximums["dihedrals"]:
                    dihedData = line.split()
                    if "#" in dihedData: del dihedData[dihedData.index("#"):]
                    dihedData = tuple(map(int,dihedData))
                    mol.addDihedral(dihedData)
                    counts["dihedrals"] += 1
                else:    
                    pass
                
        # Determine if the line is the start of a block
        if line in boolBlock: 
            boolBlock[line] = True
    
''' Method for parsing an OVL file.
Takes a rootName and molecule as arguments.
Parse the OVL file and add the info to the molecule and atoms
Returns the molecule.'''
def parseOVL(rootName, mol):
    # load the file
    ovlData = loadFile(rootName, "Input/OVLFiles/", ".ovl")
    # read the file line by line. Doesn't loop because all ovl files have the same number of lines
    # first line is which atoms the joining points are, based on position
    for jp in map(int, ovlData[0].split()): 
        mol.addJoiningPoint(jp)
    # second line are the transverse bonds, ie those going between the kept fragment and anything added on
    # given to correspond to the order of the given joining point atoms
    for bd in ovlData[1].split():
        broken = tuple(map(int,bd.split(',')))
        mol.addTransBond(str(broken[0]-1)+","+str(broken[1]-1))
    # remaining lines are the tiers for the atoms. we loop across them as the number of levels may change
    level = 1
    for line in ovlData[2:2+NUM_TIERS]:
        atoms = tuple(map(int, line.split()))
        # first value is number of atoms at that tier, so is ignored
        for at in atoms[1:]:
            mol.getAtm_Pos(at).setTier(level)
        level += 1

''' Method for tree searching through a molecule, from a given starting point to determine
a list of atmIndex that should be deleted as they are the overlapping part of a molecule.
Recurses on itself until the tree has been searched all the way.
Returns a list of atmIndexs'''
def graphSearch(rootAtm, mol, resetChecked=False, searchType="delete", prescreened=[]):
    if resetChecked: del CHECKED_ATOMS[:]
    if len(prescreened) > 0:
        for pre in prescreened: CHECKED_ATOMS.append(pre)
    returnList = [] # list of atoms for return
    # find the other atoms involved with bonds to rootAtm
    bonded = [] # list of atom index bonded to rootAtm
    for bond in mol.getBonds():
        if bond.getAtmAIndex() == rootAtm: 
            bonded.append(bond.getAtmBIndex())
        elif bond.getAtmBIndex() == rootAtm: 
            bonded.append(bond.getAtmAIndex())
        else: pass # ignore bonds not involving rootAtm
    returned = [] # list of atoms that has been recursively passed back
    if rootAtm not in CHECKED_ATOMS:
        # add rootAtm to the checked list
        CHECKED_ATOMS.append(rootAtm)
        if searchType is "delete":
            # if the root atom is not tier 1 or 2, mark it for deletion
            if not 0 < mol.getAtm_Index(rootAtm).getTier() < 3:
                returnList.append(rootAtm)
            # determine if the bonded atoms can be deleted
            for atom in bonded:
                if not 0 < mol.getAtm_Index(atom).getTier() < 3:
                    for at in graphSearch(atom, mol): returned.append(at) # recurse
        elif searchType is "all": # returns all atoms found by recursing down the road
            # add rootAtm to the returnList
            returnList.append(rootAtm)
            # determine which of the bonded atoms can be added to the list
            for atm in bonded:
                for at in graphSearch(atm, mol, searchType="all"): returned.append(at)
        # append the returned atoms to the list to be returned
        for atm in returned: returnList.append(atm)
            
    return returnList
        
''' Method for all atom, bond etc data to be copied from one Molecule to another, 
with the possibility of excluding some '''
def copyData(sourceMol, destMol, joiningPos=-1, exclusions=[], onlyOne=False):
    # copy atoms
    for atm in sourceMol.getAtms():
        if atm.getAtmIndex() not in exclusions:
            destMol.addAtm(copy(atm), True)
    # copy bonds
    for bond in sourceMol.getBonds():
        if not onlyOne: # all indexs have to not be in exclusions
            if bond.getAtmAIndex() not in exclusions and bond.getAtmBIndex() not in exclusions:
                destMol.addBond(copy(bond), True)
        else: # only one index has to not be in exclusions
            if bond.getAtmAIndex() not in exclusions or bond.getAtmBIndex() not in exclusions:
                destMol.addBond(copy(bond), True)
    # copy angles
    for angle in sourceMol.getAngles():
        if not onlyOne:
            if angle.getAtmAIndex() not in exclusions and angle.getAtmBIndex() not in exclusions and angle.getAtmCIndex() not in exclusions:
                destMol.addAngle(copy(angle), True)
        else:
            if angle.getAtmAIndex() not in exclusions or angle.getAtmBIndex() not in exclusions or angle.getAtmCIndex() not in exclusions:
                destMol.addAngle(copy(angle), True)
    # copy impropers
    for imp in sourceMol.getImpropers():
        if not onlyOne:
            if imp.getAtmAIndex() not in exclusions and imp.getAtmBIndex() not in exclusions and imp.getAtmCIndex() not in exclusions and imp.getAtmDIndex() not in exclusions:
                destMol.addImproper(copy(imp), True)
        else:
            if imp.getAtmAIndex() not in exclusions or imp.getAtmBIndex() not in exclusions or imp.getAtmCIndex() not in exclusions or imp.getAtmDIndex() not in exclusions:
                destMol.addImproper(copy(imp), True)
    # copy dihedrals
    for dihed in sourceMol.getDihedrals():
        if not onlyOne:
            if dihed.getAtmAIndex() not in exclusions and dihed.getAtmBIndex() not in exclusions and dihed.getAtmCIndex() not in exclusions and dihed.getAtmDIndex() not in exclusions:
                destMol.addDihedral(copy(dihed), True)
        else:
            if dihed.getAtmAIndex() not in exclusions or dihed.getAtmBIndex() not in exclusions or dihed.getAtmCIndex() not in exclusions or dihed.getAtmDIndex() not in exclusions:
                destMol.addDihedral(copy(dihed), True)
    # copy joining point data
    for i in range(len(sourceMol.getJoiningPoints())):
        if i != joiningPos:
            destMol.addJoiningPoint(copy(sourceMol.getJoiningPoints()[i]), True)
            destMol.addTransBond(sourceMol.getTransBonds()[i])
    # return the destMol
    return destMol


''' Method to generate a change matrix.
It gives the index shifts between two molecules'''
def genChangeMatrix(sourceMol, destMol, additionalShift=0):
    m = {}
    shift=0
    for i in range(len(sourceMol.getAtms())): # index of molA
        if sourceMol.getAtm_Index(i).getAtmIndex() == destMol.getAtm_Index(shift).getAtmIndex():
            m[i+1] = shift+1+additionalShift # store as positions,(x,y) such that x becomes y
            if shift < len(destMol.getAtms())-1: shift += 1 # stop adding to shift once it's as big as can get so don't get IndexOutOfRange errors
        else: m[i+1] = -1     # store as positions, -1 as been deleted  
    return m

''' Method to shift all atmIndexs in a molecule based on a change matrix.'''
def shiftIndices(molecule, changeMatrix):
    # shift the atom index
    for a in molecule.getAtms():
        a.setAtmPos(changeMatrix[a.getAtmPos()])
        ## shift the exclusions
        exlus = a.getNonBondExcludeAtms()
        newExcludes = []
        for i in range(len(exlus)):
            try:
                if changeMatrix[exlus[i]] != -1 and changeMatrix[exlus[i]] > a.getAtmPos():
                    newExcludes.append(changeMatrix[exlus[i]])
            except KeyError:
                pass
        a.setNonBondExcludeAtms(newExcludes)
    # shift the atom indexes referenced in bonds
    for a in molecule.getBonds():
        a.setAtmAPos(changeMatrix[a.getAtmAPos()])
        a.setAtmBPos(changeMatrix[a.getAtmBPos()])
    # shift the atom indexes referenced in angles
    for a in molecule.getAngles():
        a.setAtmAPos(changeMatrix[a.getAtmAPos()])
        a.setAtmBPos(changeMatrix[a.getAtmBPos()])
        a.setAtmCPos(changeMatrix[a.getAtmCPos()])
    # shift the atom indexes referenced in impropers
    for a in molecule.getImpropers():
        a.setAtmAPos(changeMatrix[a.getAtmAPos()])
        a.setAtmBPos(changeMatrix[a.getAtmBPos()])
        a.setAtmCPos(changeMatrix[a.getAtmCPos()])
        a.setAtmDPos(changeMatrix[a.getAtmDPos()])
    # shift the atom indexes referenced in dihedrals
    for a in molecule.getDihedrals():
        a.setAtmAPos(changeMatrix[a.getAtmAPos()])
        a.setAtmBPos(changeMatrix[a.getAtmBPos()])
        a.setAtmCPos(changeMatrix[a.getAtmCPos()])
        a.setAtmDPos(changeMatrix[a.getAtmDPos()])
    # shift the atom indexes of joining points
    l = []
    for a in molecule.getJoiningPoints():
        l.append(changeMatrix[a+1]-1)
    molecule.setJoiningPoint(l)
    # shift the bondIDs of the transverse bonds
    newIDS = []
    for b in molecule.getTransBonds():
        bits = list(map(int,b.split(',')))
        bits[0] = changeMatrix[bits[0]+1]
        bits[1] = changeMatrix[bits[1]+1]
        if bits[0] < bits[1]:
            newIDS.append(str(bits[0]-1)+","+str(bits[1]-1))
        else:
            newIDS.append(str(bits[1]-1)+","+str(bits[0]-1))
    molecule.setTransBonds(newIDS)
    
        
    return molecule    

''' Method for handling the changes in shift matrix due to the addition of two molecules.'''
def handleJoin(chngMatA, chngMatB, molA, molB, diffA, diffB, AJoinPos, BJoinPos):
    # changMat is the change matrix as generated between mol and trunc
    # diff are the atmIndexs removed from mol to get trunc
    # the A series are the important ones
    # Handling of the joining atoms, such that the bond between the two fragments is handled
    # deal with the joining point atoms. Get position of not kept atom from bond across zone
    bondA = molA.getBond(molA.getBondIDs().index(molA.getTransBonds()[AJoinPos])) # bond between zones in molA
    if bondA.getAtmAIndex() in diffA:
        lostA = bondA.getAtmAPos()
        rootA = bondA.getAtmBPos()
    else:
        lostA = bondA.getAtmBPos()
        rootA = bondA.getAtmAPos()
    # get ID of kept atom from bond across zone in B
    bondB = molB.getBond(molB.getBondIDs().index(molB.getTransBonds()[BJoinPos])) # bond between zones in molB
    if bondB.getAtmAIndex() not in diffB:
        keptB = bondB.getAtmAPos()
        rootB = bondB.getAtmBPos()
    else:
        keptB = bondB.getAtmBPos()
        rootB = bondB.getAtmAPos()
    # set molA shift matrix for lost A to be shifted of keptB
    chngMatA[lostA] = chngMatB[keptB]
    
    # Handling of the atoms one bond away from the joining atoms. Should handle the angles between the fragments
    # make a list of the atoms one bond away from the joining atoms, for each molecule
    # find the bonds that involve the joining atoms, and not the root atom
    oneBondA = []
    for bond in molA.getBonds():
        if bond.getAtmAPos() == lostA and not bond.getAtmBPos() == rootA:
            oneBondA.append(bond.getAtmBPos())
        elif bond.getAtmBPos() == lostA and not bond.getAtmAPos() == rootA:
            oneBondA.append(bond.getAtmAPos())
    oneBondB = []
    for bond in molB.getBonds():
        if bond.getAtmAPos() == keptB and not bond.getAtmBPos() == rootB:
            oneBondB.append(bond.getAtmBPos())
        elif bond.getAtmBPos() == keptB and not bond.getAtmAPos() == rootB:
            oneBondB.append(bond.getAtmAPos())
    # compare the list of atoms
    assignedBs = [] # add to this list once an atom of A has been assigned a corresponding B so that we don't get double assingment
    for atomA in oneBondA:
        for atomB in oneBondB:
            if molA.getAtm_Pos(atomA).getAtmType() == molB.getAtm_Pos(atomB).getAtmType() and atomB not in assignedBs:
                # if also the same IACM, strongly assign
                if molA.getAtm_Pos(atomA).getIACM() == molB.getAtm_Pos(atomB).getIACM():
                    chngMatA[molA.getAtm_Pos(atomA).getAtmPos()] = chngMatB[molB.getAtm_Pos(atomB).getAtmPos()]
                    assignedBs.append(atomB)
                    break
                else: # otherwise weakly assign it. Needs to be implemented somehow
                    chngMatA[molA.getAtm_Pos(atomA).getAtmPos()] = chngMatB[molB.getAtm_Pos(atomB).getAtmPos()]
                    assignedBs.append(atomB)
                    break
    # Handling the atoms two bonds away form the joining atoms. Should finish the dihedrals off
    # make a list of all the atoms two bonds away from the joining atoms, ie one bond away from oneBondX
    # for each molecule
    assignedBs = []
    twoBondA = []
    twoBondB = []
    for bond in molA.getBonds():
        if bond.getAtmAPos() in oneBondA and not bond.getAtmBPos() == lostA:
            twoBondA.append(bond.getAtmBPos())
        elif bond.getAtmBPos() in oneBondA and not bond.getAtmAPos() == lostA:
            twoBondA.append(bond.getAtmAPos())
    for bond in molB.getBonds():
        if bond.getAtmAPos() in oneBondB and not bond.getAtmBPos() == keptB:
            twoBondB.append(bond.getAtmBPos())
        elif bond.getAtmBPos() in oneBondB and not bond.getAtmAPos() == keptB:
            twoBondB.append(bond.getAtmAPos())
    # compare the list of atoms
    for atomA in twoBondA:
        for atomB in twoBondB:
            if molA.getAtm_Pos(atomA).getAtmType() == molB.getAtm_Pos(atomB).getAtmType() and atomB not in assignedBs:
                # if also the same IACM, strongly assign
                if molA.getAtm_Pos(atomA).getIACM() == molB.getAtm_Pos(atomB).getIACM():
                    chngMatA[molA.getAtm_Pos(atomA).getAtmPos()] = chngMatB[molB.getAtm_Pos(atomB).getAtmPos()]
                    assignedBs.append(atomB)
                    break
                else: # otherwise weakly assign it. Needs to be implemented somehow
                    chngMatA[molA.getAtm_Pos(atomA).getAtmPos()] = chngMatB[molB.getAtm_Pos(atomB).getAtmPos()]
                    assignedBs.append(atomB)
                    break    
    
    return chngMatA

''' Method for handling the charge balancing in a joint molecule '''
def balanceCharge(jointMolecule, targetCharge):
    '''
    Current: adds/subtracts 0.001 to each atom in order from the largest magnitude to the least, until 
    the net charge reaches the target charge.    
    Expansion plans:
    3. Add the charge to the atoms in such a manner that there is a minimal change in the molecule's
    dipole moment. The big issue here is determining what the dipole should be, or how to obtain it 
    based on the two constiuent molecules. Also coming up with an elegant means of determining how the
    charge should be distributed. Brute force will work but be slow and ugly.
    4. Possible better means of handling the charge dispersion required.
    5. Bastardise the ATB's charge group code to redetermine the charge groups on a molecule. Possibly
    only do this on the final molecule that is generated. 
    '''
    chargeVariables = {"totalCharge":0,"deltaCharge":0.001,"chargeChangeCount":0,"targetCharge":targetCharge}
    # calculate the current total charge on the molecule
    for a in jointMolecule.getAtms():
        chargeVariables["totalCharge"] += a.getCharge()
    chargeVariables["totalCharge"] = round(chargeVariables["totalCharge"],3)
    iterCharge = abs(chargeVariables["totalCharge"] - chargeVariables["targetCharge"])
    
    # get the molecule's dipole
    moleculeDipole = jointMolecule.getDipole()
    
    # get a sorted list of the atoms, based on the charges
    atmList = sorted(jointMolecule.getAtms(), key=lambda atm: abs(atm.charge), reverse=True)
    
    # balance out the charge as needed
    print("--Target charge of {targetCharge} needs to be reached from a current charge of {totalCharge}. Delta is +/- {deltaCharge}.".format(**chargeVariables))    
    while abs(abs(chargeVariables["totalCharge"]) - abs(targetCharge)) > 0.00001:
        for atom in atmList:
            currentCharge = atom.getCharge()
            if chargeVariables["totalCharge"] < targetCharge:
                newCharge = round(currentCharge + chargeVariables["deltaCharge"],3)
                chargeVariables["totalCharge"] = round(chargeVariables["totalCharge"] + chargeVariables["deltaCharge"],3)
                atom.setCharge(newCharge)
                chargeVariables["chargeChangeCount"] += 1
            elif chargeVariables["totalCharge"] > targetCharge:
                newCharge = round(currentCharge - chargeVariables["deltaCharge"],3)
                chargeVariables["totalCharge"] = round(chargeVariables["totalCharge"] - chargeVariables["deltaCharge"],3)
                atom.setCharge(newCharge)
                chargeVariables["chargeChangeCount"] += 1
            elif -0.00001 < chargeVariables["totalCharge"] < 0.00001:
                break
    print("--{chargeChangeCount} charge assignments of {deltaCharge} were needed to obtain the target charge of {totalCharge}".format(**chargeVariables))
    # some information on the changed dipole
    newDipole = jointMolecule.getDipole()
    print("---Dipole delta is ",newDipole-moleculeDipole)
    print("---Dipole magnitude delta is ",numpy.linalg.norm(newDipole-moleculeDipole))
    # reassign charge groups based on the new charges
    jointMolecule = newChargeGroups(jointMolecule)
    ######################################
    ##### Sensible charge assignment #####
    ######################################
    # calculate number of charge groups
    
    ''' Method for performing least squares fitting of excess/deficient charges to the molecule in
    order to gain the required targetCharge. A matrix of linear equation coefficients is generated,
    along with a vector of the results, to give Ax = b, which is then solved for x using a least
    squares fit. The rules used for the linear equations are:
    1. The dipole of the overall molecule must not change. This means that the dipole of the 
    distributed charge must be (0,0,0). Leads to the first 3 linear equations: sum(x_i*c_i) = 0,
    sum(y_i*c_i) = 0, and sum(z_i*c_i) = 0. Coordinates are to be reference from the geometric centre
    of the molecule.
    2. The sum of the dispersed charges must equal the excess/deficient charge present on the molecule.
    Leads to the final obvious linear equation: sum(1*c_i) = requiredCharge.
    3. If the charge on an atom is small enough (-0.005 < x < 0.005), no more charge is to be added to it.
    Leads to 1*c_i = 0 for the atom, i, of interest.
    #### STILL TO IMPLEMENT RULE 4 ####
    4. Symmetry equivalent atoms should have the same new charge assigned to them so that they remain
    symmetry equivalent. For each pair of symmetry equivalent atoms, ij, and equation is added such that
    1*c_i + -1*c_j = 0
    5. Charge groups with 3 or more atoms in them should also have no overall dipole change. For each
    such charge group, this adds another 3 equations following rule 1 but for only those atoms in the group.
    Coordinates of the atoms are modified such that they are reference from the geometric centre of the
    charge group, not the geometric centre of the molecule.
    6. If the charge on a charge group is currently an integer, force that charge group to have zero total
    charge added.
    #### POSSIBLE NEW RULE ####
    7. Given that the excess charge come from an area around the joining point of the two molecules, only
    apply the delta to those atoms within a certain cutoff distance of the joining point.
    '''
def leastSquaresCharges(molecule, targetCharge):
    # calculate the current charge on the molecule
    currentCharge = 0.
    for atm in molecule.getAtms():
        currentCharge += atm.getCharge()
    # and use it to find the excess/deficient charge on the molecule
    requiredCharge = targetCharge - currentCharge

    print(requiredCharge)
    # build the matrix and vector based on the first two, always there rules
    coefficientMatrix = numpy.ones((1,molecule.getNumAtms())) # combined with generation of sumVector
    sumVector = numpy.array((requiredCharge))     # adds rule 2 in
    moleculeCentre = determineGeoCentre(molecule,list(range(molecule.getNumAtms())))
    if CHARGE_RULES["rule1"]:
        ##### RULE 1 START #####
        chargeMatrix = numpy.zeros((3,molecule.getNumAtms()))
        for i in range(molecule.getNumAtms()):
            chargeMatrix[0,i] = molecule.getAtms()[i].getX()-moleculeCentre[0]
            chargeMatrix[1,i] = molecule.getAtms()[i].getY()-moleculeCentre[1]
            chargeMatrix[2,i] = molecule.getAtms()[i].getZ()-moleculeCentre[2]
    # add the remaining rules in with atm based rules
    for atm in molecule.getAtms():
        if CHARGE_RULES["rule3"]:
            ##### RULE 3 START #####
            if -0.005 < atm.getCharge() < 0.005:
                tempVector = numpy.zeros((1,molecule.getNumAtms()))
                tempVector[0,atm.getAtmIndex()] = 1
                coefficientMatrix = numpy.vstack((coefficientMatrix,tempVector))
                sumVector = numpy.hstack((sumVector,numpy.array((0.))))
                print("---Atom %d has too low a charge. No excess charge will be assigned to it." % atm.getAtmPos())
            ##### RULE 3 END #####
    # add in charge group based rules
    for group in molecule.getChargeGroups():
        if CHARGE_RULES["rule5"]:
            ##### RULE 5 START #####
            if len(group) >= 3:
                # calculate the geo centre of the charge group
                geoCentre = determineGeoCentre(molecule, group)
                # make the matrix needed
                chargeMatrix = numpy.zeros((3,molecule.getNumAtms()))
                for atm in group:
                    chargeMatrix[0,atm] = molecule.getAtms()[atm].getX()-geoCentre[0]
                    chargeMatrix[1,atm] = molecule.getAtms()[atm].getY()-geoCentre[1]
                    chargeMatrix[2,atm] = molecule.getAtms()[atm].getZ()-geoCentre[2]
                # add the chargeMatrix to the coefficientMatrix and 0's to the sumVector
                coefficientMatrix = numpy.vstack((coefficientMatrix,chargeMatrix))
                sumVector = numpy.hstack((sumVector,numpy.array((0.,0.,0.))))
            ##### RULE 5 END #####
        if CHARGE_RULES["rule6"]:
            ##### RULE 6 START #####
            # determine charge group charge
            totalCharge = 0.
            for atm in group:
                totalCharge += molecule.getAtms()[atm].getCharge()
            if (-1e-10 < totalCharge < 1e-10) or (-1-1e-10 < totalCharge < -1+1e-10) or (1-1e-10 < totalCharge < 1+1e-10):
                chargeMatrix = numpy.zeros((1,molecule.getNumAtms()))
                for atm in group:
                    chargeMatrix[0,atm] = 1
                # add the chargeMatrix etc
                coefficientMatrix = numpy.vstack((coefficientMatrix,chargeMatrix))
                sumVector = numpy.hstack((sumVector,numpy.array((0.))))
            ##### RULE 6 END #####
    # perform the least squares fit
    results = numpy.linalg.lstsq(coefficientMatrix, sumVector, rcond=5.e-8)
    ''' The logic from here down is a bit odd and needs looking at '''
    # assign the charges
    rollOverCharge = 0.
    numDecPoint = 6
    for atm in molecule.getAtms():
        if atm.getAtmPos() == molecule.getNumAtms(): #last atom
            additionalCharge = round(results[0][atm.getAtmIndex()] + rollOverCharge,numDecPoint)
            atm.addCharge(additionalCharge)
            realTotalCharge = 0.
            for atom in molecule.getAtms():
                realTotalCharge += atom.getCharge()
            if requiredCharge-1e-10 < realTotalCharge < requiredCharge+1e-10: pass 
            else:
                realTotalCharge = round(realTotalCharge, numDecPoint)
                atm.setCharge(atm.getCharge() - realTotalCharge)
                print("----additional %f charge added to atmId %d" %(-realTotalCharge, atm.getAtmPos()))
        else:
            additionalCharge = round(results[0][atm.getAtmIndex()] + rollOverCharge,numDecPoint)
            atm.addCharge(additionalCharge)
            rollOverCharge = results[0][atm.getAtmIndex()] + rollOverCharge - additionalCharge
    ''' End of funky logic '''
    
''' Method for reorganising charge groups due to changing charges.
    Will use bastardised ATB code to do the assignments. '''
def newChargeGroups(mol):
    print("--Charge group reassignment not yet implemented")
    return mol

''' Method for determining the geometric centre of a group of atms in a molecule.
Pass it a molecule and a list of the atms in the group. Returns a vector of the
centre coordinates,'''
def determineGeoCentre(molecule, group):
    geoCentre = numpy.array((0.,0.,0.))
    for atm in group:
        geoCentre += molecule.getAtms()[atm].getXYZ(vec=True)
    geoCentre /= len(group)
    return geoCentre

''' Method for printing out an OVL file '''  
            
''' Method for running GROMOS stuff, mainly to do an energy minimisation '''
def runGROMOS(fileName, inputString):
    os.chdir(ROOT_DIRECTORY)
    # make a temporary folder for the molecule being run through
    os.system('mkdir GROMOSData/.temp')
    os.chdir('GROMOSData/.temp')
    # copy the root IFP and MTB files and PDB, and rename the PDB file to something sensible
    os.system('cp ../54A7_old.ifp ../54A7.mtb ../pdb2g96.lib "../../OutputFiles/%s/%s.pdb" ../energymin.imd ./' %(fileName, fileName))
    # concatanate the molecule onto the end of the 54a7 MTB in temporary file
    os.system('cat "../../OutputFiles/%s/%s.mtb" >> 54A7.mtb' %(fileName, fileName))
    # generate a topology file
    print("-Making the topology")
    topologyGenerate = 'make_top @build 54A7.mtb @param 54A7_old.ifp @seq %s @solv H2O > %s.top 2> /dev/null' %(fileName, fileName)
    os.system(topologyGenerate)
    # convert the PDB into a CNF
    os.system('pdb2g96 @topo %s.top @pdb %s.pdb @lib pdb2g96.lib > %s.cnf 2> /dev/null' %(fileName, fileName, fileName))
    # modify the energymin.imd file to account for the number of atoms in the molecule
    os.system("sed -i 's/NUMOFATOMS/%d/' energymin.imd" % LOADED_MOLECULES[fileName].getNumAtms())
    os.system("sed -i 's/MOLECULESHORTCODE/%s/' energymin.imd" % inputString)
    # run the energy minimisation
    print("-Running the energy minimisation for a maximum of 2500 steps")
    os.system('/usr/md++-1.2.3/bin/md @topo %s.top @conf %s.cnf @fin %s_min.cnf @input energymin.imd > energymin.omd' %(fileName, fileName, fileName))
    # frameout the minimised structure into PDB and mv it into the outputfiles
    os.system('frameout @topo %s.top @pbc v @outformat pdb @notimeblock @traj %s_min.cnf > /dev/null' %(fileName, fileName))
    print("-Saving the minimised structure as %s_minimised.pdb" % fileName)
    os.system("sed -i '1,2d' FRAME_00001.pdb")
    os.system('mv FRAME_00001.pdb "../../OutputFiles/%s/%s_minimised.pdb"' % (fileName, fileName))
    # cleanup the temp folder
    os.chdir('../')
    os.system('rm -r .temp')
    
''' Method for performing an exclusions check on a given molecule.
It makes a topology with the molecule as is, runs check_top on said topology,
And then proceeds to add in any missing exclusions that check_top picks up. '''
def exclusionsCheck(molecule):
    # print out the MTB and PDB to a temporary file
    os.chdir(ROOT_DIRECTORY)
    os.system('mkdir .temp')
    os.chdir('.temp')
    with open(molecule.getCompndName()+".pdb", "w") as fh:
        fh.write(molecule.printPDB())
    with open(molecule.getCompndName()+".mtb", "w") as fh:
        fh.write(molecule.printMTB())
    # copy the ifp, mtb and lib file into the temp folder
    os.system('cp ../GROMOSData/54A7.ifp ../GROMOSData/54A7.mtb ../GROMOSData/pdb2g96.lib  ./')
    os.system('cat '+molecule.getCompndName()+'.mtb >> 54A7.mtb')
    # make the topology file
    os.system('make_top @build 54A7.mtb @param 54A7.ifp @seq %s @solv H2O > temp.top 2>/dev/null' %molecule.getCompndName())
    # run check_top
    os.system('check_top @topo temp.top @build 54A7.mtb @param 54A7.ifp > checked.out 2>/dev/null')
    # read the checked file into memory
    with open('checked.out', "r") as fh:
        file = fh.readlines()
    # parse for the lines of interest to exclusions and update the molecule in memory
    for line in file:
        if "neighbour of atom" in line:
            splitLine = line.split()
            excludeAtm = int(splitLine[1])
            sourceAtm = int(splitLine[10])
            molecule.getAtm_Pos(sourceAtm).addNonBondExlude(excludeAtm)
    # clean up
    os.chdir('../')
    os .system('rm -r .temp')
    
''' Method for parsing through an inputstring to remove all the special formatting.
Not the most elegant, nor efficient, but it works so who cares? '''
def stringSanitation(inputString):
    sanatisedString = inputString
    # firstly need to remove any {x}{cccc} things. They will always be in that order
    while True:
        if "{" in sanatisedString:
            charsToRemove = []
            firstPos = sanatisedString.find("{")
            secondPos = sanatisedString.find("{", firstPos+1)
            charsToRemove.append(secondPos)
            for i in range(firstPos,secondPos):
                charsToRemove.append(i)
            tempString = ""
            for i in range(len(sanatisedString)):
                if i not in charsToRemove:
                    tempString += sanatisedString[i]
            sanatisedString = tempString
        elif "<" in sanatisedString:
            firstPos = sanatisedString.find("<")
            secondPos = sanatisedString.find(">")
            tempString = ""
            for i in range(len(sanatisedString)):
                if i not in (firstPos, secondPos):
                    tempString += sanatisedString[i]
            sanatisedString = tempString
        elif "(" in sanatisedString:
            firstPos = sanatisedString.find("(")
            secondPos = sanatisedString.find(")")
            tempString = ""
            for i in range(len(sanatisedString)):
                if i not in (firstPos, secondPos):
                    tempString += sanatisedString[i]
            sanatisedString = tempString
        elif "}" in sanatisedString:
            firstPos = sanatisedString.find("}")
            tempString = ""
            for i in range(len(sanatisedString)):
                if i != firstPos:
                    tempString += sanatisedString[i]
            sanatisedString = tempString
        else:
            break
    
    return sanatisedString
    
''' Method to decide if a given string is a sensible item, ie not containing other excess things '''
def isJoinable(questionString):
    for eachChar in STRING_BREAK_CHARACTERS:
        if eachChar in questionString:
            return False
    return True

''' Method to determine if rotation is needed for two given molecules.
Rotation is needed if there is overlap between the two.'''
def rotateNeeded(molA, molB):
    for atmOne in molA.getAtms():
        for atmTwo in molB.getAtms():
            atmOneType = atmOne.getAtmType()
            atmTwoType = atmTwo.getAtmType()
            minLength = 0.
            for atm in (atmOneType,atmTwoType):
                if atm not in ATOMIC_RADII:
                    minLength += 1.0
                else:
                    minLength += ATOMIC_RADII[atm]
            actualLength = numpy.linalg.norm(atmOne.getXYZ(vec=True)-atmTwo.getXYZ(vec=True))
            if actualLength < minLength:
                return True, atmOne.getAtmIndex(), atmTwo.getAtmIndex(), minLength - actualLength
            else:
                return False, 0, 0, 0

''' Method to return all the connectivity involving a given atmID, with that atmID in a terminal position.
Returns all bonds, angles and dihedrals (regardless of if they're assigned in the force field or not) by
default, or just one type of connectivity'''
def atomConnectivity(atmID, mol, kind="all"):
    # search for all bond involving the atmID
    bonds = []
    for bond in mol.getBonds():
        if bond.getAtmAIndex() == atmID: bonds.append([atmID,bond.getAtmBIndex()])
        elif bond.getAtmBIndex() == atmID: bonds.append([atmID,bond.getAtmAIndex()])
    if kind == "bonds":
        return bonds
    angles = []
    for source in bonds:
        for bond in mol.getBonds():
            if bond.getAtmAIndex() == source[1] and bond.getAtmBIndex() not in source: angles.append([atmID,source[1],bond.getAtmBIndex()])
            elif bond.getAtmBIndex() == source[1] and bond.getAtmAIndex() not in source: angles.append([atmID,source[1],bond.getAtmAIndex()])
    if kind == "angles":
        return angles
    dihedrals = []
    for source in angles:
        for bond in mol.getBonds():
            if bond.getAtmAIndex() == source[2] and bond.getAtmBIndex() not in source: dihedrals.append([atmID,source[1],source[2],bond.getAtmBIndex()])
            elif bond.getAtmBIndex() == source[2] and bond.getAtmAIndex() not in source: dihedrals.append([atmID,source[1],source[2],bond.getAtmAIndex()])
    if kind == "dihedrals":
        return dihedrals
    return (bonds,angles,dihedrals)

''' Method for generating a graph representation of a give molecule. '''
def genGraphRep(mol):
    graph = {}
    for atm in mol.getAtms():
        graph[atm.getAtmIndex()] = []
    for bond in mol.getBonds():
        graph[bond.getAtmAIndex()].append(bond.getAtmBIndex())
        graph[bond.getAtmBIndex()].append(bond.getAtmAIndex())
    return graph

def find_all_paths(graph, start, end, path=[]):
    path = path + [start]
    if start == end: return [path]
    if start not in graph: return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end, path)
            for newpath in newpaths: paths.append(newpath)
    return paths

def find_all_paths_of_length(graph, length):
    all_paths = []
    for sourcenode in graph:
        for destnode in graph:
            start = min((sourcenode,destnode))
            end = max((sourcenode,destnode))
            paths = find_all_paths(graph, start, end)
            for path in paths: 
                if len(path) == length and path not in all_paths: all_paths.append(path)
    return all_paths
    
if __name__ == "__main__":
    words = stringSanitation("AAAA{5}{BBBB{17}{12fr}CCCC(ZZZZ<MGHTISSK>(1234)AAAA)}")
    print(words)