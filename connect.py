'''
---- Basic Use Case ----
String passed:  AAAABBBBCCCCDDDD
Requirements:   MTB file leaders of NAME_LENGTH alpha characters in length
Outcome:        Connect the given MTB files together in a linear fashion. The first joining point of AAAA is
                connected to the first joining point of BBBB. Then the second joining point of BBBB is connected
                to the first of CCCC, and so on. 
---- Connect Enclosed to Root ----
String passed:  AAAA(BBBBCCCCDDDD)
Requirements:   ( ) surrounding the MTB's to sequentially connect
Outcome:        Connect the MTB's within the ( and ) to the immediately preceeding MTB, in sequential order. AAAA
                will have BBBB connected to its first joining point, CCCC to its second and DDDD to its third. A
                new NAME_LENGTH character string will be created for the combined MTB's
---- Connect Multiples ----
String passed:  AAAA{#}{BBBBCCCC}
Requirements:   {#}{ } where # is the number of repeats to perform and the second { } set contains the MTB file
                strings to connect
Outcome:        Return the expanded string of linear connection that is required. If # = 3, returned string will 
                be AAAABBBBCCCCBBBBCCCCBBBBCCCC
---- Connect Enclosed in ring ----
String passed:  AAAA<BBBBCCCCDDDDEEEE>
Requirements:   < > surrounding the MTB's to form a ring out of
Outcome:        Connect the MTB's within the < and > such that a ring is formed. Joining is performed linearly,
                except for the first and final elements. The second joining point of BBBB is connected to the first
                of CCCC, etc; and the second joining point of EEEE is connected to the first of BBBB
'''
from config import NAME_LENGTH, LOADED_MOLECULES, MOLECULE_PREALIGNMENT
from molecule import Molecule
import genmethods, rotations
import numpy as np


''' initStringParse.
Parse initialString to pull out any special formatting that is present, and begin building up the MTB file
requested by initialString.
Working on the reverse string means that recursion is started at the lowest level, and then built up. It 
felt easier to code this way.
Firstly, account for any multiplication of substrings that needs to occur. This formatting is just for ease
of use, for example in generating a polymer. As such, all such formatting should be expanded out to a full
form before further parsing is undertaken. Returns the expanded string.
Second, parse for any root-branch connections that need to be made, denoted by ( ). Generate the MTB file
requested by the ****( ) string obtained. Give it a new, unique and random name, which is returned,
replacing the ****( ) in the initialString. Further recursion occurs until no more ( ) are present.
Thirdly, parse for any ring formating that needs to occur, denoted by < >/ Generate the MTB file of the
ring requested by the < * > string. Give it a new, unique and random name, which is returned, replacing the
< * > in the initialString. Further recursion occurs until no more < > are present.
The final string returned is one which conforms to the requirements of ---- Basic Use Case ----
'''
def initStringParse(initialString):
    # reverse the initialString, so that parsing it works easier
    reversedString = initialString[::-1]
    # if the string is purely containing alpha, call the MTB connect function on it
    if genmethods.isJoinable(reversedString):
        if len(reversedString) == NAME_LENGTH:
            return initialString
        print("Performing basicJoin on: "+initialString)
        return basicJoin(initialString)
    elif "{" in reversedString: # check for any multiplication of substrings needed
        # partition the string at the first located {
        beforeSep, _, afterSep = reversedString.partition("{")
        # beforeSep will contain the tail of the string, a } and the string to multiply. Is re-reversed incase of nested multiplies
        multString, _, tail = beforeSep[::-1].partition("}")
        # afterSep will contain the reversed number of multiples after a }, a {, and the reversedHead of the string
        numMultiples, _, reversedHead = afterSep.partition("{")
        # if the string to multiply contains any ring or root-branch strings, parse them first
        # implemented to speed up overall final MTB generation
        if "(" in multString or "<" in multString:
            multString = initStringParse(multString)
        # call the multiplyString method
        reversedBody = multiplyString(int(numMultiples[:0:-1]), multString)[::-1] # numMultiples is reversed and the leading/trailing } is discarded
        # combine the strings in the correct order
        combineStrings = [tail[::-1],reversedBody,reversedHead]
        reversedString = "".join(s for s in combineStrings)
        # recurse to take care of nesting and stuff
        return initStringParse(reversedString[::-1])
    elif "(" in reversedString: # check for any enclosed to root sections
        # partition the string at the first located (
        beforeSep, _, afterSep = reversedString.partition("(")
        # beforeSep contains the tail of the string, a ), and the string of MTBs to be added to the backBone. Is re-reversed in case of nesting
        enclString, _, tail = beforeSep[::-1].partition(")")
        # if enclString contains < recurse on it to make the loop(s)
        if "<" in enclString:
            enclString = initStringParse(enclString)
        # afterSep will contain the reversedRoot as the first NAME_LENGTH characters, and the reversedHead as the rest of the string
        reversedRoot = afterSep[0:NAME_LENGTH]
        reversedHead = afterSep[NAME_LENGTH:]
        # run the combine protocol on the Root and encl strings
        print("Performing rootBranchJoining with: "+reversedRoot[::-1]+" as root and: "+enclString+" as the branches.")
        jointString = rootBranchJoin(reversedRoot[::-1], enclString)
        # combine the strings in the correct order
        combineStrings = [tail[::-1], jointString[::-1], reversedHead]
        reversedString = "".join(s for s in combineStrings)
        # recurse if there are more (
        return initStringParse(reversedString[::-1])
    elif "<" in reversedString: # check for any rings that need to be made
        # partition the string at the first located <
        beforeSep, _, afterSep = reversedString.partition("<")
        # beforeSep contains the tail of the string, a >, and the string of MTBs to be looped. Is re-reversed in case of nesting
        ringString, _, tail = beforeSep[::-1].partition(">")
        # afterSep contains the reversedHead of the string
        reversedHead = afterSep
        # run the combine protocol on the ringString
        jointString = ringJoin(ringString)
        # combine the strings in the correct order
        combineStrings = [tail[::-1],jointString[::-1],reversedHead]
        reversedString = "".join(s for s in combineStrings)
        # recurse
        return initStringParse(reversedString[::-1])

''' Multiply out multString, numMults times. Return the resultant string '''    
def multiplyString(numMults, multString):
    print("Multiplying out the string: "+multString+", "+str(numMults)+" times.")
    tempString = ""
    i = 0
    while i < numMults:
        tempString += multString
        i += 1
    return tempString


        
''' Perform a root-branch joining between root and all the branches in branchString.
    Splits branchString into its member branches, each of NAME_LENGTH characters in length.
    Loads all the required molecules, performs the joining, and saves the final molecule as
    a new file(s), with a randomly generated name. Returns the generate random name '''
def rootBranchJoin(root, branchString):
    newName = root
    # split the branchString into its member branch names of NAME_LENGTH length
    branches = genmethods.fragmentString(branchString)
    # perform sanity check to make sure root has enough joining points for all the branches
    if LOADED_MOLECULES[root].getNumJoiningPoints() >= len(branches):
        for bra in branches:
            ## determine the long axis of the molecule
            #longAxis = rotations.longAxis(LOADED_MOLECULES[bra])
            ## determine the target axis of the molecule, ie the bond that is being rotated on to
            #atmA, atmB = LOADED_MOLECULES[newName].getTransBonds()[0].split(",")
            #targetAxis = LOADED_MOLECULES[newName].getAtms()[int(atmA)].getXYZ(vec=True) - LOADED_MOLECULES[newName].getAtms()[int(atmB)].getXYZ(vec=True)
            ## determine the angle between the longAxis and the target axis
            #theta = np.arccos(np.dot(longAxis, targetAxis)/(np.linalg.norm(longAxis)*np.linalg.norm(targetAxis)))
            ## determine the axis of rotation, given by the cross product of the target and long axis
            #rotationAxis = np.cross(longAxis,targetAxis)
            ## normalise the rotation vector
            #rotationAxis = rotationAxis/np.linalg.norm(rotationAxis)
            ## rotate the molecule
            #rotations.rotateMolecule(LOADED_MOLECULES[bra], rotationAxis, 0)
            newName = simpleJoin(newName, bra)
        # perform a check of the 1,2 and 1,3 exclusions using check_top
        print("Performing a check to ensure all 1,2 and 1,3 exclusions are present in final molecule")
        genmethods.exclusionsCheck(LOADED_MOLECULES[newName])
    else:
        print(root+" has insignificant joining points for the requested system: "+root+"("+branchString+")")
    return newName

''' Perform a basic use case joining of all the molecules in moleculeString into a single
    new molecule. Splits moleculeString into its member molecules, each of NAME_LENGTH
    characters in length. Loads all the required molecules, performs the joining, and saves
    the final molecule as a new file, with a randomly generated name. Returns the generate
    random name '''
def basicJoin(moleculeString):
    print("Basic Use Case Joining not yet fully implemented")
    newName = genmethods.genRandomString()
    # split the moleculeString into its member molecule names of NAME_LENGTH length
    molecules = genmethods.fragmentString(moleculeString)
    return newName

''' Perform a simple join of two molecules into a single.
Takes the two molecules as arguments and joins them,
connecting moleculeA's first joining point to moleculeB's first joiningPoint,
or as requested by AJoinPos and BJoinPos.
Adds the molecule to the LOADED_MOLECULES list and returns the key. '''
def simpleJoin(moleculeA, moleculeB, AJoinPos = 0, BJoinPos = 0):
    jointName = genmethods.genRandomString()
    print("Joining "+moleculeA+" with "+moleculeB+" to make "+jointName)
    # size of overlapping portions. Larger overlapping portion is given position A, ie
    # data is kept for the angles, dihedrals etc. If even, first molecule is given A
    # based on joining points determine which atoms need to be deleted from each molecule
    '''
    Current: molA, the molecule that has it's parameters kept in the overlap, is determined to be the
    molecule that has the largest number of atoms deleted.
    Expansion plans:
    Have the overlap not be based just on the number of atoms removed. Step one would be to ignore H,
    as H is really nothing much of interest. Step two would be to have the longest chain length, again
    ignoring H, of the deleted atoms from within the molecule be the deciding factor. Final, and at 
    this stage best, implementation would be to determine the length of the overlapping chain such
    that only atoms that actually have equivalent atoms in the underlying molecule are accounted for.
    '''
    print("-Determining which atoms are to be removed from each molecule")
    abSwitch = False # so can keep track of which is the actual root, ie the original A
    sizeA = genmethods.graphSearch(LOADED_MOLECULES[moleculeA].getJoiningPoints()[AJoinPos], LOADED_MOLECULES[moleculeA], resetChecked=True)
    sizeB = genmethods.graphSearch(LOADED_MOLECULES[moleculeB].getJoiningPoints()[BJoinPos], LOADED_MOLECULES[moleculeB], resetChecked=True)
    if len(sizeA) >= len(sizeB):
        molA = LOADED_MOLECULES[moleculeA]
        molB = LOADED_MOLECULES[moleculeB]
        joinA = molA.getJoiningPoints()[AJoinPos]
        joinB = molB.getJoiningPoints()[BJoinPos]
        atmsToDeleteA = sizeA
        atmsToDeleteB = sizeB
    else:
        print("--moleculeA and moleculeB have been switched")
        abSwitch = True
        molB = LOADED_MOLECULES[moleculeA]
        molA = LOADED_MOLECULES[moleculeB]
        joinA = molA.getJoiningPoints()[AJoinPos]
        joinB = molB.getJoiningPoints()[BJoinPos]
        atmsToDeleteB = sizeA
        atmsToDeleteA = sizeB
    #################################################
    ### A and B may have swapped after this point ###
    #################################################
    # centre each molecule on its joining point
    print("-Centering "+molA.getCompndName()+" on atom index "+str(joinA))
    molA.centreOn(molA.getJoiningPoints()[AJoinPos])
    print("-Centering "+molB.getCompndName()+" on atom index "+str(joinB))
    molB.centreOn(molB.getJoiningPoints()[BJoinPos])
    # align the molecules correctly
    # required bond vector for molA is the vector between the kept atom and the deleted atom immediately joint to it
    # get the required bond vectors, calculate the angle between them, get the cross product of the vectors, rotate
    if MOLECULE_PREALIGNMENT:
        print("-Aligning molecules")
        tupA = molA.getTransBondsTuple(AJoinPos)
        tupB = molB.getTransBondsTuple(BJoinPos)
        if tupA[0] not in atmsToDeleteA:
            vecA = molA.getAtm_Index(tupA[0]).getXYZ(True) - molA.getAtm_Index(tupA[1]).getXYZ(True)
        else:
            vecA = molA.getAtm_Index(tupA[1]).getXYZ(True) - molA.getAtm_Index(tupA[0]).getXYZ(True)
        if tupB[0] in atmsToDeleteB:
            vecB = molB.getAtm_Index(tupB[0]).getXYZ(True) - molB.getAtm_Index(tupB[1]).getXYZ(True)
        else:
            vecB = molB.getAtm_Index(tupB[1]).getXYZ(True) - molB.getAtm_Index(tupB[0]).getXYZ(True)
        rotAxis = rotations.rotationVectorDetermine(vecA, vecB)
        rotAngle = rotations.rotationAngleDetermine(vecA, vecB)
        if abSwitch:
            rotations.rotateMolecule(molB, molA, vect=rotAxis, ang=rotAngle)
        else:
            rotations.rotateMolecule(molA, molB, vect=rotAxis, ang=rotAngle)
    # copy the data from molA to jointMol, keeping all bonds etc with at least one atom not deleted
    print("-Copying data")
    jointMol = genmethods.copyData(molA, Molecule(), AJoinPos, atmsToDeleteA, True)
    # copy the data from molB to tempMol, keeping only bonds etc that only have not deleted atoms
    tempMol = genmethods.copyData(molB, Molecule(), BJoinPos, atmsToDeleteB)
    # determine a shift matrix for the changing atmPos
    print("-Generating change matrices")
    changeMatrixA = genmethods.genChangeMatrix(molA, jointMol)
    changeMatrixB = genmethods.genChangeMatrix(molB, tempMol, jointMol.getNumAtms())
    # apply the shift to tempMol
    print("-Applying the change matrix to tempMol")
    tempMol = genmethods.shiftIndices(tempMol, changeMatrixB)
    # handle the joining of the two molecules by modifying changeMatrixA to account for the new indices that
    # come from molB being added on.
    # Only need to worry about changeMatrixA as only A bonds etc are copied without
    print("-Handling the joining of the two molecules")
    changeMatrixA = genmethods.handleJoin(changeMatrixA, changeMatrixB, molA, molB, atmsToDeleteA, atmsToDeleteB, AJoinPos, BJoinPos)
    # apply the shift matrix to jointMol
    print("-Applying the shift matrix to jointMol")
    jointMol = genmethods.shiftIndices(jointMol, changeMatrixA)
    # determine if there is need of a rotation for the original B molecule, based on overlapping
    print("-Determining if further molecule rotation is needed...")
    rotNeedInfos = genmethods.rotateNeeded(jointMol, tempMol)
    #rotNeedInfos = (False,False)
    if rotNeedInfos[0]:
        print("--YES")
        rotCount = 0
        while rotNeedInfos[0]:
            if abSwitch:
                rotations.rotateMolecule(jointMol, tempMol, rotMolAtmIndex=rotNeedInfos[1], refMolAtmIndex=rotNeedInfos[2], displacement=rotNeedInfos[3]-rotNeedInfos[4])
                rotCount += 1
            elif not abSwitch:
                rotations.rotateMolecule(tempMol, jointMol, rotMolAtmIndex=rotNeedInfos[2], refMolAtmIndex=rotNeedInfos[1], displacement=rotNeedInfos[3]-rotNeedInfos[4])
                rotCount += 1
            rotNeedInfos = genmethods.rotateNeeded(jointMol, tempMol)
        print("---%d rotations were performed to obtain non-overlapping structures" % rotCount)
    else:
        print("--NO")
        
    # join tempMol onto jointMol
    print("-Copying tempMol into jointMol to give the joint molecule")
    jointMol = genmethods.copyData(tempMol, jointMol)
    # balance out the charges
    totalCharge = 0.
    for a in molA.getAtms():
        totalCharge += a.getCharge()
    for a in molB.getAtms():
        totalCharge += a.getCharge()
    # round to 3 dp to account for inaccurate storing of floats
    totalCharge = round(totalCharge, 3)
    print("-Balancing the charges to obtain a target charge of %1.1f" % totalCharge)
    #genmethods.sensibleCharges(jointMol, totalCharge)
    #genmethods.balanceCharge(jointMol, totalCharge)
    genmethods.leastSquaresCharges(jointMol, totalCharge)
    jointMol.setCompndName(jointName)
    for a in jointMol.getAtms():
        a.setResName(jointName)
    LOADED_MOLECULES[jointName] = jointMol
    # print out the 3 files
    #print(LOADED_MOLECULES[jointName].getJoiningPoints())
    #print(LOADED_MOLECULES[jointName].getPDBStart())
    #for a in LOADED_MOLECULES[jointName].getAtms():
    #    print(a.getPDBString())
    #for a in LOADED_MOLECULES[jointName].getBonds():
    #    if a.getAtmAIndex() > a.getAtmBIndex():
    #        temp = a.getAtmBIndex()
    #        a.setAtmBIndex(a.getAtmAIndex())
    #        a.setAtmAIndex(temp)
    #    print(a.getPDBString())
    #print("END")
    return jointName


''' Perform a ring joining for the given string of molecules.
'''
def ringJoin(ringString):
    print("Ring joining not yet fully implemented")
    newName = genmethods.genRandomString()
    # split the ringString into its member molecule names of NAME_LENGTH length
    molecules = genmethods.fragmentString(ringString)
    return newName

if __name__ == "__main__":
    mol1 = Molecule()
    genmethods.parsePDB("P8EO", mol1)
    genmethods.parseMTB("P8EO", mol1)
    genmethods.parseOVL("P8EO", mol1)
    mol2 = Molecule()
    genmethods.parsePDB("B16H", mol2)
    genmethods.parseMTB("B16H", mol2)
    genmethods.parseOVL("B16H", mol2)
    mol3 = Molecule()
    genmethods.parsePDB("MCW5", mol3)
    genmethods.parseMTB("MCW5", mol3)
    genmethods.parseOVL("MCW5", mol3)
    LOADED_MOLECULES["P8EO"] = mol1
    LOADED_MOLECULES["B16H"] = mol2
    LOADED_MOLECULES["MCW5"] = mol3
    rootBranchJoin("P8EO", "B16H")