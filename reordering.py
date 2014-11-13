from molecule import Molecule
import genmethods, sys, os
from config import ROOT_DIRECTORY
import copy


def shiftIndicesNotAtm(molecule, changeMatrix):
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
    # shift the atom index
    for a in molecule.getAtms():
        a.setAtmPos(changeMatrix[a.getAtmPos()])
        ## shift the exclusions
        exlus = a.getNonBondExcludeAtms()
        newExcludes = []
        for i in range(len(exlus)):
            try:
                #if changeMatrix[exlus[i]] != -1 and changeMatrix[exlus[i]] > a.getAtmPos():
                if changeMatrix[exlus[i]] != -1:
                    newExcludes.append(changeMatrix[exlus[i]])
            except KeyError:
                pass
        a.setNonBondExcludeAtms(newExcludes)
    
    
# load up the Molecule
mol = Molecule()
rootName = "STEM"
genmethods.parsePDB(rootName, mol)
genmethods.parseMTB(rootName, mol)
#shiftMatrix = {1:-1, 2:-2, 3:-3, 13:4, 10:5, 11:6, 12:7, 4:8, 6:9, 5:10, 7:11, 9:12, 8:13, 14:14, 15:15, 16:16, 17:17, 26:18, 23:19, 22:20, 21:21, 24:22, 25:23, 19:24, 18:25, 27:26, 28:27, 20:28, 29:29, 30:30}
shiftMatrix2 = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1, 6:-1, 7:-1, 8:-1, 9:-1, 10:-1,
                11:-1, 12:1, 13:2, 14:3, 15:4, 16:5, 17:6, 18:7, 19:8, 20:9,
                21:10, 22:11, 23:12, 24:13, 25:14, 26:15, 27:16, 28:17, 29:18, 30:19,
                31:20, 32:21, 33:22, 34:-1, 35:-1, 36:-1, 37:-1, 38:-1, 39:-1, 40:-1,
                41:-1, 42:-1, 43:-1, 44:-1, 45:-1, 46:-1, 47:-1, 48:-1, 49:-1, 50:-1,
                51:-1, 52:-1}
shiftIndicesNotAtm(mol, shiftMatrix2)
mol.setCompndName(rootName)
#copyMol = copy.deepcopy(mol)
#copyMol.atms = []
#for i in range(len(mol.getAtms())):
#    copyMol.addAtm(mol.getAtms()[shiftMatrix[i+1]-1],True)
#    copyMol.getAtms()[i].setAtmIndex(i)
#    copyMol.getAtms()[i].setResName(rootName)
#genmethods.exclusionsCheck(copyMol)
os.chdir(ROOT_DIRECTORY)
output = rootName
#os.system('mkdir OutputFiles/'+output)
# print out final file to PDB and MTB files
with open("OutputFiles/"+output+"/"+output+".pdb", "w") as fh:
    fh.write(mol.printPDB())
with open("OutputFiles/"+output+"/"+output+".mtb", "w") as fh:
    fh.write(mol.printMTB())
    

