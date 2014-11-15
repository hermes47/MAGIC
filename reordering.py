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
#mol = Molecule()
#rootName = "BEHM"
#genmethods.parsePDB(rootName, mol)
#genmethods.parseMTB(rootName, mol)
#shiftMatrix = {1:-1, 2:-2, 3:-3, 13:4, 10:5, 11:6, 12:7, 4:8, 6:9, 5:10, 7:11, 9:12, 8:13, 14:14, 15:15, 16:16, 17:17, 26:18, 23:19, 22:20, 21:21, 24:22, 25:23, 19:24, 18:25, 27:26, 28:27, 20:28, 29:29, 30:30}
shiftMatrix2 = {}
names = {23:"TCSM",24:"LIGM",25:"PCSM",26:"CERM",27:"HEPM",28:"MONM",
         29:"NCSM",30:"MELM",31:"HTCM",32:"LACM",33:"PSYM",34:"GEDM",35:"CRPM",36:"HTYM"}
for i in [23,24,25,26,27,28,29,30,31,32,33,34,35,36]:
    mol = Molecule()
    rootName = "HTYM"
    genmethods.parsePDB(rootName, mol)
    genmethods.parseMTB(rootName, mol)
    for j in range(1,i+6):
        shiftMatrix2[j] = j
    for k in range(1,100):
        if k not in shiftMatrix2: shiftMatrix2[k] = -1
    print(shiftMatrix2)
    shiftIndicesNotAtm(mol, shiftMatrix2)
    mol.setCompndName(names[i])
    mol.getAtmWithIndex(mol.getNumAtms()-1).setIACM(16)
    mol.getAtmWithIndex(mol.getNumAtms()-1).setMassNum(5)
    #copyMol = copy.deepcopy(mol)
    #copyMol.atms = []
    #for i in range(len(mol.getAtms())):
    #    copyMol.addAtm(mol.getAtms()[shiftMatrix[i+1]-1],True)
    #    copyMol.getAtms()[i].setAtmIndex(i)
    #    copyMol.getAtms()[i].setResName(rootName)
    #genmethods.exclusionsCheck(copyMol)
    os.chdir(ROOT_DIRECTORY)
    output = names[i]
    try:
        os.system('mkdir OutputFiles/'+output)
    except:
        pass
    # print out final file to PDB and MTB files
    with open("OutputFiles/"+output+"/"+output+".pdb", "w") as fh:
        fh.write(mol.printPDB())
    with open("OutputFiles/"+output+"/"+output+".mtb", "w") as fh:
        fh.write(mol.printMTB())
    with open("OutputFiles/"+output+"/"+output+".ovl", "w") as fh:
        fh.write("4\n4,6\n")
        fh.write(str(i-2)+" "+" ".join(str(x) for x in range(7,6+i)))
        fh.write("\n1 6\n2 4 5\n2 3 2\n1 1")
        
    
