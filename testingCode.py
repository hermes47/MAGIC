import genmethods, molecule, copy, numpy, rotations

mol = molecule.Molecule()
genmethods.parsePDB("../../OutputFiles/Test/Exploded", mol)
genmethods.parseMTB("FRFY", mol)
tempMol = copy.deepcopy(mol)
scaleFactor = 1/2.5

for atm in tempMol.getAtms():
    newCoords = atm.getXYZ(vec=True)*scaleFactor
    atm.setXYZ(newCoords, vec=True)
    
with open("OutputFiles/TEST/Deploded.pdb","w") as fh:
    fh.write(tempMol.printPDB())
    
print("Comparing old parameters with new.\nBonds:")
for bond in mol.getBonds():
    atmA = bond.getAtmAIndex()
    atmB = bond.getAtmBIndex()
    distOld = numpy.linalg.norm(mol.getAtmWithIndex(atmA).getXYZ(vec=True) - mol.getAtmWithIndex(atmB).getXYZ(vec=True))
    distNew = numpy.linalg.norm(tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) - tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
    oldToNew = distOld/distNew
    print("%2.4f   %2.4f   %1.1f"%(distOld,distNew,oldToNew))
print("\nAngles:")
for angle in mol.getAngles():
    atmA = angle.getAtmAIndex()
    atmB = angle.getAtmBIndex()
    atmC = angle.getAtmCIndex()
    vect1Old = mol.getAtmWithIndex(atmA).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
    vect2Old = mol.getAtmWithIndex(atmC).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
    angOld = rotations.rotationAngleDetermine(vect1Old, vect2Old)*180.0/numpy.pi
    vect1New = (tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
    vect2New = (tempMol.getAtmWithIndex(atmC).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
    angNew = rotations.rotationAngleDetermine(vect1New, vect2New)*180.0/numpy.pi
    print("%3.4f    %3.4f    %1.1f"%(angOld,angNew,angOld/angNew))
print("\nDihedrals:")
for dihedral in mol.getDihedrals():
    atmD = dihedral.getAtmDIndex()
    atmA = dihedral.getAtmAIndex()
    atmB = dihedral.getAtmBIndex()
    atmC = dihedral.getAtmCIndex()
    vect1Old = mol.getAtmWithIndex(atmA).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
    vect2Old = mol.getAtmWithIndex(atmD).getXYZ(vec=True) -  mol.getAtmWithIndex(atmC).getXYZ(vec=True)
    angOld = rotations.rotationAngleDetermine(vect1Old, vect2Old)*180.0/numpy.pi
    vect1New = (tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
    vect2New = (tempMol.getAtmWithIndex(atmD).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmC).getXYZ(vec=True))
    angNew = rotations.rotationAngleDetermine(vect1New, vect2New)*180.0/numpy.pi
    print("%3.4f    %3.4f    %1.1f"%(angOld,angNew,angOld/angNew))
print("\nImpropers:")
for dihedral in mol.getImpropers():
    atmD = dihedral.getAtmDIndex()
    atmA = dihedral.getAtmAIndex()
    atmB = dihedral.getAtmBIndex()
    atmC = dihedral.getAtmCIndex()
    vect1Old = mol.getAtmWithIndex(atmA).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
    vect2Old = mol.getAtmWithIndex(atmD).getXYZ(vec=True) -  mol.getAtmWithIndex(atmC).getXYZ(vec=True)
    angOld = rotations.rotationAngleDetermine(vect1Old, vect2Old)*180.0/numpy.pi
    vect1New = (tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
    vect2New = (tempMol.getAtmWithIndex(atmD).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmC).getXYZ(vec=True))
    angNew = rotations.rotationAngleDetermine(vect1New, vect2New)*180.0/numpy.pi
    print("%3.4f    %3.4f    %1.1f"%(angOld,angNew,angOld/angNew))
