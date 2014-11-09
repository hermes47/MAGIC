import genmethods, molecule, os, numpy, rotations, copy

print(os.getcwd())
mol = molecule.Molecule()
genmethods.parsePDB("FRFY", mol)
genmethods.parseMTB("FRFY", mol)
with open("OutputFiles/TEST/FRFY_rotating.pdb", "w") as fh:
    for dihed in mol.getDihedrals():
        # get the atoms across the dihedral we need to rotate
        atomA = dihed.getAtmBIndex()
        atomB = dihed.getAtmCIndex()
        print(atomA,atomB)
        # calculate the vector between them
        mol.centreOn(atomB)
        vect = mol.getAtmWithIndex(atomB).getXYZ(vec=True) - mol.getAtmWithIndex(atomA).getXYZ(vec=True)
        vect = vect/numpy.linalg.norm(vect)
        
        numRotations = 36
        angle = (2*numpy.pi)/numRotations
        #angle = 360/numRotations
        print(angle)
        movingAtms = genmethods.graphSearch(atomA, mol, searchType="all", prescreened=[atomB], resetChecked=True)
        #os.mkdir("OutputFiles/TEST")
        
        print(movingAtms)
        
        
        for rotation in range(numRotations+1):
            rotAngle = rotation*angle
            tempMol = copy.deepcopy(mol)
            rotations.rotateMolecule(tempMol, tempMol, vect=vect, ang=rotation*angle, applyTo=movingAtms)
            recentre = genmethods.determineGeoCentre(tempMol, list(range(tempMol.getNumAtms())))
            tempMol.centreOn(recentre)
            fh.write(tempMol.printPDB())
            fh.write("ENDMDL\n")
            del tempMol