import genmethods, molecule, os, numpy, rotations, copy

print(os.getcwd())
mol = molecule.Molecule()
genmethods.parsePDB("O1QJ", mol)
genmethods.parseMTB("O1QJ", mol)
count=0
for dihed in mol.getDihedrals():
    # get the atoms across the dihedral we need to rotate
    atomA = dihed.getAtmBIndex()
    atomB = dihed.getAtmCIndex()
    # calculate the vector between them
    mol.centreOn(atomB)
    vect = mol.getAtmWithIndex(atomB).getXYZ(vec=True) - mol.getAtmWithIndex(atomA).getXYZ(vec=True)
    vect = vect/numpy.linalg.norm(vect)
    
    numRotations = 36
    angle = (2*numpy.pi)/numRotations
    #angle = 360/numRotations
    movingAtms = genmethods.graphSearch(atomA, mol, searchType="all", prescreened=[atomB], resetChecked=True)
    #os.mkdir("OutputFiles/TEST")
    #with open("OutputFiles/TEST/O1QJ_%d_DH-%s.pdb" %(int(round(angle*180/numpy.pi,0)),str(dihed.getAtms())), "w") as fh:
    for rotation in range(numRotations+1):
        rotAngle = rotation*angle
        tempMol = copy.deepcopy(mol)
        rotations.rotateMolecule(tempMol, tempMol, vect=vect, ang=rotation*angle, applyTo=movingAtms)
        recentre = genmethods.determineGeoCentre(tempMol, list(range(tempMol.getNumAtms())))
        tempMol.centreOn(recentre)
        atms = dihed.getAtms()
        vect1 = tempMol.getAtmWithIndex(atms[1]).getXYZ(vec=True) - tempMol.getAtmWithIndex(atms[0]).getXYZ(vec=True)
        vect2 = tempMol.getAtmWithIndex(atms[3]).getXYZ(vec=True) - tempMol.getAtmWithIndex(atms[2]).getXYZ(vec=True)
        dihedralAngle = rotations.rotationAngleDetermine(vect1, vect2)
        print(round(dihedralAngle*180/numpy.pi,3))
        with open("OutputFiles/TEST/O1QJ_DH-%d_Rot-%d.pdb" %(count,int(round(rotAngle*180/numpy.pi,0))), "w") as fh:
            fh.write(tempMol.printPDB())
            fh.write("ENDMDL\n")
        name = "rigid_O1QJ_%d_%d" %(count,int(round(rotAngle*180/numpy.pi,0)))
        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
            fh.write("%chk="+name+".chk\n")
            fh.write("%mem=2gb\n")
            fh.write("%NProcShared=1\n")
            fh.write("#p b3lyp/6-31g(d)\n\n")
            fh.write("SP "+name)
            fh.write("\n\n0 1\n")
            fh.write(tempMol.printChargeData())
            fh.write("\n\n")
        name = "fixedDihedral_O1QJ_%d_%d" %(count,int(round(rotAngle*180/numpy.pi,0)))
        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
            fh.write("%chk="+name+".chk\n")
            fh.write("%mem=2gb\n")
            fh.write("%NProcShared=1\n")
            fh.write("#p b3lyp/6-31g(d) opt=modredundant\n\n")
            fh.write("SP "+name)
            fh.write("\n\n0 1\n")
            fh.write(tempMol.printChargeData())
            fh.write("\n"+str(dihed.getAtmAPos())+" "+str(dihed.getAtmBPos())+" "+str(dihed.getAtmCPos())+" "+str(dihed.getAtmDPos())+" F\n\n")
        name = "fixedBonds_O1QJ_%d_%d" %(count,int(round(rotAngle*180/numpy.pi,0)))
        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
            fh.write("%chk="+name+".chk\n")
            fh.write("%mem=2gb\n")
            fh.write("%NProcShared=1\n")
            fh.write("#p b3lyp/6-31g(d) opt=modredundant\n\n")
            fh.write("SP "+name)
            fh.write("\n\n0 1\n")
            fh.write(tempMol.printChargeData())
            fh.write("\n"+str(dihed.getAtmAPos())+" "+str(dihed.getAtmBPos())+" "+str(dihed.getAtmCPos())+" "+str(dihed.getAtmDPos())+" F\n")
            for bond in tempMol.getBonds():
                fh.write(str(bond.getAtmAPos())+" "+str(bond.getAtmBPos())+" F\n")
            fh.write("\n\n")
        del tempMol
    count +=1