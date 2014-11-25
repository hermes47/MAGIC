#!/usr/bin/env python3
import genmethods, molecule, os, rotations, copy, sys, getopt, os
import numpy as np
from config import NUM_DIHED_ROTATE
''' Implementing the dihedral improvements process. Assume the input is an optimised PDB and AA MTB '''

# wrapper function for calling from command line
def main(argv):
    ifpfile = '54A7.ifp'
    source = 'O1QJ'
    runtype = 'generate'
    try:
        opts, args = getopt.getopt(argv, "i: f: t:",["ifpfile=","source=","runtype="])
        print("opts: ",opts)
        print("args: ",args)
    except getopt.GetoptError:
        print("magic.py -f <SourceFileRootName> -i <IFPFile> -t <RunType>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifpfile"):
            ifpfile = arg
        elif opt in ("-f", "--source"):
            source = arg
        elif opt in ("-t", "--runtype"):
            runtype = arg
    if runtype == 'generate': generateDihedralFiles(source, ifpfile)
    elif runtype == 'extract': pass
    elif runtype == 'picknfit': pass

def generateDihedralFiles(source, ifp):
    mol = molecule.Molecule()
    genmethods.parsePDB(source, mol)
    genmethods.parseMTB(source, mol)
    print(genmethods.atomConnectivity(4, mol, 'dihedrals'))
    for dihed in mol.getDihedrals():
        atms = dihed.getAtms()
        name = "-".join(str(x) for x in atms)
        print("Checking if in a ring") 
        modMol = copy.deepcopy(mol) # so we don't mess up the dihedrals to come
        print("Setting dihedral angle to 0")
        modMol.centreOn(atms[1])
        vecs = [] # the three bond vectors that define the dihedral
        for i in (0,1,2):
            vecs.append(modMol.getAtmWithIndex(atms[i]).getXYZ(vec=True)-modMol.getAtmWithIndex(atms[i+1]).getXYZ(vec=True))
            vecs[i] = vecs[i]/np.linalg.norm(vecs[i])
        toRot = genmethods.graphSearch(atms[2], modMol, resetChecked=True, searchType="all", prescreened=[atms[1]])
        # rotate the dihedral to 0
        n1 = rotations.normCrossProduct(vecs[0], vecs[1])
        n2 = rotations.normCrossProduct(vecs[1], vecs[2])
        dihedAng = rotations.vectAngle(n1,n2)
        modMol.rotate(rotations.rotationMatrix(vecs[1],dihedAng),atms=toRot)
        # set up for rotating a little bit at a time
        rotAng = 2*np.pi/NUM_DIHED_ROTATE
        rotMat = rotations.rotationMatrix(vecs[1],rotAng)
        os.system('mkdir -p DihedralFiles/%s/PDB/%s/' %(source,name))
        for step in range(NUM_DIHED_ROTATE):
            os.system('mkdir -p DihedralFiles/%s/{PDB,GAMESS-RS,GAMESS-BS,GAMESS-DS,GAMESS-RN,GAMESS-BN,GAMESS-DN}/%s/' %(source,name))
            with open('DihedralFiles/'+source+'/PDB/'+name+'/%02d.pdb' % step,'w') as fh:
                fh.write(modMol.printPDB())
            with open('DihedralFiles/'+source+'/GAMESS-RN/'+name+'/%02d.inp' % step, 'w') as fh:
                fh.write(modMol.printGAMESS(dihed, step*rotAng*180/np.pi, type='rigid', solvent=False))
            with open('DihedralFiles/'+source+'/GAMESS-RS/'+name+'/%02d.inp' % step, 'w') as fh:
                fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, type='rigid'))
            #with open('DihedralFiles/'+source+'/GAMESS-BN/'+name+'/%02d.inp' % step, 'w') as fh:
            #    fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, type='bonds', solvent=False))
            #with open('DihedralFiles/'+source+'/GAMESS-BS/'+name+'/%02d.inp' % step, 'w') as fh:
            #    fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, type='bonds'))
            with open('DihedralFiles/'+source+'/GAMESS-DN/'+name+'/%02d.inp' % step, 'w') as fh:
                fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, solvent=False))
            with open('DihedralFiles/'+source+'/GAMESS-DS/'+name+'/%02d.inp' % step, 'w') as fh:
                fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi))
            modMol.rotate(rotMat, atms=toRot)
            #print()
        #for i in range(NUM_DIHED_ROTATE):
            
    print("Generating")          

#print(os.getcwd())
#mol = molecule.Molecule()
#genmethods.parsePDB("O1QJ", mol)
#genmethods.parseMTB("O1QJ", mol)
#count=0
#for dihed in mol.getDihedrals():
#    # get the atoms across the dihedral we need to rotate
#    atomA = dihed.getAtmBIndex()
#    atomB = dihed.getAtmCIndex()
#    # calculate the vector between them
#    mol.centreOn(atomB)
#    vect = mol.getAtmWithIndex(atomB).getXYZ(vec=True) - mol.getAtmWithIndex(atomA).getXYZ(vec=True)
#    vect = vect/numpy.linalg.norm(vect)
#    
#    numRotations = 36
#    angle = (2*numpy.pi)/numRotations
#    #angle = 360/numRotations
#    movingAtms = genmethods.graphSearch(atomA, mol, searchType="all", prescreened=[atomB], resetChecked=True)
#    #os.mkdir("OutputFiles/TEST")
#    #with open("OutputFiles/TEST/O1QJ_%d_DH-%s.pdb" %(int(round(angle*180/numpy.pi,0)),str(dihed.getAtms())), "w") as fh:
#    for rotation in range(numRotations+1):
#        rotAngle = rotation*angle
#        tempMol = copy.deepcopy(mol)
#        rotations.rotateMolecule(tempMol, tempMol, vect=vect, ang=rotation*angle, applyTo=movingAtms)
#        recentre = genmethods.determineGeoCentre(tempMol, list(range(tempMol.getNumAtms())))
#        tempMol.centreOn(recentre)
#        atms = dihed.getAtms()
#        vect1 = tempMol.getAtmWithIndex(atms[1]).getXYZ(vec=True) - tempMol.getAtmWithIndex(atms[0]).getXYZ(vec=True)
#        vect2 = tempMol.getAtmWithIndex(atms[3]).getXYZ(vec=True) - tempMol.getAtmWithIndex(atms[2]).getXYZ(vec=True)
#        dihedralAngle = rotations.rotationAngleDetermine(vect1, vect2)
#        print(round(dihedralAngle*180/numpy.pi,3))
#        with open("OutputFiles/TEST/O1QJ_DH-%d_Rot-%d.pdb" %(count,int(round(rotAngle*180/numpy.pi,0))), "w") as fh:
#            fh.write(tempMol.printPDB())
#            fh.write("ENDMDL\n")
#        name = "rigid_O1QJ_%d_%d" %(count,int(round(rotAngle*180/numpy.pi,0)))
#        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
#            fh.write("%chk="+name+".chk\n")
#            fh.write("%mem=2gb\n")
#            fh.write("%NProcShared=1\n")
#            fh.write("#p b3lyp/6-31g(d)\n\n")
#            fh.write("SP "+name)
#            fh.write("\n\n0 1\n")
#            fh.write(tempMol.printChargeData())
#            fh.write("\n\n")
#        name = "fixedDihedral_O1QJ_%d_%d" %(count,int(round(rotAngle*180/numpy.pi,0)))
#        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
#            fh.write("%chk="+name+".chk\n")
#            fh.write("%mem=2gb\n")
#            fh.write("%NProcShared=1\n")
#            fh.write("#p b3lyp/6-31g(d) opt=modredundant\n\n")
#            fh.write("SP "+name)
#            fh.write("\n\n0 1\n")
#            fh.write(tempMol.printChargeData())
#            fh.write("\n"+str(dihed.getAtmAPos())+" "+str(dihed.getAtmBPos())+" "+str(dihed.getAtmCPos())+" "+str(dihed.getAtmDPos())+" F\n\n")
#        name = "fixedBonds_O1QJ_%d_%d" %(count,int(round(rotAngle*180/numpy.pi,0)))
#        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
#            fh.write("%chk="+name+".chk\n")
#            fh.write("%mem=2gb\n")
#            fh.write("%NProcShared=1\n")
#            fh.write("#p b3lyp/6-31g(d) opt=modredundant\n\n")
#            fh.write("SP "+name)
#            fh.write("\n\n0 1\n")
#            fh.write(tempMol.printChargeData())
#            fh.write("\n"+str(dihed.getAtmAPos())+" "+str(dihed.getAtmBPos())+" "+str(dihed.getAtmCPos())+" "+str(dihed.getAtmDPos())+" F\n")
#            for bond in tempMol.getBonds():
#                fh.write(str(bond.getAtmAPos())+" "+str(bond.getAtmBPos())+" F\n")
#            fh.write("\n\n")
#        del tempMol
#    count +=1
#    
#
if __name__ == "__main__":
    main(sys.argv[1:])