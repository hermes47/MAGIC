'''
Created on 20/09/2014

@author: iwelsh
'''
import genmethods
import molecule
import os

os.chdir("/Users/iwelsh/tmp/")

fileName = "/Users/iwelsh/tmp/pdbfiles.dat"

with open(fileName,"r") as fh:
    PDBFiles = fh.readlines()

fileCount=0
    
for pdb in PDBFiles:
    mol = molecule.Molecule()
    passedData = genmethods.loadFile(pdb[:-5], "/Users/iwelsh/tmp/", ".pdb", True)
    genmethods.parsePDB("rootName", mol, False, passedData)
    newName="MTMH_File_"+str(fileCount)
    with open(newName+".gau","w") as f:
        f.write("%chk="+newName+".chk\n")
        f.write("%mem=8gb\n")
        f.write("%NProcShared=4\n")
        f.write("#p b3lyp/6-311+g(d,p) scrf=(smd,solvent=water)\n\n")
        f.write("SP for deltaG in solvent for "+newName)
        f.write("\n\n1 1\n")
        f.write(mol.printChargeData())
        f.write("\n\n--Link1--\n%chk="+newName+".chk\n%mem=8gb\n%NProcShared=4\n")
        f.write("#p b3lyp/6-311+g(d,p) geom=checkpoint\n\n")
        f.write("SP for deltaG in vacuum for "+newName)
        f.write("\n\n1 1\n\n\n")
    fileCount += 1
    print("file written")
    

