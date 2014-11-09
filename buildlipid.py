# import the needed files
import sys, os
import connect
import genmethods
from ifp import IFP
from molecule import Molecule
from config import *
import pickle

# main loop for joining mtb files. Pass it an IFP file, followed by the string detailing how you want things to connect together
# and finally all the MTB/PDB root names that are passed
def MAGIC(ifpFile, connectString, *args):
    # load the IFP file into memory and pass it through the IFP parser
    #with open(ifpFile, "r") as fh:
    #    file = fh.readlines()
    #ifpData = IFP(file)
    print("Working from: %s" % ROOT_DIRECTORY)
    ## LOAD ALL MOLECULES
    '''
    Current: molecules present in the sanatised connect string are loaded from either
    a pickled file or the PDB, MTB and OVL files.
    '''
    # sanatise the input string. because formatting is an annoyance
    sanity = genmethods.stringSanitation(connectString)
    # split into a list, based on NAME_LENGTH
    for c in range(0,len(sanity),NAME_LENGTH):
        rootName = sanity[c:c+NAME_LENGTH]
        if rootName+".dat" in SAVED_PICKLED_FILES and rootName not in LOADED_MOLECULES:
            with open("%s/PickledFiles/%s.dat" %(ROOT_DIRECTORY, rootName), "rb") as fh:
                print("Loading "+rootName+" from pickled data")
                LOADED_MOLECULES[rootName] = pickle.load(fh)
        elif rootName not in LOADED_MOLECULES: # only add unique molecules
            # make the new molecule
            mol = Molecule()
            # go through the PDB, MTB, and OVL data, in order
            print("Loading "+rootName+" from files")
            genmethods.parsePDB(rootName, mol)
            genmethods.parseMTB(rootName, mol)
            genmethods.parseOVL(rootName, mol)
            # add the complete molecule to the LOADED_MOLECULES list
            LOADED_MOLECULES[rootName] = mol
            # as well as saving it as pickle data for later
            #with open("%s/PickledFiles/%s.dat" %(ROOT_DIRECTORY, rootName), "wb") as fh:
            #    pickle.dump(mol, fh, pickle.HIGHEST_PROTOCOL)
    ## END LOAD ALL MOLECULES
    # make sure we aren't going to hit the recursion limit
    recurseCount = 0
    for c in STRING_BREAK_CHARACTERS:
        recurseCount += connectString.count(c)
    if not recurseCount >= sys.getrecursionlimit():
        # parse the connectString
        output = connect.initStringParse(connectString)
        print("Final structure output as %s" % output)
        os.chdir(ROOT_DIRECTORY)
        os.system('mkdir OutputFiles/'+output)
        # print out final file to PDB and MTB files
        LOADED_MOLECULES[output].centreOn(list(genmethods.determineGeoCentre(LOADED_MOLECULES[output], range(LOADED_MOLECULES[output].getNumAtms()))))
        with open("OutputFiles/"+output+"/"+output+".pdb", "w") as fh:
            fh.write(LOADED_MOLECULES[output].printPDB())
        with open("OutputFiles/"+output+"/"+output+".mtb", "w") as fh:
            fh.write(LOADED_MOLECULES[output].printMTB())
        with open("OutputFiles/"+output+"/"+output+".dat", "w") as fh:
            fh.write(LOADED_MOLECULES[output].printChargeData())
        # run an energy minimisation on the PDB file, using GROMOS
        print("Running an energy minimisation on the joint PDB structure")
        genmethods.runGROMOS(output, connectString)
        os.chdir(ROOT_DIRECTORY)
        #with open("PickledFiles/"+output+".dat", "wb") as fh:
        #    pickle.dump(LOADED_MOLECULES[output], fh, pickle.HIGHEST_PROTOCOL)
        #with open("PickledFiles/"+output+".dat", "rb") as fh:
        #    tempMol = pickle.load(fh)
        #print(tempMol.printPDB())
    else:
        print("Recursion limit could be breached. Please submit string sequence in  smaller parts.")
    




''' run the MTBConnect script '''
if __name__ == "__main__":
    argv = ""
    ifpFile = "54a7.ifp"
    inputString = "PbEO(BigH{2}{MCWs})"
    MAGIC(ifpFile,inputString)