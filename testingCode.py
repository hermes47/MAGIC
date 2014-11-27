#!/usr/bin/env python3

import genmethods, molecule, copy, numpy, rotations, statistics

print(statistics.mean([1,2,3,4,5,6,7,8,9,10,0]),statistics.stdev([0,1,2,3,4,5,6,7,8,9,10]))

mol = molecule.Molecule()
genmethods.parsePDB("O1QJ", mol)
genmethods.parseMTB("O1QJ", mol)
tempMol = copy.deepcopy(mol)
scaleFactor = 2.5

for atm in tempMol.getAtms():
    newCoords = atm.getXYZ(vec=True)*scaleFactor
    atm.setXYZ(newCoords, vec=True)

graph = genmethods.genGraphRep(mol)

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
    
find_all_paths_of_length(graph, 4)
#print("Comparing old parameters with new.\nBonds:")
#for bond in mol.getBonds():
#    atmA = bond.getAtmAIndex()
#    atmB = bond.getAtmBIndex()
#    distOld = numpy.linalg.norm(mol.getAtmWithIndex(atmA).getXYZ(vec=True) - mol.getAtmWithIndex(atmB).getXYZ(vec=True))
#    distNew = numpy.linalg.norm(tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) - tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
#    oldToNew = distOld/distNew
#    print("%2.4f   %2.4f   %1.1f"%(distOld,distNew,oldToNew))
#print("\nAngles:")
#for angle in mol.getAngles():
#    atmA = angle.getAtmAIndex()
#    atmB = angle.getAtmBIndex()
#    atmC = angle.getAtmCIndex()
#    vect1Old = mol.getAtmWithIndex(atmA).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
#    vect2Old = mol.getAtmWithIndex(atmC).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
#    angOld = rotations.rotationAngleDetermine(vect1Old, vect2Old)*180.0/numpy.pi
#    vect1New = (tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
#    vect2New = (tempMol.getAtmWithIndex(atmC).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
#    angNew = rotations.rotationAngleDetermine(vect1New, vect2New)*180.0/numpy.pi
#    print("%3.4f    %3.4f    %1.1f"%(angOld,angNew,angOld/angNew))
#print("\nDihedrals:")
#for dihedral in mol.getDihedrals():
#    atmD = dihedral.getAtmDIndex()
#    atmA = dihedral.getAtmAIndex()
#    atmB = dihedral.getAtmBIndex()
#    atmC = dihedral.getAtmCIndex()
#    vect1Old = mol.getAtmWithIndex(atmA).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
#    vect2Old = mol.getAtmWithIndex(atmD).getXYZ(vec=True) -  mol.getAtmWithIndex(atmC).getXYZ(vec=True)
#    angOld = rotations.rotationAngleDetermine(vect1Old, vect2Old)*180.0/numpy.pi
#    vect1New = (tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
#    vect2New = (tempMol.getAtmWithIndex(atmD).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmC).getXYZ(vec=True))
#    angNew = rotations.rotationAngleDetermine(vect1New, vect2New)*180.0/numpy.pi
#    print("%3.4f    %3.4f    %1.1f"%(angOld,angNew,angOld/angNew))
#print("\nImpropers:")
#for dihedral in mol.getImpropers():
#    atmD = dihedral.getAtmDIndex()
#    atmA = dihedral.getAtmAIndex()
#    atmB = dihedral.getAtmBIndex()
#    atmC = dihedral.getAtmCIndex()
#    vect1Old = mol.getAtmWithIndex(atmA).getXYZ(vec=True) -  mol.getAtmWithIndex(atmB).getXYZ(vec=True)
#    vect2Old = mol.getAtmWithIndex(atmD).getXYZ(vec=True) -  mol.getAtmWithIndex(atmC).getXYZ(vec=True)
#    angOld = rotations.rotationAngleDetermine(vect1Old, vect2Old)*180.0/numpy.pi
#    vect1New = (tempMol.getAtmWithIndex(atmA).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmB).getXYZ(vec=True))
#    vect2New = (tempMol.getAtmWithIndex(atmD).getXYZ(vec=True) -  tempMol.getAtmWithIndex(atmC).getXYZ(vec=True))
#    angNew = rotations.rotationAngleDetermine(vect1New, vect2New)*180.0/numpy.pi
#    print("%3.4f    %3.4f    %1.1f"%(angOld,angNew,angOld/angNew))

#import os
#count = 0
#for length in range(5,37):
#    maxDBs = int(((length-5)/3) + 1)
#    for numDBs in range(1,maxDBs+1):
#        highest = length - 3*numDBs + 1
#        for startPos in range(3,highest+1):
#            precedingAlk = startPos - 3
#            proceedingAlk = length - 2 - precedingAlk - numDBs*3
#            if startPos == 3:
#                os.system('./magic.py -s "CIST{%d}{CISI}{%d}{ALKI}ENDG"' %(numDBs-1,length-4-3*(numDBs-1)))
#                os.system('./magic.py -s "TRAT{%d}{TRAI}{%d}{ALKI}ENDG"' %(numDBs-1,length-4-3*(numDBs-1)))
#                count += 2
#            else:
#                os.system('./magic.py -s "ALKT{%d}{ALKI}{%d}{CISI}{%d}{ALKI}ENDG"' %(precedingAlk,numDBs,proceedingAlk))
#                os.system('./magic.py -s "ALKT{%d}{ALKI}{%d}{TRAI}{%d}{ALKI}ENDG"' %(precedingAlk,numDBs,proceedingAlk))
#                count += 2
#print(count)
#import numpy as np
#import scipy.optimize as sp
#
#angles = np.array([0,0.174532925,0.34906585,0.523598776,0.698131701,0.872664626,1.047197551,1.221730476,1.396263402,1.570796327,1.745329252,
#                   1.919862177,2.094395102,2.268928028,2.443460953,2.617993878,2.792526803,2.967059728,3.141592654,3.316125579,3.490658504,
#                   3.665191429,3.839724354,4.01425728,4.188790205,4.36332313,4.537856055,4.71238898,4.886921906,5.061454831,5.235987756,
#                   5.410520681,5.585053606,5.759586532,5.934119457,6.108652382,6.283185307])
#
#energies = np.array([0,0.713721171,2.938630257,6.017711638,9.163094769,11.54534482,12.47505537,11.77926374,9.666812697,6.660124229,3.519264833,
#                     1.11286806,0.110026828,0.81282592,3.029351784,6.107495861,9.231798853,11.6052456,12.53574906,11.80469696,9.697972131,6.671613417,
#                     3.519154563,1.115317651,0.1081706,0.815653583,3.05046343,6.13676231,9.270674632,11.62124015,12.52827163,11.78622919,
#                     9.644572087,6.576213249,3.408831053,1.002510418,0])
#
#def func(x, k, m, phase, shift):
#    return k*(1+np.cos(m*x+phase)) + shift
#results  = sp.curve_fit(func, angles, energies, p0=(5.0,3,1,12))
#print(results)
#Fk = np.fft.fft(energies)
#print(Fk)