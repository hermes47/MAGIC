#!/usr/bin/env python3
import genmethods, molecule, rotations, copy, sys, getopt, os, openbabel
import numpy as np
from config import NUM_DIHED_ROTATE, K_RMSD_TOL, PHI_ABS_TOL, K_PER_VAR_TOL, FFT_TARGET_DP, PHI_ROUND_VAL
import matplotlib.pyplot as plt
''' Implementing the dihedral improvements process. Assume the input is an optimised PDB and AA MTB '''

# wrapper function for calling from command line
def main(argv):
    ifpfile = '54A7'
    source = 'O1QJ'
    runtype = 'full'
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
    print(os.getcwd())
    mol, diheds, ifpfile = loadData(source, ifpfile)
    if runtype == 'generate': generateDihedralFiles(source, ifpfile, mol, diheds)
    elif runtype == 'extract': extractEnergies(source, ifpfile, mol, diheds)
    elif runtype == 'fit': fitDihedralTerms(source, mol, diheds)
    elif runtype == 'set': setDihedralTerms(source, ifpfile, mol, diheds)
    elif runtype == 'full':
        extractEnergies(source, ifpfile, mol, diheds)
        os.chdir('/Users/iwelsh/GitHub/MAGIC/')
        fitDihedralTerms(source, mol, diheds)
        os.chdir('/Users/iwelsh/GitHub/MAGIC/')
        setDihedralTerms(source, ifpfile, mol, diheds)
        os.chdir('/Users/iwelsh/GitHub/MAGIC/')

def loadData(source, ifp):
    mol = molecule.Molecule()
    genmethods.parsePDB(source, mol)
    genmethods.parseMTB(source, mol)
    uniqueDiheds = []
    for dihed in mol.getDihedrals():
        #name = "-".join(str(x) for x in dihed.getAtms())
        name = dihed.getAtms()
        if name not in uniqueDiheds: uniqueDiheds.append(name)
        # check if the dihed is in a ring
    ifpFile = genmethods.parseIFP(ifp, os.getcwd()+'/DihedralFiles/')
    return mol, uniqueDiheds, ifpFile

def setDihedralTerms(source, ifp, mol, diheds):
    mol.clearDihedrals()
    newDihedralTypesAdded = False
    for dihed in diheds:
        name = "-".join(str(x) for x in dihed)
        atmIDs = name.split('-')
        with open('DihedralFiles/%s/ANAFiles/%s_AllPeaks.txt'%(source,name),'r') as fh:
            lines = fh.readlines()
        peaks = []
        for lin in lines:
            if lin.startswith('#'): pass
            else:
                lin = lin.split()
                peaks.append({'k':float(lin[0]),'sign':int(lin[1]),'m':int(lin[2]),'phi0':float(lin[3]),'rmsd':float(lin[4])})
        currentrmsd = 100
        for peak in peaks:
            if currentrmsd < K_RMSD_TOL: continue
            multioptions = ifp.getTorsWithMultiplicity(peak['m']) # multiplicity must match exactly
            highlypossible = []
            for opt in multioptions:
                if not peak['k']*(1-K_PER_VAR_TOL) < opt['CP'] < peak['k']*(1+K_PER_VAR_TOL): continue
                opt['kdiff'] = np.absolute(opt['CP']-peak['k'])
                if not peak['phi0'] - PHI_ABS_TOL < opt['PD'] < peak['phi0'] + PHI_ABS_TOL: continue 
                opt['phidiff'] = np.absolute(opt['PD']-peak['phi0'])
                highlypossible.append(opt)
            # actually make the choice now that it's been narrowed down
            if len(highlypossible) == 0: # needs a new parameter
                currentMax = ifp.MaxTorsType()
                newDihedral = {'N':currentMax+1,'PD':peak['phi0'],'CP':peak['k'],'MP':peak['m']}
                ifp.addTorsDihedralType(newDihedral)
                dihedralData = (int(atmIDs[0])+1,int(atmIDs[1])+1,int(atmIDs[2])+1,int(atmIDs[3])+1,currentMax+1)
                mol.addDihedral(dihedralData)
                newDihedralTypesAdded = True
                print('Made new dihedral parameter type %d'%(currentMax+1))
            elif len(highlypossible) == 1: # only one possible choice
                dihedralData = (int(atmIDs[0])+1,int(atmIDs[1])+1,int(atmIDs[2])+1,int(atmIDs[3])+1,highlypossible[0]['N'])
                mol.addDihedral(dihedralData)
                print('Only choice was %d'%highlypossible[0]['N'])
            elif len(highlypossible) > 1: # need to pick the best option
                currentMax = 100
                bestChoice = {}
                for opt in highlypossible:
                    if opt['kdiff'] < currentMax:
                        bestChoice = opt
                        currentMax = opt['kdiff']
                dihedralData = (int(atmIDs[0])+1,int(atmIDs[1])+1,int(atmIDs[2])+1,int(atmIDs[3])+1,bestChoice['N'])
                mol.addDihedral(dihedralData)
                print('Picked dihedral %s from choice of: '%bestChoice['N'] + ', '.join(str(x['N']) for x in highlypossible))
            currentrmsd = peak['rmsd']
    with open('DihedralFiles/%s/%s_mod.mtb'%(source,source),'w') as fh:
        fh.write(mol.printMTB())
    if newDihedralTypesAdded: 
        with open('DihedralFiles/%s.ifp'%ifp.getForcefield(),'w') as fh:
            fh.write(ifp.printIFP())

def generateDihedralFiles(source, ifp, mol, diheds):
    #graphRep = genmethods.genGraphRep(mol)
    for dihed in diheds:
        name = "-".join(str(x) for x in dihed)
        modMol = copy.deepcopy(mol) # so we don't mess up the dihedrals to come
        print("Setting dihedral angle to 0")
        modMol.centreOn(dihed[1])
        vecs = [] # the three bond vectors that define the dihedral
        for i in (0,1,2):
            vecs.append(modMol.getAtmWithIndex(dihed[i]).getXYZ(vec=True)-modMol.getAtmWithIndex(dihed[i+1]).getXYZ(vec=True))
            vecs[i] = vecs[i]/np.linalg.norm(vecs[i])
        toRot = genmethods.graphSearch(dihed[2], modMol, resetChecked=True, searchType="all", prescreened=[dihed[1]])
        # rotate the dihedral to 0
        n1 = rotations.normCrossProduct(vecs[0], vecs[1])
        n2 = rotations.normCrossProduct(vecs[1], vecs[2])
        m1 = rotations.normCrossProduct(n1,vecs[1])
        x = np.dot(n1,n2)
        y = np.dot(m1,n2)
        dihedAng = np.arctan2(y,x)
        modMol.rotate(rotations.rotationMatrix(vecs[1],dihedAng),atms=toRot)
        # set up for rotating a little bit at a time
        rotAng = 2*np.pi/NUM_DIHED_ROTATE
        rotMat = rotations.rotationMatrix(vecs[1],rotAng)
        for step in range(NUM_DIHED_ROTATE):
            os.system('mkdir -p DihedralFiles/%s/{PDB,GAMESS-RS,GAMESS-BS,GAMESS-DS,GAMESS-RN,GAMESS-BN,GAMESS-DN,Gaussian}/%s/' %(source,name))
            with open('DihedralFiles/'+source+'/PDB/'+name+'/%02d.pdb' % step,'w') as fh:
                fh.write(modMol.printPDB())
            with open('DihedralFiles/'+source+'/GAMESS-RN/'+name+'/%02d.inp' % step, 'w') as fh:
                fh.write(modMol.printGAMESS(dihed, step*rotAng*180/np.pi, kind='rigid', solvent=False))
            #with open('DihedralFiles/'+source+'/GAMESS-RS/'+name+'/%02d.inp' % step, 'w') as fh:
            #    fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, type='rigid'))
            #with open('DihedralFiles/'+source+'/GAMESS-BN/'+name+'/%02d.inp' % step, 'w') as fh:
            #    fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, type='bonds', solvent=False))
            #with open('DihedralFiles/'+source+'/GAMESS-BS/'+name+'/%02d.inp' % step, 'w') as fh:
            #    fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, type='bonds'))
            #with open('DihedralFiles/'+source+'/GAMESS-DN/'+name+'/%02d.inp' % step, 'w') as fh:
            #    fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi, solvent=False))
            #with open('DihedralFiles/'+source+'/GAMESS-DS/'+name+'/%02d.inp' % step, 'w') as fh:
            #    fh.write(modMol.printGAMESS(dihed,step*rotAng*180/np.pi))
            with open('DihedralFiles/'+source+'/Gaussian/'+name+'/%02d.gau' % step, 'w') as fh:
                fh.write("%"+"chk=%02d.chk\n"%step)
                fh.write("%mem=4gb\n")
                fh.write("%NProcShared=2\n")
                fh.write("#p b3lyp/6-31g(d) opt=modredundant\n\n")
                fh.write("SP "+name)
                fh.write("\n\n0 1\n")
                fh.write(modMol.printChargeData())
                fh.write("\n"+str(dihed[0]+1)+" "+str(dihed[1]+1)+" "+str(dihed[2]+1)+" "+str(dihed[3]+1)+" F\n\n")
            modMol.rotate(rotMat, atms=toRot)
            #print()
        #for i in range(NUM_DIHED_ROTATE):
            
def extractEnergies(source, ifp, mol, diheds, kind='GAMESS-RN'):
    folders = os.listdir(path='/Users/iwelsh/GitHub/MAGIC/DihedralFiles/'+source+'/'+kind+'/')
    try:
        os.mkdir('/Users/iwelsh/GitHub/MAGIC/DihedralFiles/'+source+'/ANAFiles/')
    except: pass
    for dihed in diheds:
        qmenergies = {}
        mdenergies = {}
        name = "-".join(str(x) for x in dihed)
        if name not in folders: continue
        for i in range(NUM_DIHED_ROTATE):
            os.chdir('/Users/iwelsh/GitHub/MAGIC/DihedralFiles/'+source+'/')
            filename = '%02d' % i
            with open(kind+'/'+name+'/'+filename+'.log', 'r') as fh:
                logdata = fh.readlines()
            extractLines = False
            finalData = ''
            if kind == 'Gaussian':
                for line in logdata:
                    if extractLines: finalData += line[1:-1]
                    if 'l9999.exe' in line: extractLines = True
                    if '@' in line: extractLines = False
                finalData = finalData.split('\\')
            elif kind == 'GAMESS-RN':
                for line in logdata:
                    if 'TOTAL ENERGY' in line: finalData = line.split()
            if kind == 'Gaussian':
                obConversion = openbabel.OBConversion()
                obConversion.SetInAndOutFormats("g03","pdb")
                obMol = openbabel.OBMol()
                #obConversion.ReadFile(obMol, kind+'/'+name+'/'+filename+'.log')
                obConversion.ReadString(obMol, logdata)
                pdbData = obConversion.WriteString(obMol).splitlines()[2:]
                for line in pdbData:
                    if line.startswith('HETATM') or line.startswith('ATOM'):
                        atmData = line.split()
                        mol.getAtmWithIndex(int(atmData[1])-1).setXYZ((float(atmData[5]),float(atmData[6]),float(atmData[7])))
                with open('PDB/'+name+'/'+filename+'.pdb', 'w') as fh:
                    fh.write(mol.printPDB())
            os.system('cp ../../GROMOSData/singlepoint.imd ./PDB/'+name+'/'+filename+'.imd')
            os.system("sed -i 's/NUMOFATOMS/%d/' ./PDB/%s/%s.imd" % (mol.getNumAtms(),name,filename))
            os.system('make_top @build ../../Input/MTBFiles/'+source+'.mtb ../../GROMOSData/54A7.mtb @param ../../GROMOSData/54A7.ifp @seq '+source+' @solv H2O > ./PDB/'+name+'/'+filename+'.top')
            os.system('pdb2g96 @topo ./PDB/'+name+'/'+filename+'.top @pdb ./PDB/'+name+'/'+filename+'.pdb @lib ../../GROMOSData/pdb2g96.lib @out ./PDB/'+name+'/'+filename+'.cnf 2>/dev/null')
            os.chdir('./PDB/'+name)
            os.system('md @topo %s.top @conf %s.cnf @fin %s_min.cnf @input %s.imd > %s.omd' %(filename, filename, filename,filename,filename))
            with open('./%s.omd' % filename, 'r') as fh:
                mdData = fh.readlines()
            for line in mdData:
                if line.startswith('E_Non-bonded'): mdenergies[filename] = float(line.split()[2])
            for item in finalData:
                if item.startswith('HF=') and kind == 'Gaussian': qmenergies[filename] = float(item[3:])
                if item.startswith('-') and kind == 'GAMESS-RN': qmenergies[filename] = float(item)
        mdReference = mdenergies['00']
        qmReference = qmenergies['00']
        coenergies = {}
        for energy in mdenergies: mdenergies[energy] = round(mdenergies[energy] - mdReference,10)
        for energy in qmenergies: qmenergies[energy] = round((qmenergies[energy] - qmReference)*2625.5,10)
        for energy in qmenergies: coenergies[energy] = qmenergies[energy]-mdenergies[energy]
        with open('../../ANAFiles/%s_QMEnergies.txt'%name,'w') as fh:
            for config in sorted(qmenergies): fh.write('%s     %f\n' %(config,qmenergies[config]))
        with open('../../ANAFiles/%s_MDEnergies.txt'%name,'w') as fh:
            for config in sorted(mdenergies): fh.write('%s     %f\n' %(config,mdenergies[config]))
        with open('../../ANAFiles/%s_CombinedEnergies.txt'%name,'w') as fh:
            for config in sorted(mdenergies): fh.write('%s     %f\n' %(config,qmenergies[config]-mdenergies[config]))
        x= []
        for i in range(NUM_DIHED_ROTATE): x.append(i*360/NUM_DIHED_ROTATE)
        tempvals = []
        for i in sorted(mdenergies): tempvals.append(mdenergies[i])
        plt.plot(x, tempvals, label='MD')
        del tempvals
        tempvals = []
        for i in sorted(qmenergies): tempvals.append(qmenergies[i])
        plt.plot(x, tempvals, label='QM')
        del tempvals
        tempvals = []
        for i in sorted(coenergies): tempvals.append(coenergies[i])
        plt.plot(x, tempvals, label='CB')
        plt.xlabel('phi')
        plt.ylabel('Energy kJ/mol')
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)
        plt.savefig('../../ANAFiles/%s_EnergiesPlot.eps'%name, format='eps')
        plt.close()
    print('extracting energies from calculations')
    
def fitDihedralTerms(source, mol, diheds, kind='GAMESS-RN'):
    folders = os.listdir(path='/Users/iwelsh/GitHub/MAGIC/DihedralFiles/'+source+'/'+kind+'/')
    os.chdir('/Users/iwelsh/GitHub/MAGIC/DihedralFiles/'+source+'/ANAFiles')
    for dihed in diheds:
        name = "-".join(str(x) for x in dihed)
        if name not in folders: continue
        with open('./%s_CombinedEnergies.txt'%name,'r') as fh:
            dataInput = fh.readlines()
        energies = []
        for lin in dataInput: energies.append(float(lin.split()[1]))
        for _ in range(int(round(FFT_TARGET_DP/NUM_DIHED_ROTATE)-1)):
            for j in range(NUM_DIHED_ROTATE): energies.append(energies[j])
        energyMean = np.mean(energies)
        energies -= energyMean # gets rid of DC shift
        n = len(energies)
        dx = 2*np.pi/NUM_DIHED_ROTATE
        x = dx*np.arange(0,n)
        Fk = np.fft.rfft(energies)/n # Fourier coefficients (divided by n) 
        nu = np.fft.rfftfreq(n,dx) # Natural frequencies
        Fk = np.fft.fftshift(Fk) # Shift zero freq to center
        nu = np.fft.fftshift(nu) # Shift zero freq to center
        energies += energyMean
        n=len(Fk)
        peaks = []
        for i in range(1,n-1):
            if np.absolute(Fk[i]) > np.absolute(Fk[i-1]) and np.absolute(Fk[i]) > np.absolute(Fk[i+1]): 
                peaks.append([np.absolute(np.real(Fk[i]))*2,np.sign(np.real(Fk[i])),round(nu[i]*2*np.pi,0),np.angle(Fk)[i]])
        for peak in peaks:
            peak[0] = peak[0]*peak[1]/np.cos(peak[3]) # adjusts k for phase
            #if peak[1] < 0: print(peak[3], np.cos(np.pi), np.cos(180))
            #peak[3] = peak[1]*peak[3]
        peaks = sorted(peaks,reverse=True)
        rmsd = 100.
        keptPeaks = []
        for i in range(len(peaks)):
            keptPeaks.append(peaks[i])
            calcValues = []
            for j in range(NUM_DIHED_ROTATE):
                jvalue = 0.
                for para in keptPeaks:
                    jvalue += para[0]*(1+np.cos(para[2]*x[j]+para[1]*para[3]))-para[0]*(1+np.cos(para[2]*x[0]+para[1]*para[3]))
                    #jvalue += para[0]*(1+np.cos(para[2]*x[j]))-para[0]*(1+np.cos(para[2]*x[0]))
                calcValues.append(jvalue)
            rmsd = np.sqrt(np.mean((calcValues-energies[0:NUM_DIHED_ROTATE])**2))
            keptPeaks[i].append(round(rmsd,4))
            ax = plt.subplot(111)
            ax.plot(x[0:NUM_DIHED_ROTATE]*180/np.pi, energies[0:NUM_DIHED_ROTATE], label='data')
            ax.plot(x[0:NUM_DIHED_ROTATE]*180/np.pi,calcValues[0:NUM_DIHED_ROTATE], label='%d parameter fit'%(i+1))
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2, fancybox=True, shadow=True)
            plt.xlabel('phi')
            plt.ylabel('Energy kJ/mol')
            plt.savefig('./%s_%03d_ParamPlot.eps'%(name,i+1), format='eps')
            plt.close()
        with open('./%s_AllPeaks.txt'%name,'w') as fh:
            fh.write('# k    sign    m    phi0    rmsd\n')
            for peak in keptPeaks:
                fh.write('%5.2f  %2d  %2d  %5d.0  %6.4f\n' %(peak[0],int(peak[1]),int(peak[2]),int(round(peak[3]*180/np.pi/PHI_ROUND_VAL)*PHI_ROUND_VAL),peak[4]))
                #fh.write('%f  %d  %d  %f  %f\n' %(peak[0],int(peak[1]),int(peak[2]),peak[3]*180/np.pi,peak[4]))
        plt.plot(nu*2*np.pi, np.real(Fk)*2)
        plt.xlabel('Multiplicity')
        plt.ylabel('Amplitude x2')
        plt.savefig('./%s_FFTPlot.eps'%name, format='eps')
        plt.close()
    print('fft fitting dihedral profiles')
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
#    vect = vect/np.linalg.norm(vect)
#    
#    numRotations = 36
#    angle = (2*np.pi)/numRotations
#    #angle = 360/numRotations
#    movingAtms = genmethods.graphSearch(atomA, mol, searchType="all", prescreened=[atomB], resetChecked=True)
#    #os.mkdir("OutputFiles/TEST")
#    #with open("OutputFiles/TEST/O1QJ_%d_DH-%s.pdb" %(int(round(angle*180/numpy.pi,0)),str(dihed.getAtms())), "w") as fh:
#    for rotation in range(numRotations+1):
#        rotAngle = rotation*angle
#        tempMol = copy.deepcopy(mol)
#        rotations.rotateMolecule(tempMol, tempMol, vect=vect, ang=rotation*angle, applyTo=movingAtms)
#        tempMol.rotate(rotMat, atms=toRot)
#        atms = dihed.getAtms()
#        vect1 = tempMol.getAtmWithIndex(atms[1]).getXYZ(vec=True) - tempMol.getAtmWithIndex(atms[0]).getXYZ(vec=True)
#        vect2 = tempMol.getAtmWithIndex(atms[3]).getXYZ(vec=True) - tempMol.getAtmWithIndex(atms[2]).getXYZ(vec=True)
#        dihedralAngle = rotations.vectAngle(vect1, vect2)
#        print(round(dihedralAngle*180/np.pi,3))
#        with open("OutputFiles/TEST/O1QJ_DH-%d_Rot-%d.pdb" %(count,int(round(rotAngle*180/np.pi,0))), "w") as fh:
#            fh.write(tempMol.printPDB())
#            fh.write("ENDMDL\n")
#        name = "rigid_O1QJ_%d_%d" %(count,int(round(rotAngle*180/np.pi,0)))
#        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
#            fh.write("%chk="+name+".chk\n")
#            fh.write("%mem=2gb\n")
#            fh.write("%NProcShared=1\n")
#            fh.write("#p b3lyp/6-31g(d)\n\n")
#            fh.write("SP "+name)
#            fh.write("\n\n0 1\n")
#            fh.write(tempMol.printChargeData())
#            fh.write("\n\n")
#        name = "fixedDihedral_O1QJ_%d_%d" %(count,int(round(rotAngle*180/np.pi,0)))
#        with open("OutputFiles/TEST/"+name+".gau", "w") as fh:
#            fh.write("%chk="+name+".chk\n")
#            fh.write("%mem=2gb\n")
#            fh.write("%NProcShared=1\n")
#            fh.write("#p b3lyp/6-31g(d) opt=modredundant\n\n")
#            fh.write("SP "+name)
#            fh.write("\n\n0 1\n")
#            fh.write(tempMol.printChargeData())
#            fh.write("\n"+str(dihed.getAtmAPos())+" "+str(dihed.getAtmBPos())+" "+str(dihed.getAtmCPos())+" "+str(dihed.getAtmDPos())+" F\n\n")
#        name = "fixedBonds_O1QJ_%d_%d" %(count,int(round(rotAngle*180/np.pi,0)))
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