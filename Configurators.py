from ConfiguratorBase import Configurator
from Configuration import Configuration
import os
import logging
import utils as utils
import json
import numpy as np

def runCommand(configuration, setting, command):
    if configuration.getSetting(setting)['verbose']:
        print(command)
    if not configuration.getSetting(setting)["dryRun"]:
        status = os.system(command)
        if status != 0: logging.warning("Failed to run " + command)

#create reduced proteins
def reduceFunction(config,setting):
    pathAttract = config.getSetting('attractToolPath')
    reduceSettings = config.getSetting(setting)
    output = config.getOutputFile(setting,'out')
    inPb = config.getInputFile(setting,'protein')
    chain = reduceSettings['chain']
    pythonBinary = config.getSetting('pythonBinary')
    bash_command = "{} {}/reduce.py {} {} --chain {} > /dev/null".format(pythonBinary,pathAttract, inPb, output, chain)
    runCommand(config,setting, bash_command)

#create reduced proteins
def allAtomFunction(config,setting):
    pathAttract = config.getSetting('attractToolPath')
    reduceSettings = config.getSetting(setting)
    output = config.getOutputFile(setting,'out')
    protein = config.getInputFile(setting,'protein')
    chain = reduceSettings['chain']
    pythonBinary = config.getSetting('pythonBinary')
    mapping = config.getOutputFile(setting,'mapping')
    attractBinPath = config.getSetting('attractBinPath')

    bash_command = "{} {}/../allatom/aareduce.py {} {} --chain {} --dumppatch --pdb2pqr > {}".format(pythonBinary,attractBinPath,protein,output,chain,mapping)
    runCommand(config,setting, bash_command)

#create reduced proteins
def heavyFunction(config,setting):
    pathAttract = config.getSetting('attractToolPath')
    reduceSettings = config.getSetting(setting)
    output = config.getOutputFile(setting,'out')
    protein = config.getInputFile(setting,'protein')
    chain = reduceSettings['chain']
    pythonBinary = config.getSetting('pythonBinary')
    attractBinPath = config.getSetting('attractBinPath')

    bash_command = "{} {}/../allatom/aareduce.py {} {} --heavy --chain {} --readpatch  > /dev/null".format(pythonBinary,attractBinPath,protein,output,chain)
    runCommand(config,setting, bash_command)

#create modes
def modeFunction(config,setting):
    toolPath = config.getSetting('attractToolPath')
    modeSettings = config.getSetting(setting)
    numModes = modeSettings['numModes']
    output = config.getOutputFile(setting,'out')
    inPb = config.getInputFile(setting,'protein')
    pythonBinary = config.getSetting('pythonBinary')

    bash_command = "{} {}/modes.py {} {} > {} ".format(pythonBinary,toolPath, inPb, numModes, output)
    runCommand(config, setting,bash_command)

#create modes evaluation 
def modeEvalFunction(config,setting):
    pdb_bound = config.getInputFile(setting,'protein_bound')
    pdb_unbound = config.getInputFile(setting,'protein_unbound')
    mode_file = config.getInputFile(setting,'mode_file')
    
    output = config.getOutputFile(setting,'out')

    if config.getSetting(setting)['verbose']:
        print("evaluating modes for" ,mode_file)

    if not config.getSetting(setting)["dryRun"]:
        bound_list = utils.readFileToList(pdb_bound)
        unbound_list = utils.readFileToList(pdb_unbound)

        unbound_residues = utils.getResidueFromPDBlines(unbound_list)
        
        bound_CA = utils.getCAOnlyFromPDBLines(bound_list)
        unbound_CA = utils.getCAOnlyFromPDBLines(unbound_list)

        bound_CA_pos = utils.getCoordinatesFromPDBlines(bound_CA)
        unbound_CA_pos = utils.getCoordinatesFromPDBlines(unbound_CA)


   
        modes = utils.read_modes(mode_file)
        cumulative_overlap = 0
        eval_dict = {}
        for modeIdx, mode in modes.items(): 
            ca_modes = utils.getCAModes(unbound_residues,mode['evec'])
            overlap = utils.getOverlap (unbound_CA_pos,bound_CA_pos, ca_modes)
            cumulative_overlap += overlap**2
            contributionCA = utils.getModeContribution(bound_CA_pos - unbound_CA_pos, ca_modes).tolist()
            norm = utils.getModeNorm(mode['evec'])
            contribution = contributionCA * norm
            magnitude = utils.getModeMagnitude(ca_modes)
            maximaIndices = utils.getIndexMaxima(magnitude)
            maxima = magnitude[maximaIndices]
            eval_dict[modeIdx] = { 'overlap':overlap, 'cum_overlap': np.sqrt(cumulative_overlap),'eigenvalue':mode['eval'],'norm':norm,'contribution':contribution,'contribution_ca':contributionCA,'maxima_indices':maximaIndices.tolist(), 'maxima_values':maxima.tolist() }

        utils.saveToJson(output, {'bound':pdb_bound, 'unbound':pdb_unbound, 'mode_file':mode_file, 'modes': eval_dict})


def createBoundModesFunction(config,setting):
    pdb_bound = config.getInputFile(setting,'protein_bound')
    pdb_unbound = config.getInputFile(setting,'protein_unbound')

    output = config.getOutputFile(setting,'out')
    if config.getSetting(setting)['verbose']:
        print("Create bound modes for " + pdb_unbound)
    if not config.getSetting(setting)["dryRun"]:
        bound_list = utils.readFileToList(pdb_bound)
        unbound_list = utils.readFileToList(pdb_unbound)

        bound_CA_pos = utils.getCoordinatesFromPDBlines(bound_list)
        unbound_CA_pos = utils.getCoordinatesFromPDBlines(unbound_list)

        pos_delta = bound_CA_pos - unbound_CA_pos
        norm = utils.getModeNorm(pos_delta)
    
        
        utils.writeModeFile(output, pos_delta.T[0],pos_delta.T[1],pos_delta.T[2], 1.0) #/norm**2


def manipulateModesFunction(config,setting):
    secondary_file = config.getInputFile(setting,'secondary')
    pdb = config.getInputFile(setting,'protein')
    mode_file = config.getInputFile(setting,'mode_file')

    output = config.getOutputFile(setting,'out')
    settings = config.getSetting(setting)
    if config.getSetting(setting)['verbose']:
        print("Manipulating modes for" ,mode_file)
    if not config.getSetting(setting)["dryRun"]:
        pdb_list = utils.readFileToList(pdb)    
        resIndices = utils.getResidueIndicesFromPDBLines(pdb_list)   

    
        modes = utils.read_modes(mode_file)
        sec = utils.readSecondaryStructure(secondary_file)

        for mode in modes.values():
            size = len(mode['evec'])
            for i in range(size):
                if sec[resIndices[i]] in settings['manipulate']:
                    mode['evec'][i] = np.zeros(3)

        utils.writeModeFileFromDict(modes,output)

#create alphatbetFiles
def alphabetFunction(config,setting):
    output = config.getOutputFile(setting,'out')
    inPb = config.getInputFile(setting,'protein')
    
    bash_command = "awk '{}' {} | sort -nu > {}".format('{print substr($0,58,2)}', inPb, output)
    runCommand(config,setting, bash_command)

#TODO substitute with real path
#create grid Files
def gridFunction(config,setting):
    output = config.getOutputFile(setting,'out')
    inPb = config.getInputFile(setting,'protein')
    alphabet = config.getInputFile(setting,'alphabet')
   
    attractBinPath = config.getSetting('attractBinPath')
    param = config.getSetting("attractParFile")
    bash_command = "/home/glenn/Downloads/attract_untouched/attract/bin/make-grid-omp {} {} 10.0 12.0 {}  --alphabet {} > /dev/null".format(inPb,param, output, alphabet)#attractBinPath
    runCommand(config, setting,bash_command)
       
#create dof File
def dofFunction(config,setting):
    dofSettings = config.getSetting(setting)

    receptor = config.getInputFile(setting,'receptor')
    ligand = config.getInputFile(setting,'ligand')
    output = config.getOutputFile(setting,'out')
    translate = config.getOutputFile(setting,'translate')
    pathout = config.getOutputFolder(setting)
    attractBinPath = config.getSetting('attractBinPath')

    systsearch = dofSettings['dofbinary']

    bash_command = "cp {}/../rotation.dat {}/{}".format(attractBinPath, pathout, "rotation.dat")
    runCommand(config, setting,bash_command)

    bash_command = "{}/translate {} {} > {}".format( attractBinPath, receptor, ligand, translate)
    runCommand(config, setting,bash_command)

    bash_command = "{} {}/rotation.dat  {} > {}".format( systsearch, pathout, translate, output )
    runCommand(config, setting,bash_command)
  
#create dof File
def joinModesFunction(config,setting):
    receptorModes = config.getInputFile(setting,'receptorModes')
    ligandModes = config.getInputFile(setting,'ligandModes')
    output = config.getOutputFile(setting,'out')
  
    bash_command = "cat /dev/null > {}".format(output)
    runCommand(config, setting,bash_command)
    bash_command = "cat {}  >> {}".format(receptorModes, output)
    runCommand(config, setting,bash_command)
    bash_command = "cat {}  >> {}".format(ligandModes, output)
    runCommand(config, setting,bash_command)

#create CA only pdb
def CAFunction(config,setting):
    protein = config.getInputFile(setting,'protein')
    output = config.getOutputFile(setting,'out')
    
    bash_command = "grep CA {}  > {}".format(protein, output)
    runCommand(config, setting,bash_command)

#docking
def dockingFunction(config,setting):
    dof = config.getInputFile(setting,'dof')
    dockSettings = config.getSetting(setting)

    attractBinary = dockSettings['GPUattractBinary']
    receptor = config.getInputFile(setting,'receptor')
    ligand = config.getInputFile(setting,'ligand')
    gridRec = config.getInputFile(setting,'gridRec')
    gridLig = config.getInputFile(setting,'gridLig')

    alphabetRec = config.getInputFile(setting,'alphabetRec')
    alphabetLig = config.getInputFile(setting,'alphabetLig')
    modesRec = config.getInputFile(setting,'modesRec')
    modesLig = config.getInputFile(setting,'modesLig')
    result = config.getOutputFile(setting,'out')
    devices = dockSettings["GPUdeviceIds"]
    d_string = ""
    if len(devices) > 0:
        for d in devices:
            d_string += " -d {}".format(d)

    mode_string = ""
    if dockSettings["numModesRec"] > 0 or dockSettings["numModesLig"] > 0 :
        mode_string = "--numModesRec {} --numModesLig {} --modesr {} --modesl {} --evscale {}".format(dockSettings["numModesRec"],dockSettings["numModesLig"], modesRec, modesLig, dockSettings['evScale'])
    bash_command = "{} em --dof {} -p {} --alphabetrec {} --alphabetlig {} --gridrec {} --gridlig {} -r {} -l {} {} {}  > {}".format(attractBinary,dof,  dockSettings["attractParFile"], alphabetRec, alphabetLig, gridRec, gridLig, receptor, ligand,mode_string ,d_string, result  )
    runCommand(config, setting, bash_command)

#scoring
def scoringFunction(config,setting):
    scoringSettings = config.getSetting(setting)

    dof = config.getInputFile(setting,      'dof')
    attractBinary = scoringSettings[        'attractBinary']
    receptor = config.getInputFile(setting, 'receptor')
    ligand = config.getInputFile(setting,   'ligand')
    modes = config.getInputFile(setting,    'joinedModes')
    result = config.getOutputFile(setting,  'out')
    
    mode_string = ""
    if scoringSettings["numModesRec"] > 0 or scoringSettings["numModesLig"] > 0 :
        mode_string = "--numModesRec {} --numModesLig {} --modes {} --evscale {}".format(scoringSettings["numModesRec"],scoringSettings["numModesLig"], modes, scoringSettings['evScale'])
    rcut = scoringSettings["cutoff"]
    rcutstr = "--rcut {}".format(rcut)  if rcut != None and rcut > 0  else  ""

    bash_command = "{} {} {} {} {} --fix-receptor --score {} {}  > {}".format(attractBinary,dof,scoringSettings["attractParFile"],receptor,ligand,rcutstr,mode_string,result)
    runCommand(config, setting, bash_command)

#Sorting
def SortingFunction(config,setting):
    toolPath = config.getSetting(                   'attractToolPath')
    pythonBinary = config.getSetting(               'pythonBinary')

    input_dof = config.getInputFile(setting,        'input_dof')
    sorted_output = config.getOutputFile(setting,   'out')
    
    bash_command = "{} {}/sort.py {} > {}".format(pythonBinary, toolPath, input_dof, sorted_output)
    runCommand(config, setting, bash_command)


#FillEnergies
def FillEnergyFunction(config,setting):
    toolPath = config.getSetting(                   'attractToolPath')
    pythonBinary = config.getSetting(               'pythonBinary')


    docking = config.getInputFile(setting,          'dockingResult')
    scoring = config.getInputFile(setting,          'scoringResult')
    filled_output = config.getOutputFile(setting,   'out')
    
    bash_command = "{} {}/fill-energies.py {} {} > {}".format(pythonBinary,toolPath, docking,scoring, filled_output)
    runCommand(config, setting, bash_command)

#Redundant
def RedundantFunction(config,setting):
    RedundantSetting = config.getSetting(       setting)
    attractBinPath = config.getSetting(         'attractBinPath')
    inputDofFile = config.getInputFile(setting, 'inputDof')
    toolPath = config.getSetting(               'attractToolPath')
    output = config.getOutputFile(setting,"out")
    pythonBinary = config.getSetting('pythonBinary')

    mode_string = ""
    if RedundantSetting['numModesRec'] > 0 or RedundantSetting['numModesLig'] > 0:
        mode_string = "--modes {} {} ".format( RedundantSetting['numModesRec'], RedundantSetting['numModesLig'])

    bash_command = "{}/deredundant {} 2 {} | {} {}/fill-deredundant.py /dev/stdin {} > {}".format(
        attractBinPath, inputDofFile, mode_string, pythonBinary, toolPath,
        inputDofFile, output)
    runCommand(config, setting, bash_command)

#Demode
def demodeFunction(config,setting):
    toolPath =          config.getSetting(              'attractToolPath')
    inputDofFile =      config.getInputFile(setting,    'inputDof')
    output =            config.getOutputFile(setting,   'out')
    pythonBinary = config.getSetting('pythonBinary')
    
    bash_command = "{} {}/demode.py {} > {}".format(pythonBinary, toolPath, inputDofFile, output)
    runCommand(config, setting, bash_command)

#Top
def topFunction(config,setting):
    TopSetting =    config.getSetting(        setting)
    toolPath =      config.getSetting(              'attractToolPath')
    inputDofFile =  config.getInputFile(setting,    'inputDof')
    output =        config.getOutputFile(setting,   'out')

    bash_command = "{}/top {} {} > {}".format(toolPath, inputDofFile, TopSetting['numDof'], output)
    runCommand(config, setting, bash_command)

#IRMSD
def IRMSDFunction(config,setting):
    irmsdSetting =    config.getSetting(        setting)
    toolPath =      config.getSetting(          'attractBinPath')
    inputDofFile =  config.getInputFile(setting,'inputDof')
    receptor =      config.getInputFile(setting,'receptor')
    receptorRef =   config.getInputFile(setting,'receptorRef')
    ligand =        config.getInputFile(setting,'ligand')
    lignadRef =     config.getInputFile(setting,'ligandRef')
    modes =         config.getInputFile(setting,'modes')
    output =        config.getOutputFile(setting,'out')
    pythonBinary = config.getSetting(            'pythonBinary')

    mode_string = ""
    #TODO remoce
    if irmsdSetting['numModesRec'] > 0 or irmsdSetting['numModesLig'] > 0:
        mode_string = "--modes {}".format(modes)

    bash_command = "{} {}/irmsd.py {} {} {} {} {}  {} > {}".format(pythonBinary,toolPath, inputDofFile , receptor, receptorRef,ligand, lignadRef,mode_string,output)
    runCommand(config, setting, bash_command)

#RMSD
def RMSDFunction(config,setting):
    rmsdSetting =    config.getSetting(        setting)
    toolPath =          config.getSetting(          'attractBinPath')
    inputDofFile =      config.getInputFile(setting,'inputDof')
    receptorRef =       config.getInputFile(setting,'receptorRef')
    ligand =            config.getInputFile(setting,'ligand')
    ligandRef =         config.getInputFile(setting,'ligandRef')
    modes =             config.getInputFile(setting,'modes')
    output =            config.getOutputFile(setting,'out')
    pythonBinary = config.getSetting(                'pythonBinary')

    mode_string = ""
    if rmsdSetting['numModesRec'] > 0  or rmsdSetting['numModesLig'] > 0 :
        mode_string = " --modes {} ".format(modes)

    bash_command = "{} {}/lrmsd.py {} {} {} {} --receptor {} > {}".format(pythonBinary,toolPath, inputDofFile , ligand, ligandRef,mode_string,receptorRef,output)
    runCommand(config, setting, bash_command)

#fnat
def FNATFunction(config,setting):
    fnatSetting =   config.getSetting(setting)
    toolPath =      config.getSetting(          'attractBinPath')
    inputDofFile =  config.getInputFile(setting,'inputDof')
    receptor =      config.getInputFile(setting,'receptor')
    receptorRef =   config.getInputFile(setting,'receptorRef')
    ligand =        config.getInputFile(setting,'ligand')
    ligandRef =     config.getInputFile(setting,'ligandRef')
    modes =         config.getInputFile(setting,'modes')
    output =        config.getOutputFile(setting,'out')
    pythonBinary = config.getSetting(            'pythonBinary')

    mode_string = ""
    #TODO remoce

    if fnatSetting['numModesRec'] > 0  or fnatSetting['numModesLig'] > 0 :
        mode_string = " --modes {} ".format(modes)

    bash_command = "{} {}/fnat.py {} 5 {} {} {} {} {} > {}".format(pythonBinary,toolPath, inputDofFile , receptor, receptorRef,ligand, ligandRef,mode_string,output)
    runCommand(config, setting, bash_command)

#collect pdbs from the resulting dofs
def CollectFunction(config,setting):
    """Saves the best scored structures to a pdb file. The number of structures has to be specified"""
    collectSetting = config.getSetting(setting)
    attractBinPath = config.getSetting(          'attractBinPath')
    inputDofFile =  config.getInputFile(setting, 'inputDof')
    receptor =      config.getInputFile(setting, 'receptor')
    ligand =        config.getInputFile(setting, 'ligand')
    modes =         config.getInputFile(setting, 'modes')
    output =        config.getOutputFile(setting,'out')

    mode_string =""
    if collectSetting['numModesRec'] > 0 or collectSetting["numModesLig"] > 0:
        mode_string =" --modes {}".format(modes)
    bash_command = "{}/collect {} {} {} {} > {}".format(attractBinPath, inputDofFile, receptor, ligand,mode_string, output)
    runCommand(config, setting, bash_command)

#save the configuration
def saveSettings(config, setting):
    """ saves the configuration to the specified filename"""
    configFilename = config.getOutputFile(setting, 'out')
    if config.getSetting(setting)['verbose']:
        print("saving configFile to " + configFilename)
    if not config.getSetting(setting)["dryRun"]:
        config.save(configFilename)

def FindTermini(config, setting):
    cutSetting = config.getSetting(setting)
    inputPdb = config.getInputFile(setting,                 "pdb")
    looseTerminiLog = config.getOutputFile(setting,         "out")
    cutoff =                        cutSetting[             'cutoff']

    if config.getSetting(setting)['verbose']:
        print("Find Termini from  " + inputPdb)
    if not config.getSetting(setting)["dryRun"]:
        # log = utils.findAndCutLooseTermini(inputPdb, cutPdb, cutoff)
        log = utils.FindLooseTermini(inputPdb, cutoff=cutoff)
        log['cutoff'] = cutoff
        utils.saveToJson(looseTerminiLog, log)

def cutTermini(config, setting):
    cutSetting = config.getSetting(setting)
    inputPdb = config.getInputFile(setting,     "pdb")
    cutlog = config.getInputFile(setting,       "cutlog")

    cutPdb = config.getOutputFile(setting,      "out")

    if config.getSetting(setting)['verbose']:
        print("Cut Termini from  " + inputPdb + " and output to " + cutPdb)
    if not config.getSetting(setting)["dryRun"]:
        log = utils.loadFromJson(cutlog)
        residues = log['looseTerminiFront'] + log['looseTerminiBack']
        pdblines = utils.readFileToList(inputPdb)
        utils.cutTerminiAndWriteToPdb(residues,pdblines, cutPdb)


#creates a secondary structure file from an input pdb file
def CreateSecondary(config, setting):
    inputPdb = config.getInputFile(setting,         "pdb")
    secdondaryFile = config.getOutputFile(setting,  "out")

    bash_command ="stride {} > {}".format(inputPdb, secdondaryFile)
    runCommand(config, setting, bash_command)

def GetInterface(config, setting):
    #receptor = config.getInputFile(setting, 'receptor')
    #ligand = config.getInputFile(setting,'ligand')
    pdb = config.getInputFile(setting,              'pdb')
    interfaceFile = config.getOutputFile(setting,   'out')
    cutoff = config.getSetting(setting)['cutoff']

    if config.getSetting(setting)['verbose']:
        print("Get interface from pdb " + pdb )
    if not config.getSetting(setting)["dryRun"]:
        structures = utils.parseBIOPdbToStructure(pdb)
        interfaces = []
        for struct in structures:
            receptor =  struct['A']       
            ligand =    struct['B']                
        
            contactResiduesRec, contactResiduesLig = utils.getInterfaceResidues(receptor,ligand,cutoff )
            recinterfaceResidues = utils.getResidueIds(contactResiduesRec)
            liginterfaceResidues = utils.getResidueIds(contactResiduesLig)

            interfaces.append({'model': struct.id,"recInterfaceResidues": recinterfaceResidues,"ligInterfaceResidues": liginterfaceResidues})

        utils.saveToJson(interfaceFile, {'file': pdb, 'cutoff': cutoff,'interfaces': interfaces})

def SuperimposeStructures(config, setting):
    inputPdb = config.getInputFile(setting, "pdb")
    refPdb = config.getInputFile(setting, "refpdb")
    superPdb = config.getOutputFile(setting, "out")
    
    if config.getSetting(setting)['verbose']:
        print("Superimpose pdb " + inputPdb )
    if not config.getSetting(setting)["dryRun"]:
        utils.superimposePdb(inputPdb, refPdb, superPdb)



def evaluateModeDOFS(config, setting):
    input_dof_file = config.getInputFile(setting, "input_dof")
    mode_evaluation_rec = config.getInputFile(setting, "mode_evaluation_rec")
    mode_evaluation_lig = config.getInputFile(setting, "mode_evaluation_rec")
    output = config.getOutputFile(setting,"out")
    
    dof_eval_settings = config.getSetting(setting)
    num_eval = dof_eval_settings['num_eval']
    numModesRec = dof_eval_settings['numModesRec']
    numModesLig = dof_eval_settings['numModesLig']

    if config.getSetting(setting)['verbose']:
        print("evaluating dofs for" , input_dof_file)
    if not config.getSetting(setting)["dryRun"]:
        dof_dict = utils.read_Dof(input_dof_file)
        sorted_keys = np.sort(np.asarray(list(dof_dict.keys()),dtype=np.int))

        contributions_rec = {}
        for key, val in json.load(open(mode_evaluation_rec, 'r'))['modes'].items():
            contributions_rec[int(key)] = val['contribution']
        
        contributions_lig = {}
        for key, val in json.load(open(mode_evaluation_lig, 'r'))['modes'].items():
            contributions_lig[int(key)] = val['contribution']

        result = {}
        for key in sorted_keys[:num_eval]:
            dof = dof_dict[key]
            modes_rec = dof['rec'][6:]
            modes_lig = dof['lig'][6:]
            rec = {}
            lig = {}
            for i,mode in enumerate(modes_rec):
                rec[str(i+1)] = {'ratio':np.float64(mode)/contributions_rec[i+1] -1, 'dof':mode, 'mode':contributions_rec[i+1]}
            for i,mode in enumerate(modes_lig):
                lig[str(i+1)] = {'ratio':np.float64(mode)/contributions_lig[i+1] -1, 'dof':mode,  'mode':contributions_lig[i+1]}
            result[str(key)] = {'rec': rec, 'lig':lig}
        utils.saveToJson(filename =output, data=result)
    
def createTestDof(config, setting):
    output = config.getOutputFile(setting, 'out')
    testDof = """#pivot auto\n#centered receptor: false\n#centered ligands: false\n#1\n           0           0           0           0           0           0\n    0           0           0           0           0           0"""
    
    testDof = """#pivot auto
#centered receptor: false
#centered ligands: false
#1
           0           0           0           0           0           0
           0           0           0           0           0           0\r\n"""
    if not config.getSetting(setting)["dryRun"]:
        with open(output, 'w+') as f:
            f.write(testDof)
        

def prunePDB(config, setting):

    inputPdb = config.getInputFile(setting, 'pdb')
    output = config.getOutputFile(setting, 'out')


    if config.getSetting(setting)['verbose']:
        print("pruning pdb: " , inputPdb)
    if not config.getSetting(setting)["dryRun"]:
        currId = ""
        pdb_lines = utils.readFileToList(inputPdb)

        for i in range(len(pdb_lines)):
            line = pdb_lines[i]
            resId = utils.pdbentry('residueIdChain',line)
            if currId != resId:
                currId = resId
                if i-3 > 0 and utils.pdbentry('atom',pdb_lines[i-3] ).strip() == 'O':
                    lineRestStart = pdb_lines[i-2]
                    lineRestEnd = pdb_lines[i-1]
                    x = utils.pdbentry('posX',lineRestStart)
                    y = utils.pdbentry('posY',lineRestStart)
                    z = utils.pdbentry('posZ',lineRestStart)

                    lineRestEnd = utils.setpdbentry('posX',lineRestEnd,x )
                    lineRestEnd = utils.setpdbentry('posY',lineRestEnd,y )
                    lineRestEnd = utils.setpdbentry('posZ',lineRestEnd,z )
                    pdb_lines[i-1] = lineRestEnd

        with open(output,'w+') as f:
            for line in pdb_lines:
                f.write(line)



########################specify configurators########################

#create coarsegrained representation of protein
ReduceConfigurator      = Configurator('reduce', reduceFunction)

#create all Atom Model of protein
AllAtomConfigurator     = Configurator('allAtom', allAtomFunction)

#create heavy atom model of protein
HeavyConfigurator       = Configurator('heavy', heavyFunction)

#create alphabet file of protein. Needed to calculate thr grid of the partner

AlphabetConfigurator    = Configurator('alphabet', alphabetFunction)

#create modes for protein
ModeConfigurator        = Configurator('modes', modeFunction)

#create grid for a given protein
GridConfigurator        = Configurator('grid', gridFunction)

#create dofs for a receptor ligand pair
DofConfigurator         = Configurator('dof', dofFunction)

#join two modefiles together into one file
JoinModesConfigurator   = Configurator('joinModes', joinModesFunction)

#create a CA- only represenetation of a protein
CAConfigurator          = Configurator('CA', CAFunction)

#dock receptor and ligand
DockingConfigurator     = Configurator('docking', dockingFunction)

#score results
ScoringConfigurator     = Configurator('scoring', scoringFunction)

#create a file with filled in energies using the docking and scoring results
FillEnergyConfigurator     = Configurator('fill_energy', FillEnergyFunction)

#sort results by energy
SortingConfigurator     = Configurator('sorting', SortingFunction)

#removes redundant results 
DeRedundantConfigurator = Configurator('deredundant', RedundantFunction)

#remove modes from results
DemodeConfigurator      = Configurator('demode', demodeFunction)

#create a new file with the specified number of top results e.g. 50
TopConfigurator         = Configurator('top',topFunction)

#calculate the rmsd for a dof -fil input
RMSDConfigurator        = Configurator('rmsd', RMSDFunction)

#reate interface rmsd for a dof -file and the receptor/ligand
IRMSDConfigurator       = Configurator('irmsd', IRMSDFunction)

#calculate the native contact for a given dof file and receptor/ligand
FNATConfigurator        = Configurator('fnat', FNATFunction)

#create pdbfiles for given dofs and receptor/ligand
CollectConfigurator     = Configurator('collect', CollectFunction)

#allows to save the configuration to file
SaveSettingConfigurator = Configurator('saveSettings', saveSettings )

#cuts of all residues of a protein that is have no second next neighbor below a cutoff (~5.5Angstroem)
FindTerminiConfigurator  = Configurator('findtermini', FindTermini)

#cuts of all residues of a protein that is have no second next neighbor below a cutoff (~5.5Angstroem)
CutTerminiConfigurator  = Configurator('cut', cutTermini)

#creates a secondary structure file of the inputpdb according to stride
SecondaryConfigurator   = Configurator('secondary', CreateSecondary)

#analyses the interface of two structures
InterfaceCondigurator   = Configurator('interface', GetInterface)

#superimposes a structe to a certain target and saves it to a pdb
SuperimposeConfigurator = Configurator('superimpose', SuperimposeStructures)














configuratorFunctions = {
    'reduce':       reduceFunction,
    'allAtom':      allAtomFunction,
    'heavy':        heavyFunction,
    'alphabet':     alphabetFunction,
    'modes':        modeFunction,
    'grid':         gridFunction,
    'dof':          dofFunction,
    'joinModes':    joinModesFunction,
    'CA':           CAFunction,
    'docking':      dockingFunction,
    'scoring':      scoringFunction,
    'fill_energy':  FillEnergyFunction,
    'sorting':      SortingFunction,
    'deredundant':  RedundantFunction,
    'demode':       demodeFunction,
    'top':          topFunction,
    'rmsd':         RMSDFunction,
    'irmsd':        IRMSDFunction,
    'fnat':         FNATFunction,
    'collect':      CollectFunction,
    'saveSettings': saveSettings ,
    'findtermini':  FindTermini,
    'cut':          cutTermini,
    'secondary':    CreateSecondary,
    'interface':    GetInterface,
    'superimpose':  SuperimposeStructures,
    'mode_evaluation':modeEvalFunction,
    'bound_mode': createBoundModesFunction,
    'dof_evaluation':evaluateModeDOFS,
    'dof_test':     createTestDof,
    'mode_manipulate':manipulateModesFunction,
    'prune':        prunePDB
    }
