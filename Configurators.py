from ConfiguratorBase import Configurator
from Configuration import Configuration
import os
import logging
import utils as utils
import json

#TODO superimposing 
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
    toolPath = config.getSetting(           'attractToolPath')
    docking = config.getInputFile(setting,  'dockingResult')
    scoring = config.getInputFile(setting,  'scoringResult')
    sortedOutput = config.getOutputFile(setting, 'out')
    pythonBinary = config.getSetting('pythonBinary')
    
    bash_command = "{} {}/fill-energies.py {} {} > {}".format(pythonBinary,toolPath, docking,scoring, sortedOutput)
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
        print(log)
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
    utils.superimposePdb(inputPdb, refPdb, superPdb)

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