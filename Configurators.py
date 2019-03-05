from ConfiguratorBase import Configurator
from Configuration import Configuration


import os

def runCommand(configuration, setting, command):
    if configuration.getSetting(setting)['verbose']:
        print(command)
    if not configuration.getSetting(setting)["dryRun"]:
        os.system(command)


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
    inPb = config.getInputFile(setting,'protein')
    chain = reduceSettings['chain']
    pythonBinary = config.getSetting('pythonBinary')

    bash_command = "{} {}/reduce.py {} {} --chain {} > /dev/null".format(pythonBinary,pathAttract, inPb, output, chain)
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
  
#create grid Files
def gridFunction(config,setting):
    output = config.getOutputFile(setting,'out')
    inPb = config.getInputFile(setting,'protein')
    alphabet = config.getInputFile(setting,'alphabet')
   
    attractBinPath = config.getSetting('attractBinPath')
    param = config.getSetting("attractParFile")
    bash_command = "{}/make-grid-omp {} {} 10.0 12.0 {}  --alphabet {}".format(attractBinPath,inPb,param, output, alphabet)
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

#cut free termini
def cutFunction(config,setting):
    setting = "cut"
    protein = config.getInputFile(setting,'protein')
    output = config.getOutputFile(setting,'out')
    
    #bash_command = "grep CA {}  > {}".format(protein, output)
    #runCommand(config, setting,bash_command)

#docking
def dockingFunction(config,setting):
    dof = config.getInputFile(setting,'dof')
    dockSettings = config.getSetting(setting)

    attractBinary = dockSettings['GPUattractBinary']
    receptor = config.getInputFile(setting,'receptor')
    ligand = config.getInputFile(setting,'ligand')
    gridRec = config.getInputFile(setting,'alphabetRec')
    gridLig = config.getInputFile(setting,'alphabetLig')

    alphabetRec = config.getInputFile(setting,'alphabetRec')
    alphabetLig = config.getInputFile(setting,'alphabetLig')
    modesRec = config.getInputFile(setting,'modesRec')
    modesLig = config.getInputFile(setting,'modesLig')
    result = config.getOutputFile(setting,'out')
    
    bash_command = "{} em --dof {} -p {} --alphabetrec {} --alphabetlig {} --gridrec {} --gridlig {} -r {} -l {} --numModesRec {} --numModesLig {} --modesr {} --modesl {} -d 0  > {}".format(attractBinary,dof,  dockSettings["attractParFile"], alphabetRec, alphabetLig, gridRec, gridLig, receptor, ligand, dockSettings["numModesRec"],dockSettings["numModesLig"], modesRec, modesLig, result  )
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
    
    bash_command = "{} {} {} {} {} --fix-receptor --modes {} --numModesRec {} --numModesLig  {} --vmax 1000  --rcut 50 --score > {}".format(attractBinary,dof,scoringSettings["attractParFile"],receptor,ligand,modes,scoringSettings["numModesRec"],scoringSettings["numModesLig"],result)
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
    
    bash_command = "{}/deredundant {} 2 --modes {} {} | {} {}/fill-deredundant.py /dev/stdin {} > {}".format(
        attractBinPath, inputDofFile, RedundantSetting['numModesRec'], RedundantSetting['numModesLig'], pythonBinary, toolPath,
        inputDofFile, output)
    runCommand(config, setting, bash_command)

#Demode
def demodeFunction(config,setting):
    toolPath =          config.getSetting(          'attractToolPath')
    inputDofFile =      config.getInputFile(setting,'inputDof')
    output =            config.getOutputFile(setting, 'out')
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
    toolPath =      config.getSetting(          'attractToolPath')
    inputDofFile =  config.getInputFile(setting,'inputDof')
    receptor =      config.getInputFile(setting,'receptor')
    receptorRef =   config.getInputFile(setting,'receptorRef')
    ligand =        config.getInputFile(setting,'ligand')
    lignadRef =     config.getInputFile(setting,'ligandRef')
    modes =         config.getInputFile(setting,'modes')
    output =        config.getOutputFile(setting,'out')
    pythonBinary = config.getSetting('pythonBinary')

    bash_command = "{} {}/irmsd.py {} {} {} {} {} --modes {} > {}".format(pythonBinary,toolPath, inputDofFile , receptor, receptorRef,ligand, lignadRef,modes,output)
    runCommand(config, setting, bash_command)

#RMSD
def RMSDFunction(config,setting):
    toolPath =          config.getSetting(          'attractToolPath')
    inputDofFile =      config.getInputFile(setting,'inputDof')
    receptorRef =       config.getInputFile(setting,'receptorRef')
    ligand =            config.getInputFile(setting,'ligand')
    ligandRef =         config.getInputFile(setting,'ligandRef')
    modes =             config.getInputFile(setting,'modes')
    output =            config.getOutputFile(setting,'out')
    pythonBinary = config.getSetting('pythonBinary')

    bash_command = "{} {}/lrmsd.py {} {} {} --modes {} --receptor {} > {}".format(pythonBinary,toolPath, inputDofFile , ligand, ligandRef,modes,receptorRef,output)
    runCommand(config, setting, bash_command)

#fnat
def FNATFunction(config,setting):
    toolPath =      config.getSetting(          'attractToolPath')
    inputDofFile =  config.getInputFile(setting,'inputDof')
    receptor =      config.getInputFile(setting,'receptor')
    receptorRef =   config.getInputFile(setting,'receptorRef')
    ligand =        config.getInputFile(setting,'ligand')
    ligandRef =     config.getInputFile(setting,'ligandRef')
    modes =         config.getInputFile(setting,'modes')
    output =        config.getOutputFile(setting,'out')
    pythonBinary = config.getSetting('pythonBinary')

    bash_command = "{} {}/fnat.py {} 5 {} {} {} {} --modes {} > {}".format(pythonBinary,toolPath, inputDofFile , receptor, receptorRef,ligand, ligandRef,modes,output)
    runCommand(config, setting, bash_command)

#collect pdbs from the resulting dofs
def CollectFunction(config,setting):
    attractBinPath = config.getSetting(         'attractBinPath')
    inputDofFile =  config.getInputFile(setting,'inputDof')
    receptor =      config.getInputFile(setting,'receptor')
    ligand =        config.getInputFile(setting,'ligand')
    output =        config.getOutputFile(setting,'out')

    bash_command = "{}/collect {} {} {} > {}".format(attractBinPath, inputDofFile, receptor, ligand, output)
    runCommand(config, setting, bash_command)

#specify configurators

#create coarsegrained representation of protein
ReduceConfigurator      = Configurator('reduce',reduceFunction)

#create all Atom Model of protein
AllAtomConfigurator     = Configurator('allAtom',allAtomFunction)

#create heavy atom model of protein
HeavyConfigurator       = Configurator('heavy',heavyFunction)

#create alphabet file of protein. Needed to calculate thr grid of the partner

AlphabetConfigurator    = Configurator('alphabet',alphabetFunction)

#create modes for protein
ModeConfigurator        = Configurator('modes',modeFunction)

#create grid for a given protein
GridConfigurator        = Configurator('grid',gridFunction)

#create dofs for a receptor ligand pair
DofConfigurator         = Configurator('dof',dofFunction)

#join two modefiles together into one file
JoinModesConfigurator   = Configurator('joinModes',joinModesFunction)

#create a CA- only represenetation of a protein
CAConfigurator          = Configurator('CA',CAFunction)

#cut free termini
CutConfigurator         = Configurator('cut',cutFunction)

#dock receptor and ligand
DockingConfigurator     = Configurator('docking',dockingFunction)

#score results
ScoringConfigurator     = Configurator('scoring',scoringFunction)

#sort results by energy
SortingConfigurator     = Configurator('sorting',SortingFunction)

#remove redundant results
DeRedundantConfigurator   = Configurator('deredundant',RedundantFunction)

#remoce modes from results
DemodeConfigurator      = Configurator('demode',demodeFunction)

#create a new file with the specified number of top results e.g. 50
TopConfigurator         = Configurator('top',topFunction)

#calculate the rmsd for a dof -fil input
RMSDConfigurator        = Configurator('rmsd',RMSDFunction)

#reate interface rmsd for a dof -file and the receptor/ligand
IRMSDConfigurator       = Configurator('irmsd',IRMSDFunction)

#calculate the native contact for a given dof file and receptor/ligand
FNATConfigurator        = Configurator('fnat',FNATFunction)

#create pdbfiles for given dofs and receptor/ligand
CollectConfigurator     = Configurator('collect',CollectFunction)
