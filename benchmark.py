from Configurators import *
from Settings import *
from Configuration import Configuration
from Pipeline import *
from ConfiguratorBase import createLoggingFile
import math
from datetime import datetime
import sys


F = False
T = True

protType =      "unbound"
protTypeRef =   "refe"
pathfile =      "paths_local.json"
paths = json.load(open(pathfile,'r'))
basePath = paths["basePath"]
utils.strideBinary = paths["strideBinary"]
utils.rmscaBinary = paths["rmscaBinary"]
protlist_file = basePath + "/protlist"
dry = False
verbose = True
overwrite = False
save_consoleOutput = T

configs = [
{"numModesRec":1,"numModesLig":1, "scale": 1,'prune':F,'cut':F,'modesOnly':F,'bound':F,'boundModes':F,'singleDof':F,'manipulateModes':F,'rigidStart':F,'overwrite': F,'extension':'_hin99_7'},
{"numModesRec":1,"numModesLig":1, "scale": 1,'prune':F,'cut':T,'modesOnly':F,'bound':F,'boundModes':F,'singleDof':F,'manipulateModes':F,'rigidStart':F,'overwrite': F,'extension':'_cut1'},
{"numModesRec":1,"numModesLig":1, "scale": 1,'prune':T,'cut':F,'modesOnly':F,'bound':F,'boundModes':F,'singleDof':F,'manipulateModes':F,'rigidStart':F,'overwrite': F,'extension':'_pruned1'},
{"numModesRec":1,"numModesLig":1, "scale": 1,'prune':T,'cut':T,'modesOnly':F,'bound':F,'boundModes':F,'singleDof':F,'manipulateModes':F,'rigidStart':F,'overwrite': F,'extension':'_cut_pruned1'},
]

with open(basePath + "/configs_{}.json".format(str(datetime.now())), "w") as write_file:
    json.dump({'configs':configs,'description': "recalculate the mode evaluation for the receptor modes (bound)."}, write_file, indent=4)

deleted_proteins = ['1I9R','1BJ1', '1FSK', '1IQD', '1K4C', '1KXQ', '1NCA', '1NSN', '1QFW', '2HMI', '2JEL', '9QFW','1DE4','4FQI','4GAM','4GXU','1EXB','4FQI', '1EER','4H30' ,'2OOB','1N2C']
proteins = [line.split()[0].strip() for line in open(protlist_file ,'r').readlines() ]
try:
    for dele in deleted_proteins:
        proteins.remove(dele)
except:
    pass

#proteins =[ 
#'7CEI',
#'1AHW',
#'1AK4',
#'1AKJ',
#'1ATN',
#'1AVX',
#'1B6C',
#'3MXW'
#]

#create Modes, coarse grain Proteins and alphabet files
singleConfiguration_01 = [{'conf': [ 
        {'setting': 'superimpose' ,         'configurator': 'superimpose'},
        {'setting': 'findtermini',          'configurator': 'findtermini'},
        {'setting': 'cut',                  'configurator': 'cut'},
        {'setting': 'allAtom',              'configurator': 'allAtom'},
        {'setting': 'reduce',               'configurator': 'reduce'},
        {'setting': 'heavy',                'configurator': 'heavy'},
        {'setting': 'prune',                'configurator': 'prune'},
        {'setting': 'modes',                'configurator': 'modes'},
        {'setting': 'modes_heavy',          'configurator': 'modes'},
        {'setting': 'alphabet',             'configurator': 'alphabet'},
        {'setting': 'secondary',            'configurator': 'secondary'},
        {'setting': 'mode_manipulate',      'configurator': 'mode_manipulate'},
        {'setting': 'mode_manipulate_heavy','configurator': 'mode_manipulate'},
        {'setting': 'saveSettings',         'configurator': 'saveSettings'},
        {'setting': 'protein_evaluation',   'configurator': 'protein_evaluation'},
        ],'numThreads': 8}]

singleConfigurationRef = [{'conf':[
        { 'setting': 'cut' ,                'configurator': 'cut'},
        #{ 'setting': 'superimpose' ,        'configurator': 'superimpose'},
        { 'setting': 'allAtom' ,            'configurator': 'allAtom'},
        { 'setting': 'reduce' ,             'configurator': 'reduce'},
        { 'setting': 'heavy' ,              'configurator':  'heavy'},
        { 'setting': 'saveSettings',        'configurator': 'saveSettings'} ,
        {'setting': 'prune',                'configurator': 'prune'},
        ] , 'numThreads': 8}
]
#create Grids. this has to be a seperate step since the alphabet file of the partner is needed
singleConfiguration_02 = [ {'conf':[
        {'setting': 'bound_mode',           'configurator': 'bound_mode'},
        {'setting': 'bound_mode_heavy',     'configurator': 'bound_mode'},
        {'setting': 'mode_evaluation',      'configurator': 'mode_evaluation'},
        { 'setting': 'grid',                'configurator':'grid'}
        ],'numThreads': 2}
]
#All files are created which needed receptor and ligand
pairConfiguration = [{'conf':[
       { 'setting' : 'dof_test',           'configurator': 'dof_test'      },
        { 'setting': 'dof',                 'configurator': 'dof'           },
        { 'setting': 'joinModes',           'configurator': 'joinModes'     },
        { 'setting': 'joinModes_heavy',     'configurator': 'joinModes'     } 
        ],'numThreads':8}
]
#perform docking and scoring as well analysis
runConfiguration = [ {'conf':
         { 'setting':'docking',              'configurator':'docking'},
     'numThreads': 1},
     {'conf':
         { 'setting':'scoring',             'configurator':'scoring'},
      'numThreads':2},
    {'conf':[
        { 'setting': 'fill_energy',         'configurator': 'fill_energy'   },
        { 'setting': 'sorting',             'configurator': 'sorting'       },
        { 'setting': 'deredundant',         'configurator': 'deredundant'   },
        { 'setting': 'top',                 'configurator': 'top'           },
        { 'setting': 'demode',              'configurator': 'demode'        },
        { 'setting': 'irmsd',               'configurator': 'irmsd'         },
        ###{ 'setting': 'irmsd_nomodes',     'configurator': 'irmsd'         },
        { 'setting': 'rmsd',                'configurator': 'rmsd'          },
        ## { 'setting': 'rmsd_unsorted',      'configurator': 'rmsd'          },
        ##{ 'setting': 'rmsd_nomodes',      'configurator': 'rmsd'          },
        { 'setting': 'fnat',                'configurator': 'fnat'          },
        ###{ 'setting': 'fnat_nomodes',      'configurator': 'fnat'          },
        { 'setting': 'saveSettings',        'configurator': 'saveSettings'  },
        { 'setting': 'collect'  ,           'configurator': 'collect'       },
         { 'setting': 'dof_evaluation',      'configurator': 'dof_evaluation'},
        { 'setting': 'dof_extractionLow',   'configurator': 'dof_extraction'},
        { 'setting': 'dof_extractionHigh',  'configurator': 'dof_extraction'},
         { 'setting': 'dof_evaluation_low',  'configurator': 'dof_evaluation'},
         { 'setting': 'dof_evaluation_high', 'configurator': 'dof_evaluation'},
        { 'setting': 'collect_Low'  ,       'configurator': 'collect'       },
        { 'setting': 'collect_High'  ,      'configurator': 'collect'       },
        ##{ 'setting': 'interface',          'configurator':'interface'      },
        { 'setting': 'interface_low',       'configurator':'interface'      },
        { 'setting': 'interface_high',      'configurator':'interface'      }
    ],'numThreads': 4}
    ]

bufferSize = 4*len(proteins) * len(configs)
prefix = "00_"
for i,c in enumerate (configs[:1]):
    prefix += "_BM{}_mr{}ml{}s{}_{}".format(i, c['numModesRec'], c['numModesLig'], c['scale'],c['extension'])


createLoggingFile(basePath+"/{}_loggingFile_{}.log".format(prefix,str(datetime.now().date())))
consoleOutputFile = basePath+"/{}_ConsoleOutput_{}.log".format(prefix,str(datetime.now().date()))
if save_consoleOutput:
    sys.stdout = open(consoleOutputFile, 'w')
#proteins = ['3MXW']

num = len(proteins) * len(configs)
#configure receptor and ligand in two steps
pipelineRec         = createPipeline( bufferSize, singleConfiguration_01,  num)
pipelineLig         = createPipeline( bufferSize, singleConfiguration_01, num)
#configure grids and joined Modes
pipelineRec2        = createPipeline( bufferSize, singleConfiguration_02, num)
pipelineLig2        = createPipeline( bufferSize, singleConfiguration_02, num)
#configure referenceStructures
pipelineRecRef      = createPipeline( bufferSize, singleConfigurationRef,  num)
pipelineLigRef      = createPipeline( bufferSize, singleConfigurationRef,  num)
#runpairconfiguration
pipelinePairConfig  = createPipeline( bufferSize, pairConfiguration,  num)
#run docking, scoring and analysis
pipelinePairRun     = createPipeline( bufferSize, runConfiguration,  num)


for i,config in enumerate(configs):

    use_cut = config['cut']
    use_bound = config['bound']
    use_boundModes = config['boundModes']
    use_singleDof = config['singleDof']
    use_pruning = config['prune']
    use_manipulated = config['manipulateModes']
    extension = config['extension']
    modeType = "hin99"
    protType = 'unbound'
    protTypeRef = 'refe'
    if use_bound:
        protType = 'refe'
    if use_cut:
        protType = "unbound-cut"
        protTypeRef = "refe-cut"
    if use_boundModes:
        modeType = "bound"
    if use_manipulated:
        modeType = "hin99"

    numModesRec = config['numModesRec']
    numModesLig = config['numModesLig']
    numModes = 20
    if use_boundModes :
        numModes = 1

    
    scale = config['scale']

    frac, whole = math.modf(scale)
    bm = "02_bm_dG_mr{}_ml{}_s{}p{}_sO_c50_mr{}_ml{}_s{}p{}{}".format(numModesRec, numModesLig, int(whole),'{:6f}'.format(frac)[2:],numModesRec, numModesLig, int(whole),'{:6f}'.format(frac)[2:],extension)
    print(" START NEW BENCHMARK NUMBER {} mr{} ml{} scale{} name {}".format(i,numModesRec,numModesLig, scale,bm))
        
    print("Create configrations".upper())

    logging.warning("\n\n---------------------CREATE CONFIGURATIONS--------------------------\n".upper())

    for protein in proteins:
        receptorConfig  = Configuration(getDefaultSingleSetting(    protein=protein, chain="A", protType = protType, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite, 
            pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"],
            modeType = modeType))
        
        ligandConfig    = Configuration(getDefaultSingleSetting(    protein=protein, chain="B", protType = protType, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite,
            pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"],
            modeType = modeType))
        
        receptorRefConfig = Configuration(getDefaultSingleSetting(  protein=protein, chain="A", protType = protTypeRef, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite,
            pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"],
            modeType = modeType))
        
        ligandRefConfig = Configuration(getDefaultSingleSetting(    protein=protein, chain="B", protType = protTypeRef, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite,
            pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"],
            modeType = modeType))

        pairConfig      = Configuration(getDefaultPairSetting(benchmarkName = bm,protein = protein, protType = protType, protTypeRef = protTypeRef, 
            numModesRec = numModesRec,
            numModesLig = numModesLig,
            basePath = basePath + "/{}".format(protein), 
            dry=dry, verbose = verbose, overwrite = overwrite,
            attractBinary   =  paths["attractBinary"] , 
            attractBinPath  =  paths["attractBinPath"], 
            attractParFile  =  paths["attractParFile"],
            dofBinary       =  paths["dofBinary"],
            pythonBinary    =  paths["python2Binary"],
            attractToolPath =  paths["attractToolPath"],
            attractBinaryGPU = paths["attractBinaryGPU"],
            deviceIds = [0,1], evScale=scale,
            modeType = modeType))

        #receptorRefConfig.settings['allAtom']['in']['protein'] =    'superimpose'
        #ligandRefConfig.settings['allAtom']['in']['protein'] =      'superimpose'


        receptorConfig.settings['cut']['in']['pdb'] =    'superimpose'
        ligandConfig.settings['cut']['in']['pdb']   =    'superimpose'

        receptorConfig.settings['allAtom']['in']['protein'] =    'superimpose'
        ligandConfig.settings['allAtom']['in']['protein']   =    'superimpose'

        # receptorRefConfig.files['refpdb']   = receptorConfig.files['pdb'] 
        # ligandRefConfig.files['refpdb']     =   ligandConfig.files['pdb'] 

        receptorConfig.files['refpdb']   = receptorRefConfig.files['pdb'] 
        ligandConfig.files['refpdb']     =   ligandRefConfig.files['pdb'] 


        receptorRefConfig.files['cutlog']   = receptorConfig.files['cutlog'] 
        ligandRefConfig.files['cutlog']     = ligandConfig.files['cutlog'] 

        receptorConfig.files['alphabetPartner'] = ligandConfig.files['alphabet'] 
        ligandConfig.files['alphabetPartner'] = receptorConfig.files['alphabet'] 

        receptorConfig.files['partner_bound'] = receptorRefConfig.files['reduce']
        ligandConfig.files['partner_bound'] = ligandRefConfig.files['reduce']

        ligandConfig.files['partner_bound_heavy'] = ligandRefConfig.files['heavy']
        receptorConfig.files['partner_bound_heavy'] = receptorRefConfig.files['heavy']

        pairFiles = pairConfig.files
        pairFiles['mode_evaluation_rec'] = receptorConfig.files['mode_evaluation']
        pairFiles['mode_evaluation_lig'] = ligandConfig.files['mode_evaluation']

        pairFiles['receptor'] =             receptorConfig.files['reduce']
        pairFiles['receptorSec'] =             receptorConfig.files['secondary']
        pairFiles['receptor_heavy'] =       receptorConfig.files['heavy']

        pairFiles['modesRec'] =     receptorConfig.files['modes']
        pairFiles['modesRec_heavy'] =     receptorConfig.files['modes_heavy']

        pairFiles['gridRec'] =      receptorConfig.files['grid']
        pairFiles['alphabetRec'] =  receptorConfig.files['alphabetPartner']

        pairFiles['ligand'] =       ligandConfig.files['reduce']
        pairFiles['ligandSec'] =             receptorConfig.files['secondary']

        pairFiles['ligand_heavy'] = ligandConfig.files['heavy']
        pairFiles['modesLig'] =     ligandConfig.files['modes']
        pairFiles['modesLig_heavy'] = ligandConfig.files['modes_heavy']
        pairFiles['gridLig'] =      ligandConfig.files['grid']
        pairFiles['alphabetLig'] =  ligandConfig.files['alphabetPartner']

        pairFiles['receptorRef'] = receptorRefConfig.files['reduce']
        pairFiles['ligandRef'] = ligandRefConfig.files['reduce']
        pairFiles['receptorRef_heavy'] = receptorRefConfig.files['heavy']
        pairFiles['ligandRef_heavy'] = ligandRefConfig.files['heavy']


        if use_cut:
            receptorConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "A","unbound")
            ligandConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "B","unbound")

            receptorConfig.files['superimpose']['name'] = '{}{}-{}'.format(protein, "A","unbound")
            ligandConfig.files['superimpose']['name'] = '{}{}-{}'.format(protein, "B","unbound")

            receptorConfig.files['cutlog']['name'] = '{}{}-{}'.format(protein, "A","unbound")
            ligandConfig.files['cutlog']['name'] = '{}{}-{}'.format(protein, "B","unbound")

            receptorConfig.files['cut']['name'] = '{}{}-{}'.format(protein, "A","unbound")
            ligandConfig.files['cut']['name'] = '{}{}-{}'.format(protein, "B","unbound")

            receptorRefConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "A","refe")
            ligandRefConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "B","refe")

            receptorRefConfig.files['cut']['name'] = '{}{}-{}'.format(protein, "A","refe")
            ligandRefConfig.files['cut']['name'] = '{}{}-{}'.format(protein, "B","refe")

            receptorConfig.settings['allAtom']['in']['protein']  = 'cut'
            ligandConfig.settings['allAtom']['in']['protein'] = 'cut'

            receptorConfig.settings['secondary']['in']['pdb']  = 'cut'
            ligandConfig.settings['secondary']['in']['pdb'] = 'cut'

            receptorRefConfig.files['cut']['name']  = '{}{}-{}'.format(protein, "A","refe")
            ligandRefConfig.files['cut']['name'] = '{}{}-{}'.format(protein, "B","refe")
            
            # receptorRefConfig.settings['superimpose']['in']['pdb'] =    'cut'
            # ligandRefConfig.settings['superimpose']['in']['pdb'] =      'cut'
            
            receptorRefConfig.settings['allAtom']['in']['protein'] =    'cut'
            ligandRefConfig.settings['allAtom']['in']['protein'] =      'cut'
            
            #receptorRefConfig.files['refpdb'] = receptorConfig.files['cut'] 
            #ligandRefConfig.files['refpdb'] =   ligandConfig.files['cut'] 

            # receptorConfig.files['refpdb'] = receptorRefConfig.files['cut'] 
            # ligandConfig.files['refpdb'] =   ligandRefConfig.files['cut'] 

            receptorRefConfig.files['cutlog'] = receptorConfig.files['cutlog'] 
            ligandRefConfig.files['cutlog'] = ligandConfig.files['cutlog'] 

        if use_pruning:
            receptorConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "A","unbound")
            ligandConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "B","unbound")
            receptorRefConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "A","refe")
            ligandRefConfig.files['pdb']['name'] = '{}{}-{}'.format(protein, "B","refe")            
            
            pairFiles['receptor'] =             receptorConfig.files['prune']
            pairFiles['ligand'] =             ligandConfig.files['prune']

            pairFiles['receptorRef'] =             receptorRefConfig.files['prune']
            pairFiles['ligandRef'] =             ligandRefConfig.files['prune']

            receptorConfig.files['grid']['extension']  = "-grid-pruned.grid"
            ligandConfig.files['grid']['extension']  = "-grid-pruned.grid"
            receptorConfig.settings['grid']['in']['protein']  = 'prune'
            ligandConfig.settings['grid']['in']['protein']  = 'prune'

            # receptorConfig.settings['modes']['in']['protein']  = 'prune'
            # ligandConfig.settings['modes']['in']['protein']  = 'prune'

            pairFiles['gridRec'] =      receptorConfig.files['grid']
            pairFiles['gridLig'] =      ligandConfig.files['grid']
            pairFiles['dof']["extension"] =     "-dof-pruned.dat"
            pairFiles['dof_test']["extension"] =     "-dof_test-pruned.dat"

        if use_boundModes:
            # ligandConfig.files['modes'] = ligandConfig.files['bound_modes']
            # receptorConfig.files['modes'] = receptorConfig.files['bound_modes']
            
            # ligandConfig.files['modes_heavy'] = ligandConfig.files['bound_modes_heavy']
            # receptorConfig.files['modes_heavy'] = receptorConfig.files['bound_modes_heavy']

            receptorConfig.settings['mode_evaluation']['in']['mode_file'] =  'bound_modes'
            ligandConfig.settings['mode_evaluation']['in']['mode_file'] =  'bound_modes'

            pairFiles['modesRec_heavy'] =     receptorConfig.files['bound_modes_heavy']
            pairFiles['modesLig_heavy'] =     ligandConfig.files['bound_modes_heavy']
            pairFiles['modesRec'] =     receptorConfig.files['bound_modes']
            pairFiles['modesLig'] =     ligandConfig.files['bound_modes']

        

        if use_singleDof:
            pairConfig.settings['docking']['in']['dof'] = 'dof_test'

        if use_manipulated:
            # receptorConfig.settings['modes']['dryRun'] =  True
            # ligandConfig.settings['modes']['dryRun'] =  True

            # receptorConfig.settings['modes_heavy']['dryRun'] =  True
            # ligandConfig.settings['modes_heavy']['dryRun'] =  True


            # ligandConfig.files['modes'] = ligandConfig.files['modes_manipulate']
            # receptorConfig.files['modes'] = receptorConfig.files['modes_manipulate']
            
            # ligandConfig.files['modes_heavy'] = ligandConfig.files['modes_manipulate_heavy']
            # receptorConfig.files['modes_heavy'] = receptorConfig.files['modes_manipulate_heavy']

            # receptorConfig.settings['mode_evaluation']['in']['mode_file'] =  'modes_manipulate'
            # ligandConfig.settings['mode_evaluation']['in']['mode_file'] =  'modes_manipulate'

            pairFiles['modesRec'] =     receptorConfig.files['modes_manipulate']
            pairFiles['modesLig'] =     ligandConfig.files['modes_manipulate']        
            pairFiles['modesRec_heavy'] =     receptorConfig.files['modes_manipulate_heavy']
            pairFiles['modesLig_heavy'] =     ligandConfig.files['modes_manipulate_heavy']    
            pairFiles['modesLig_heavy'] =     ligandConfig.files['modes_manipulate_heavy']    
            pairFiles['joinedModes']['extension'] = "-joinedModes-r{}-l{}-{}.dat".format(20,20,'manipulated')
            pairFiles['joinedModes_heavy']['extension'] = "-joinedModes-heavy-r{}-l{}-{}.dat".format(20,20,'manipulated')

        # pairConfig.settings['interface_low']['overwrite'] = True
        # pairConfig.settings['interface_high']['overwrite'] = True
        # receptorConfig.settings['mode_evaluation']['overwrite'] = True
        # ligandConfig.settings['mode_evaluation']['overwrite'] = True
        # receptorConfig.settings['protein_evaluation']['overwrite'] = True
        # ligandConfig.settings['protein_evaluation']['overwrite'] = True


        if config['rigidStart']:
            pairFiles['rigidStart'] = {}
            pairFiles['rigidStart']['folder'] = "bm_dG_mr0_ml0_s1p000000_sO_c50_mr0_ml0_s1p000000_hin99/result"
            pairFiles['rigidStart']['name'] = "{}".format(protein)
            pairFiles['rigidStart']['extension'] = "-unbound-docking.dat"
            pairConfig.settings['docking']['in']['dof'] = 'rigidStart'

        if config['modesOnly']:
            pairConfig.settings['docking']['modesOnly'] = True

        #pairConfig.settings['dof_evaluation_low']['overwrite'] = True
        #pairConfig.settings['dof_evaluation_high']['overwrite'] = True
        #pairConfig.settings['dof_evaluation']['overwrite'] = True

        pipelineRec.put((copy.deepcopy(     receptorConfig),protein))         
        pipelineLig.put((copy.deepcopy(     ligandConfig),protein))                  
        pipelineRec2.put((copy.deepcopy(    receptorConfig),protein))                 
        pipelineLig2.put((copy.deepcopy(    ligandConfig),protein))                  
        pipelineRecRef.put((copy.deepcopy(  receptorRefConfig),protein))      
        pipelineLigRef.put((copy.deepcopy(  ligandRefConfig),protein))      
        pipelinePairConfig.put((copy.deepcopy(pairConfig),protein))  
        pipelinePairRun.put((copy.deepcopy( pairConfig),protein))     


print("\n DO FIRST CONFIGURATION")
logging.warning("\n\n---------------------DO FIRST CONFIGURATION-------------------------\n")

pipelineRec.start()
pipelineLig.start()
pipelineRec.join()
pipelineLig.join()


print("\n DO REFERENCE CONFIGURATION")
logging.warning("\n\n---------------------DO REFERENCE CONFIGURATION---------------------\n")

pipelineRecRef.start()
pipelineLigRef.start()
pipelineRecRef.join()
pipelineLigRef.join()

print("\nDO SECOND CONFIGURATION")
logging.warning("\n\n---------------------DO SECOND CONFIGURATION------------------------\n")

pipelineRec2.start()
pipelineLig2.start()
pipelineRec2.join()
pipelineLig2.join()

print("\nDO PAIR CONFIGURATION")
logging.warning("\n\n---------------------DO PAIR CONFIGURATION--------------------------\n")

pipelinePairConfig.start()
pipelinePairConfig.join()

# print("\nDO RUN")
# logging.warning("\n\n---------------------DO RUN-----------------------------------------\n")

# pipelinePairRun.start()
# pipelinePairRun.join()






  













# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.5, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new'},
# {"numModesRec":1,  "numModesLig":1, "scale": 2, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new'},






# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut_fresh'},

# {"numModesRec":0,  "numModesLig":0, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':True, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_bound_fresh'},

# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_fresh'},

# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':True, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_manipulated_fresh'},

# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': True, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_prune_fresh'},

# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_fresh'},








# {"numModesRec":0,  "numModesLig":1, "scale": 1, 'prune': True, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned_2'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.5, 'prune': True, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned'},

# {"numModesRec":1,  "numModesLig":1, "scale": 1.5, 'prune': True, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned'},

# {"numModesRec":1,  "numModesLig":1, "scale": 2.0, 'prune': True, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned'},

# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': True, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned'},

# {"numModesRec":1,  "numModesLig":1, "scale": 10, 'prune': True, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned'},


# #cutPrunded

# {"numModesRec":1,  "numModesLig":1, "scale": 0.1, 'prune': True, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut_pruned'},


# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': True, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut_pruned_2'},


# {"numModesRec":1,  "numModesLig":1, "scale": 2.5, 'prune': True, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut_pruned'},



# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': True, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut_pruned'},


# #cut

# {"numModesRec":1,  "numModesLig":1, "scale": 0.5, 'prune': False, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut'},

# {"numModesRec":1,  "numModesLig":1, "scale": 1.5, 'prune': False, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut'},

# {"numModesRec":1,  "numModesLig":1, "scale": 2.0, 'prune': False, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut'},

# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': False, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut'},

# {"numModesRec":1,  "numModesLig":1, "scale": 10, 'prune': False, 'cut':True,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_cut'},


# #boundModes

# {"numModesRec":1,  "numModesLig":1, "scale": 0.01, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.8, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes'},

# {"numModesRec":1,  "numModesLig":1, "scale": 2.0, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes'},

# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':False, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes'},


# #doftest

# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_modesOnly_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.5, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_modesOnly_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 2, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_modesOnly_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_modesOnly_singleDof'},



# #doftest 

# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': True, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned_modesOnly_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.5, 'prune': True, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned_modesOnly_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 2, 'prune': True, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned_modesOnly_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': True, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_pruned_modesOnly_singleDof'},











# {"numModesRec":1,  "numModesLig":1, "scale": 0.1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 2.5, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 10, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_new_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 0.5, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 2.5, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 10, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},

# {"numModesRec":1,  "numModesLig":1, "scale": 0.5, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_modesOnly_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 0.1, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_modesOnly_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 2.5, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_modesOnly_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 5, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_modesOnly_singleDof'},
# {"numModesRec":1,  "numModesLig":1, "scale": 10, 'prune': False, 'cut':False,'modesOnly':True,
# 'bound':False, 'boundModes': True,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_boundModes_modesOnly_singleDof'},


# {"numModesRec":1,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},
# {"numModesRec":1,  "numModesLig":0, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},
# {"numModesRec":0,  "numModesLig":1, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},
# {"numModesRec":0,  "numModesLig":0, "scale": 1, 'prune': False, 'cut':False,'modesOnly':False,
# 'bound':False, 'boundModes': False,'singleDof':True, 'manipulateModes':False, 
# 'rigidStart': False, 'overwrite': False, 'extension':'_hin99_singleDof'},