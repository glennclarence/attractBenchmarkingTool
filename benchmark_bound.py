from Configurators import *
from Settings import *
from Configuration import Configuration
from Pipeline import *

from ConfiguratorBase import createLoggingFile
import math
from datetime import datetime
import sys



#create Modes, coarse grain Proteins and alphabet files
singleConfiguration_01 = [
    {
    'conf': [ 
        #{'setting': 'findtermini',  'configurator': 'findtermini'},
        #{'setting': 'cut',          'configurator': 'cut'},
        {'setting': 'allAtom',      'configurator': 'allAtom'},
        {'setting': 'reduce',       'configurator': 'reduce'},
        {'setting': 'heavy',        'configurator': 'heavy'},
        {'setting': 'modes',        'configurator': 'modes'},
        {'setting': 'modes_heavy',  'configurator': 'modes'},
        {'setting': 'alphabet',     'configurator': 'alphabet'},
        {'setting': 'secondary',    'configurator': 'secondary'},
        {'setting': 'mode_manipulate','configurator': 'mode_manipulate'},
        {'setting': 'saveSettings', 'configurator': 'saveSettings'}
        ],
    'numThreads': 1}
]
#for reference structure only apply coarse grain conversion
singleConfigurationRef = [
    {'conf':[
        #{ 'setting': 'cut' ,        'configurator': 'cut'},
        { 'setting': 'superimpose' ,'configurator': 'superimpose'},
        { 'setting': 'allAtom' ,    'configurator': 'allAtom'},
        { 'setting': 'reduce' ,     'configurator': 'reduce'},
        { 'setting': 'heavy' ,     'configurator':  'heavy'},
        { 'setting': 'saveSettings','configurator': 'saveSettings'} 
        ] ,
    'numThreads': 2}
]

#create Grids. this has to be a seperate step since the alphabet file of the partner is needed
singleConfiguration_02 = [
    {'conf':
        { 'setting': 'grid',        'configurator':'grid'},
    'numThreads': 1}
]

modeEvaluation = [
    {'conf':
        [
        { 'setting': 'mode_evaluation',        'configurator':'mode_evaluation'},
       # { 'setting': 'bound_mode',        'configurator':'bound_mode'}
        ],
    'numThreads': 1}
]

#All files are created which needed receptor and ligand
pairConfiguration = [
    {'conf':
        [
       #{ 'setting' : 'dof_test',    'configurator': 'dof_test'},
        { 'setting': 'dof',         'configurator': 'dof'},
        { 'setting': 'joinModes',   'configurator': 'joinModes'} ,
        { 'setting': 'joinModes_heavy',   'configurator': 'joinModes'} 
        ], 
    'numThreads':1
    }
]

#perform docking and scoring as well analysis
runConfiguration = [
    {'conf':
        { 'setting':'docking',     'configurator':'docking'},
    'numThreads': 1},
    {'conf':
         { 'setting':'scoring','configurator':'scoring'},
    'numThreads': 2},
    {'conf':[
        { 'setting': 'fill_energy',     'configurator': 'fill_energy'  },
        { 'setting': 'sorting',         'configurator': 'sorting'      },
        { 'setting': 'deredundant',     'configurator': 'deredundant'  },
        { 'setting': 'top',             'configurator': 'top'          },
        { 'setting': 'demode',          'configurator': 'demode'       },
        { 'setting': 'irmsd',           'configurator': 'irmsd'        },
        ##{ 'setting': 'irmsd_nomodes',   'configurator': 'irmsd'        },
        { 'setting': 'rmsd',            'configurator': 'rmsd'          },
        ##{ 'setting': 'rmsd_nomodes',    'configurator': 'rmsd'         },
        { 'setting': 'fnat',            'configurator': 'fnat'          },
        ##{ 'setting': 'fnat_nomodes',    'configurator': 'fnat'         },
        { 'setting': 'saveSettings',    'configurator': 'saveSettings'  },
        { 'setting': 'collect'  ,       'configurator': 'collect'    },
        { 'setting': 'dof_evaluation',  'configurator': 'dof_evaluation'},
        #{ 'setting': 'interface',       'configurator':'interface'}
                ],
    'numThreads': 1}
    ]


deleted_proteins = ['1BJ1', '1FSK', '1IQD', '1K4C', '1KXQ', '1NCA', '1NSN', '1QFW', '2HMI', '2JEL', '9QFW','1DE4','4FQI','4GAM','4GXU','1EXB','4FQI', '1EER','4H30' ]

protType = "refe"
protTypeRef = "refe"
pathfile = "paths_local.json"



paths = json.load(open(pathfile,'r'))
basePath = paths["basePath"]
utils.strideBinary = paths["strideBinary"]
utils.rmscaBinary = paths["rmscaBinary"]

protlist_file = basePath + "/protlist"

proteins =[]
with open(protlist_file ,'r') as f:
    lines = f.readlines()
    for line in lines:
        proteins.append(line.split()[0])
try:
    for dele in deleted_proteins:
        proteins.remove(dele)
except:
    pass
#proteins =['7CEI']



configs = [
{"numModesRec":0,  "numModesLig":0, "scale": 1},
#{"numModesRec":1,  "numModesLig":0, "scale": 1},
#{"numModesRec":0,  "numModesLig":1,"scale": 1},
#{"numModesRec":1, "numModesLig":1,"scale": 1}
]


extension = "_hin99_bound"
prefix = "00_"
for i,c in enumerate (configs):
    prefix += "_BM{}_mr{}ml{}s{}".format(i, c['numModesRec'], c['numModesLig'], c['scale'])
prefix += extension

createLoggingFile(basePath+"/{}_loggingFile_{}.log".format(prefix,str(datetime.now().date())))
consoleOutputFil = basePath+"/{}_ConsoleOutput_{}.log".format(prefix,str(datetime.now().date()))
sys.stdout = open(consoleOutputFil, 'w')

dry = False
verbose = True
overwrite = False
num = len(proteins)

for i,config in enumerate(configs):
    numModesRec = config['numModesRec']
    numModesLig = config['numModesLig']
    scale = config['scale']
    bufferSize = 1000

    frac, whole = math.modf(scale)
    bm = "bm_dG_mr{}_ml{}_s{}p{}_sO_c50_mr{}_ml{}_s{}p{}{}".format(numModesRec, numModesLig, int(whole),'{:6f}'.format(frac)[2:],numModesRec, numModesLig, int(whole),'{:6f}'.format(frac)[2:],extension)
    print(" START NEW BENCHMARK NUMBER {} mr{} ml{} scale{} name {}".format(i,numModesRec,numModesLig, scale,bm))
        
    print("\nset up buffer queues")
    logging.warning("\nset up buffer queues")
    inputRecQ = queue.Queue(bufferSize)
    inputLigQ = queue.Queue(bufferSize)
    inputRecQ_2 = queue.Queue(bufferSize)
    inputLigQ_2 = queue.Queue(bufferSize)
    inputRecRefQ = queue.Queue(bufferSize)
    inputLigRefQ = queue.Queue(bufferSize)
    outputRecRefQ = queue.Queue(bufferSize)
    outputLigRefQ = queue.Queue(bufferSize)

    outputRecQ = queue.Queue(bufferSize)
    outputLigQ = queue.Queue(bufferSize)

    outputRecQ2 = queue.Queue(bufferSize)
    outputLigQ2 = queue.Queue(bufferSize)


    inputRecModeEvalQ = queue.Queue(bufferSize)
    inputLigModeEvalQ = queue.Queue(bufferSize)
    outputRecModeEvalQ = queue.Queue(bufferSize)
    outputLigModeEvalQ = queue.Queue(bufferSize)

    pairQ = queue.Queue(bufferSize)
    pairQOutConfig = queue.Queue(bufferSize)
    pairQOutRes = queue.Queue(bufferSize)


    print("create configrations".upper())
    logging.warning("\n\n---------------------CREATE CONFIGURATIONS--------------------------\n".upper())

    for protein in proteins:
        numModes = 20
        receptorConfig  = Configuration(getDefaultSingleSetting(    protein=protein, chain="A", protType = protType, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite, 
            pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"]))
        ligandConfig    = Configuration(getDefaultSingleSetting(    protein=protein, chain="B", protType = protType, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite,
        pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"]))
        receptorRefConfig = Configuration(getDefaultSingleSetting(  protein=protein, chain="A", protType = protTypeRef, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite,
            pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"]))
        ligandRefConfig = Configuration(getDefaultSingleSetting(    protein=protein, chain="B", protType = protTypeRef, numModes = numModes,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose, overwrite = overwrite,
            pythonBinary = paths['python2Binary'],
            attractToolPath = paths["attractToolPath"],
            attractBinPath = paths["attractBinPath"],
            attractParFile = paths["attractParFile"]))


        receptorRefConfig.settings['allAtom']['in']['protein'] =    'superimpose'
        ligandRefConfig.settings['allAtom']['in']['protein'] =      'superimpose'

        receptorRefConfig.files['refpdb']   = receptorConfig.files['pdb'] 
        ligandRefConfig.files['refpdb']     =   ligandConfig.files['pdb'] 

        receptorRefConfig.files['cutlog']   = receptorConfig.files['cutlog'] 
        ligandRefConfig.files['cutlog']     = ligandConfig.files['cutlog'] 

        receptorConfig.files['alphabetPartner'] = ligandConfig.files['alphabet'] 
        ligandConfig.files['alphabetPartner'] = receptorConfig.files['alphabet'] 

        pairConfig      = Configuration(getDefaultPairSetting(benchmarkName = bm,protein = protein, protType = protType, protTypeRef = protTypeRef, 
        numModesRec = numModesRec,numModesLig = numModesLig,basePath = basePath + "/{}".format(protein), dry=dry, verbose = verbose, overwrite = overwrite,
        attractBinary   =  paths["attractBinary"] , 
        attractBinPath  =  paths["attractBinPath"], 
        attractParFile  =  paths["attractParFile"],
        dofBinary       =  paths["dofBinary"],
        pythonBinary    =  paths["python2Binary"],
        attractToolPath =  paths["attractToolPath"],
        attractBinaryGPU = paths["attractBinaryGPU"],
        deviceIds = [0,1], evScale=scale))

        pairFiles = pairConfig.files

        pairFiles['mode_evaluation_rec'] = receptorConfig.files['mode_evaluation']
        pairFiles['mode_evaluation_lig'] = ligandConfig.files['mode_evaluation']

        pairFiles['receptor'] =             receptorConfig.files['reduce']
        pairFiles['receptor_heavy'] =       receptorConfig.files['heavy']

        pairFiles['modesRec'] =     receptorConfig.files['modes']
        pairFiles['modesRec_heavy'] =     receptorConfig.files['modes_heavy']

        pairFiles['gridRec'] =      receptorConfig.files['grid']
        pairFiles['alphabetRec'] =  receptorConfig.files['alphabetPartner']

        pairFiles['ligand'] =       ligandConfig.files['reduce']
        pairFiles['ligand_heavy'] = ligandConfig.files['heavy']
        pairFiles['modesLig'] =     ligandConfig.files['modes']
        pairFiles['modesLig_heavy'] =     ligandConfig.files['modes_heavy']
        pairFiles['gridLig'] =      ligandConfig.files['grid']
        pairFiles['alphabetLig'] =  ligandConfig.files['alphabetPartner']

        pairFiles['receptorRef'] =  receptorRefConfig.files['reduce']
        pairFiles['ligandRef'] =    ligandRefConfig.files['reduce']
        pairFiles['receptorRef_heavy'] =  receptorRefConfig.files['heavy']
        pairFiles['ligandRef_heavy'] =    ligandRefConfig.files['heavy']
        

        ligand_modeEval = Configuration(getDefaultSingleSetting(    protein=protein, chain="B", protType = protType, numModes = numModesLig,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose))
        receptor_modeEval = Configuration(getDefaultSingleSetting(    protein=protein, chain="A", protType = protType, numModes = numModesRec,basePath = basePath + "/{}".format(protein), dry=dry ,verbose = verbose))
        ligand_modeEval.files['partner_bound'] = ligandRefConfig.files['reduce']
        receptor_modeEval.files['partner_bound'] = receptorRefConfig.files['reduce']



        inputRecModeEvalQ.put((receptor_modeEval,protein))
        inputLigModeEvalQ.put((ligand_modeEval,protein))

        inputRecQ.put((copy.deepcopy(receptorConfig),protein))
        inputRecQ_2.put((copy.deepcopy(receptorConfig),protein))
        inputLigQ.put((copy.deepcopy(ligandConfig),protein))
        inputLigQ_2.put((copy.deepcopy(ligandConfig),protein))

        inputRecRefQ.put((receptorRefConfig,protein))
        inputLigRefQ.put((ligandRefConfig,protein))
        pairQ.put((pairConfig,protein))


    #configure receptor and ligand in two steps
    pipelineRec         = createPipeline( inputRecQ, outputRecQ,bufferSize, singleConfiguration_01,  num)
    pipelineLig         = createPipeline( inputLigQ, outputLigQ,bufferSize, singleConfiguration_01, num)
    #configure grids and joined Modes
    pipelineRec2        = createPipeline( inputRecQ_2, outputRecQ2,bufferSize, singleConfiguration_02, inputRecQ_2.qsize())
    pipelineLig2        = createPipeline( inputLigQ_2, outputLigQ2,bufferSize, singleConfiguration_02, inputLigQ_2.qsize())
    #configure referenceStructures
    pipelineRecRef      = createPipeline( inputRecRefQ, outputRecRefQ,bufferSize, singleConfigurationRef,  num)
    pipelineLigRef      = createPipeline( inputLigRefQ, outputLigRefQ,bufferSize, singleConfigurationRef,  num)
    #runpairconfiguration
    pipelinePairConfig  = createPipeline( pairQ, pairQOutConfig,bufferSize, pairConfiguration,  num)
    #run docking, scoring and analysis
    pipelinePairRun     = createPipeline( pairQOutConfig, pairQOutRes,bufferSize, runConfiguration,  num)

    pipelineModeEvalLig= createPipeline( inputLigModeEvalQ, outputLigModeEvalQ,bufferSize, modeEvaluation,  num)
    pipelineModeEvalRec= createPipeline( inputRecModeEvalQ, outputRecModeEvalQ,bufferSize, modeEvaluation,  num)

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

    print("\nDO RUN")
    logging.warning("\n\n---------------------DO RUN-----------------------------------------\n")

    pipelinePairRun.start()
    pipelinePairRun.join()

    # logging.warning("\n\n---------------------MODE EVALUATION--------------------------------\n")

    pipelineModeEvalRec.start()
    pipelineModeEvalRec.join()

    pipelineModeEvalLig.start()
    pipelineModeEvalLig.join()

