from Configurators import *
from Settings import *
from Configuration import Configuration
from Pipeline import *

#create Modes, coarse grain Proteins and alphabet files
singleConfiguration_01 = [
    {'conf':
        [
        ReduceConfigurator,
        ModeConfigurator,
        AlphabetConfigurator
        ],
    'numThreads': 1}
]

#for reference structure only apply coarse grain conversion
singleConfigurationRef = [
    {'conf':
        ReduceConfigurator  ,
    'numThreads': 1}
]

#create Grids. this has to be a seperate step since the alphabet file of the partner is needed
singleConfiguration_02 = [
    {'conf':
        GridConfigurator, 
    'numThreads': 1}
]

#All files are created which needed receptor and ligand
pairConfiguration = [
    {'conf':
        [DofConfigurator,
        JoinModesConfigurator], 
    'numThreads':1
    }

]

#perform docking and scoring as well analysis
runConfiguration = [
    {'conf':
        DockingConfigurator,
    'numThreads': 2},
    {'conf':
        ScoringConfigurator,
    'numThreads': 1},
    {'conf':[
        SortingConfigurator,
        DeRedundantConfigurator,
        DemodeConfigurator,
        TopConfigurator,
        IRMSDConfigurator,
        RMSDConfigurator,
        FNATConfigurator,
        CollectConfigurator],
    'numThreads': 1}
    ]




protein="3MXW"
chain="A"
protType = "unbound"
protTypeRef = "refe"
basePath = "basePath"
numModesRec = 5
numModesLig = 5
bm = "test"
dry = True
verbose = True

pairConfig      = Configuration(getDefaultPairSetting(benchmarkName = bm,protein = protein, protType = protType, protTypeRef = protTypeRef, numModesRec = numModesRec,numModesLig = numModesLig,basePath = basePath, dry=dry, verbose = verbose ))


receptorConfig  = Configuration(getDefaultSingleSetting(    protein=protein, chain="A", protType = protType, numModes = numModesRec,basePath = basePath, dry=dry ,verbose = verbose))
ligandConfig    = Configuration(getDefaultSingleSetting(    protein=protein, chain="B", protType = protType, numModes = numModesLig,basePath = basePath, dry=dry ,verbose = verbose))
receptorRefConfig = Configuration(getDefaultSingleSetting(  protein=protein, chain="A", protType = protTypeRef, numModes = numModesRec,basePath = basePath, dry=dry ,verbose = verbose))
ligandRefConfig = Configuration(getDefaultSingleSetting(    protein=protein, chain="B", protType = protTypeRef, numModes = numModesLig,basePath = basePath, dry=dry ,verbose = verbose))

receptorConfig.files['alphabetPartner'] = ligandConfig.files['alphabet'] 
ligandConfig.files['alphabetPartner'] = receptorConfig.files['alphabet'] 

pairFiles = pairConfig.files
pairFiles['receptor'] =     receptorConfig.files['reduce']
pairFiles['modesRec'] =     receptorConfig.files['modes']
pairFiles['gridRec'] =      receptorConfig.files['grid']
pairFiles['alphabetRec'] =  receptorConfig.files['alphabetPartner']


pairFiles['ligand'] =       ligandConfig.files['reduce']
pairFiles['modesLig'] =     ligandConfig.files['modes']
pairFiles['gridLig'] =      ligandConfig.files['grid']
pairFiles['alphabetLig'] =  ligandConfig.files['alphabetPartner']

pairFiles['receptorRef'] =  receptorRefConfig.files['reduce']
pairFiles['ligandRef'] =    ligandRefConfig.files['reduce']




bufferSize = 10
inputRecQ = queue.Queue(bufferSize)
inputLigQ = queue.Queue(bufferSize)

inputRecRefQ = queue.Queue(bufferSize)
inputLigRefQ = queue.Queue(bufferSize)
outputRecRefQ = queue.Queue(bufferSize)
outputLigRefQ = queue.Queue(bufferSize)

outputRecQ = queue.Queue(bufferSize)
outputLigQ = queue.Queue(bufferSize)

outputRecQ2 = queue.Queue(bufferSize)
outputLigQ2 = queue.Queue(bufferSize)
pairQ = queue.Queue(bufferSize)
pairQOutConfig = queue.Queue(bufferSize)
pairQOutRes = queue.Queue(bufferSize)

num = 10
for i in range(num):
    inputRecQ.put((receptorConfig,i))
    inputLigQ.put((ligandConfig,i))
    inputRecRefQ.put((receptorRefConfig,i))
    inputLigRefQ.put((ligandRefConfig,i))
    pairQ.put((pairConfig,i))

#configure receptor and ligand in two steps
pipelineRec         = createPipeline( inputRecQ, outputRecQ,bufferSize, singleConfiguration_01,  num)
pipelineLig         = createPipeline( inputLigQ, outputLigQ,bufferSize, singleConfiguration_01, num)

#configure grids and joined Modes
pipelineRec2        = createPipeline( outputRecQ, outputRecQ2,bufferSize, singleConfiguration_02,  outputRecQ.qsize())
pipelineLig2        = createPipeline( outputLigQ, outputLigQ2,bufferSize, singleConfiguration_02, outputLigQ.qsize())

#configure referenceStructures
pipelineRecRef      = createPipeline( inputRecRefQ, outputRecRefQ,bufferSize, singleConfigurationRef,  num)
pipelineLigRef      = createPipeline( inputLigRefQ, outputLigRefQ,bufferSize, singleConfigurationRef,  num)

#runpairconfiguration
pipelinePairConfig  = createPipeline( pairQ, pairQOutConfig,bufferSize, pairConfiguration,  num)
#run docking, scoring and analysis
pipelinePairRun     = createPipeline( pairQOutConfig, pairQOutRes,bufferSize, runConfiguration,  num)

print("\n DO REFERENCE CONFIGURATION")
pipelineRecRef.start()
pipelineLigRef.start()
pipelineRecRef.join()
pipelineLigRef.join()

print("\n DO FIRST CONFIGURATION")
pipelineRec.start()
pipelineLig.start()
pipelineRec.join()
pipelineLig.join()

print("\nDO SECOND CONFIGURATION")
pipelineRec2.start()
pipelineLig2.start()
pipelineRec2.join()
pipelineLig2.join()

print("\n DO PAIR CONFIGURATION")
pipelinePairConfig.start()
pipelinePairConfig.join()

print("\n DO RUN")
pipelinePairRun.start()
pipelinePairRun.join()


# print("\nRUN RECEPTOR CONFIG")
# for configurator in singleConfiguration_01:
#     configurator.setConfig(receptorConfig)
#     configurator.run()


# print("\nRUN RECEPTOR REF CONFIG")
# for configurator in singleConfigurationRef:
#     configurator.setConfig(receptorRefConfig)
#     configurator.run()


# for configurator in singleConfiguration_02:
#     configurator.setConfig(receptorConfig)
#     configurator.run()

# print("\nRUN LIGAND CONFIG")
# for configurator in singleConfiguration_01:
#     configurator.setConfig(ligandConfig)
#     configurator.run()

# for configurator in singleConfiguration_02:
#     configurator.setConfig(ligandConfig)
#     configurator.run()

# print("\nRUN RECEPTOR LIG CONFIG")
# for configurator in singleConfigurationRef:
#     configurator.setConfig(ligandRefConfig)
#     configurator.run()

# print("\nRUN PAIR CONFIG")
# for configurator in pairConfiguration:
#     configurator.setConfig(pairConfig)
#     configurator.run()

# for configurator in runConfiguration:
#     configurator.setConfig(pairConfig)
#     configurator.run()
