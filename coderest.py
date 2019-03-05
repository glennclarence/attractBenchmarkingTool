

# class ReduceConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = reduceFunction
#         self.setting = 'reduce'
#         super(ReduceConfigurator,self).__init__(self.setting)
# class ModeConfigurator(ConfiguratorBase):
#     def __init__(self,):
#         self.configFunc = modeFunction
#         self.setting = 'modes'
#         super(ModeConfigurator,self).__init__( self.setting)
# class AlphabetConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = alphabetFunction
#         self.setting = 'alphabet'
#         super(AlphabetConfigurator,self).__init__(self.setting)
 # class GridConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = gridFunction
#         self.setting = 'grid'
#         super(GridConfigurator,self).__init__( self.setting)
# class DofConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = dofFunction
#         self.setting = 'dof'
#         super(DofConfigurator,self).__init__( self.setting)
 # class JoinModesConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = joinModesFunction
#         self.setting = 'joinModes'
#         super(JoinModesConfigurator,self).__init__( self.setting)
# class CAConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = CAFunction
#         self.setting = 'CA'
#         super(CAConfigurator,self).__init__( self.setting)
# class DockingConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = dockingFunction
#         self.setting = 'docking'
#         super(DockingConfigurator,self).__init__( self.setting)
# class ScoringConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = scoringFunction
#         self.setting = 'scoring'
#         super(ScoringConfigurator,self).__init__( self.setting)
 

# class SortingConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = SortingFunction
#         self.setting = 'sorted'
#         super(SortingConfigurator,self).__init__( self.setting)



# class RedundantConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = RedundantFunction
#         self.setting = 'redundant'
#         super(RedundantConfigurator,self).__init__( self.setting)

# class DemodeConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = DemodeFunction
#         self.setting = 'demode'
#         super(DemodeConfigurator,self).__init__( self.setting)

# class TopConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = TopFunction
#         self.setting = 'top'
#         super(TopConfigurator,self).__init__( self.setting)

# class IRMSDConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = IRMSDFunction
#         self.setting = 'irmsd'
#         super(IRMSDConfigurator,self).__init__( self.setting)

# class RMSDConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = RMSDFunction
#         self.setting = 'rmsd'
#         super(RMSDConfigurator,self).__init__( self.setting)

# class FNATConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = FNATFunction
#         self.setting = 'fnat'
#         super(FNATConfigurator,self).__init__( self.setting)

# class CollectConfigurator(ConfiguratorBase):
#     def __init__(self):
#         self.configFunc = CollectFunction
#         self.setting = 'collect'
#         super(CollectConfigurator,self).__init__( self.setting)






# modeConfig= {
#     "id":"{}{}-m".format(protein, chain ),
#     "files":
#         {

#             {
#                 'pdb': '{}{}-refe.pdb'.format(protein, chain)
#             }, 
#         'output':
#             {'modefile':'{}{}-modes-{}.dat'.format(protein, chain,numModes)}
#         },
#     "folders":
#         {
#         'input':'', 
#         'output':'modes'
#         },
#     "settings": {
#         "attractToolPath": "/home/glenn/Documents/attract/tools",
#         "numModes":numModes,
#         "dryRun": dry,
#         "verbose": verbose
#     },
#     "basePath":basePath,
#     "overwrite": overwrite,
#     'extensions':{
#         'modefile':'modes.dat'
#     }
# }

# reduceConfig= {
#     "id":"{}{}-r".format(protein, chain ),
#     "files":
#         {
#         'input':
#             {
#                 'pdb': '{}{}-refe.pdb'.format(protein, chain)
#             }, 
#         'output':
#             {'reduce':'{}{}-r.dat'.format(protein, chain)}
#         },
#     "folders":{'input':'', 'output':'reduced'},
#     "settings": {
#         "attractToolPath": "/home/glenn/Documents/attract/tools",
#         "dryRun": dry,
#         "verbose": verbose
#     },
#     "basePath":basePath,
#     "overwrite": overwrite,
#     'extensions':{
#         'modefile':'modes.dat'
#     }
# }

# alphConfig= {
#     "id":"{}{}-alph".format(protein, chain ),
#     "files":
#         {
#         'input':
#             {
#                 'pdb': '{}{}-refe.pdb'.format(protein, 'A' if chain == 'B' else 'B')
#             }, 
#         'output':
#             {'alphabet':'{}{}-alphabet.grid'.format(protein, chain)}
#         },
#     "folders":{'input':'', 'output':'grid'},
#     "settings": {
#         "dryRun": dry,
#         "verbose": verbose
#     },
#     "basePath":basePath,
#     "overwrite": overwrite,
#     'extensions':{
#         'modefile':'alphabet.grid'
#     }
# }


# gridConfig= {
#     "id":"{}{}-alph".format(protein, chain ),
#     "files":
#         {
#         'input':
#             {
#                 'pdb': '{}{}-r.pdb'.format(protein, chain)
#             }, 
#         'output':
#             {'grid':'{}{}-grid.grid'.format(protein, chain)}
#         },
#     "folders":{'input':'reduced', 'output':'grid'},
#     "settings": {
#         "dryRun": dry,
#         "verbose": verbose
#     },
#     "basePath":basePath,
#     "overwrite": overwrite,
#     'extensions':{
#         'modefile':'alphabet.grid'
#     }
# }



dofConfig ={
    "id":"{}-dof".format(protein),
    "files":
        {
        'receptor':
            {
                'folder':       "",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    ".pdb"
            },
        'ligand':
            {
                'folder':       "",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    ".pdb"
            },
        'translate':
            {
                'folder':       "dof",
                "name":         'translate',
                "extension":    ".dat"
            },
        'dof':
            {
                'folder':       "dof",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-dof.dat"
            }
        },
    "settings": {
        "attractToolPath": attractToolPath,
        "attractBinPath": attractBinPath,
        "attractParFile": attractParFile,
        "basePath":basePath,

        "dof":{
            "in": {"receptor": "receptor","ligand": "ligand"  },
            "out": {"out": "dof", "translate":"translate" },
            "numModes":numModes,
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose,
            "dofbinary": "/home/glenn/Documents/Masterarbeit/git/Attract_benchmark/tools/systseach"
        }
    }
}


joinModesConfig = {
    "id":"{}-joinedModes".format(protein ),
    "files":
        {
        'receptorModes':
            {
                'folder':       "modes",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,modeType)
            },
        'ligandModes':
            {
                'folder':       "modes",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,modeType)
            },
        'joinedModes':
            {
                'folder':       "modes",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-joinedModes-{}-{}.dat".format(numModes,modeType)
            }
        },
    "settings": {
        "attractToolPath": attractToolPath,
        "attractBinPath": attractBinPath,
        "attractParFile": attractParFile,
        "basePath":basePath,

        "joinModes":{
            "in": {"receptorModes": "receptorModes","ligandModes": "ligandModes"  },
            "out": {"out": "joinedModes" },
            "numModes":numModes,
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        }
    }
}




dockingConfig ={
    "id":"{}-dock".format(protein),
    "files":
        {
        'alphabetRec':
            {
                'folder':       "grid",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":     "-alphabet.grid"
            },
        'alphabetLig':
            {
                'folder':       "grid",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":     "-alphabet.grid"
            },
        'gridRec':
            {
                'folder':       "grid",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    "-grid.grid"
            },
        'gridLig':
            {
                'folder':       "grid",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    "-grid.grid"
            },
        'modesRec':
            {
                'folder':       "modes",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,modeType)
            },
        'modesLig':
            {
                'folder':       "modes",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,modeType)
            },
        'receptor':
            {
                'folder':       "",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    ".pdb"
            },
        'ligand':
            {
                'folder':       "",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    ".pdb"
            },
        'dof':
            {
                'folder':       "dof",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-dof.dat"
            },
        'dockingResult':
            {
                'folder':       "result",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-docking.dat"
            },
        'scoringResult':
            {
                'folder':       "result",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-scoring.dat"
            },
        'joinedModes':
            {
                'folder':       "modes",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-joinedModes-{}-{}.dat".format(numModes,modeType)
            }
        },
    "settings": {

        "attractToolPath": attractToolPath,
        "attractBinPath": attractBinPath,
        "basePath":basePath,

        "docking":{
            "in": 
            {
                "dof":"dof",
                "receptor": "receptor",
                "ligand": "ligand",
                "gridRec": "gridRec",
                "gridLig": "gridLig",
                "alphabetRec":"alphabetRec",
                "alphabetLig":"alphabetLig",
                "modesRec": "modesRec",
                "modesLig": "modesLig",
                
                },
            "out": 
            {"out": "dockingResult"},
            "numModesRec":numModesRec,
            "numModesLig":numModesLig,
            "dryRun": dry,
            "overwrite": overwrite,
            "GPUattractBinary":attractBinaryGPU,
            "verbose": verbose,
            "attractParFile": attractParFile,
        },
         "scoring":{
            "in": 
            {
                "dof":"dof",
                "receptor": "receptor",
                "ligand": "ligand",
                "gridRec": "gridRec",
                "gridLig": "gridLig",
                "joinedModes": "joinedModes"               
                },
            "out": 
            {"out": "scoringResult"},
            "numModesRec":numModesRec,
            "numModesLig":numModesLig,
            "dryRun": dry,
            "overwrite": overwrite,
            "attractBinary":attractBinary,
            "verbose": verbose,
            "attractParFile": attractParFile,

        }
    }
}












CONFIGURATION

# def outputFilesExist(self, setting):
    #     """returns true if at least one output file exists"""
    #     # if len(self.outputFiles) == 0:
    #     #     raise ValueError ("CONFIGURATION no outputfiles existing for config {} ".format(self.id))
    #     exists = False

    #     for file in self.files['output'].keys():
    #         if self.outputExists(file):
    #             exists = True
    #     return exists    

    # def getNotExistingInput(self):
    #     """returns all the inputfiles that do not exist"""
    #     notExisting = []
    #     for file in self.files['input'].keys():
    #         if not os.path.isfile(self.getInputFile(file)):
    #             notExisting.append(self.getInputFile(file))
    #     return notExisting

    #  def addOrChangeFile(self, identifier, filename):
    #     self.files['input'][identifier] = filename

    # def addOrChangeOutputFile(self, identifier, filename):
    #     self.files['output'][identifier] = filename


        # def addOrChangeSetting(self, settingName, settingValue):
        # self.settings[settingName] = settingValue

        #     def getExtension(self,key):
        # return self.extensions[key]
     # def getInputFolder(self):
    #     return os.path.join(self.basePath,self.folders['input'])

    # def getOutputFolder(self):
    #     return os.path.join(self.basePath,self.folders['output'])

    # def getOutputFile(self,output):
    #     return os.path.join(self.getOutputFolder(),self.files['output'][output])
    
    # def getInputFile(self,input):
    #     return os.path.join(self.getInputFolder(),self.files['input'][input])


    WORKER ///////////////////////

                    #print(self.inputQueue.qsize())
            #if self.counter  >= self.threshold:
                            #print("kill",self.name,self.finishCounter.value(),self.inputQueue.qsize())


        #while self.dorun:
            #time.sleep(0.01)




class addObj():
    def __init__(self,base):
        self.base = base


    def compute(self, item):
        return self.base + item



BENCHMARK



# singleConfigurationQueues_01 = []
# for i in range(len(singleConfiguration_01)+1):
#     singleConfigurationQueues_01.append(queue.Queue(BUF_SIZE))

# singleConfigurationQueues_02 = []
# for i in range(len(singleConfiguration_02)+1):
#     singleConfigurationQueues_02.append(queue.Queue(BUF_SIZE))

# pairConfigurationQueues = []
# for i in range(len(pairConfiguration)+1):
#     pairConfigurationQueues.append(queue.Queue(BUF_SIZE))

# runConfigurationQueues = []
# for i in range(len(runConfiguration)+1):
#     runConfigurationQueues.append(queue.Queue(BUF_SIZE))



# numThreads=1
# threads = []
# for i in range(len(configurators)):
#     c = ConsumerThread(configurator = configurators[i], threshold = 1,inputQueue=queues[i],resultQueue = queues[i+1], name='consumer{}'.format(i))
#     threads.append(c)


# for thread in threads:
#     thread.start()
# for thread in threads:
#     thread.join()
