# numModes = 5
# protein = "3MXW"
# chain = "A"
# dry = True
# verbose = True
# overwrite = False
# numModesRec = 5
# numModesLig = 5
# basePath = "/home/glenn/Documents/{}/{}".format(protein, protein)
# attractToolPath = "/home/glenn/Documents/attract/tools"
# attractBinPath = "/home/glenn/Documents/attract/bin"
# attractParFile = "/home/glenn/Documents/attract/attract.par"
# attractBinaryGPU = "/home/glenn/Documents/attract/bin"
# attractBinary = "/home/glenn/Documents/attract/bin"
# pythonBinary = "python2"

# modeType = "hin99"
# protType = "refe"

def getDefaultSingleSetting(protein, chain, protType, numModes,basePath, 
        modeType = "hin99",
        pythonBinary = "python2",
        attractToolPath = "$ATTRACTTOOLS",
        attractBinPath = "$ATTRACTDIR",
        attractParFile = "$ATTRACTDIR/../attract.par",
        dry = False,
        verbose = False,
        overwrite = False
    ):
    singleSetting = {
    "id":"{}{}-single".format(protein, chain ),
    "files":
        {
        'pdb':
            {
                'folder':       "",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    ".pdb"
            },
        'reduce':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-r.pdb"
            },
        'allAtom':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-aa.pdb"
            },
        'heavy':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-heavy.pdb"
            },
        'alphabet':
            {
                'folder':       "input/grid",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":     "-alphabet.grid"
            },
        'alphabetPartner':
            {
                'folder':       "input/grid",
                "name":         '{}{}-{}'.format(protein, "A" if chain == "B" else "B",protType),
                "extension":     "-alphabet.grid"
            },
        'grid':
            {
                'folder':       "input/grid",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-grid.grid"
            },
        'ca':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-ca.pdb"
            },
        'cut':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-cut.pdb"
            },
        'cutlog':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-cutlog.json"
            },
            
        'modes':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,modeType)
            },
        'secondary':
            {
                'folder':       "input/secondary",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    ".stride".format(numModes,modeType)
            },
        'settingsfile':
            {
                'folder':       "settings",
                "name":         '{}{}-{}-mn{}-mt{}'.format(protein, chain,protType,numModes,modeType),
                "extension":    "-settings.json"
            }
        },
    "settings": {
        "pythonBinary":pythonBinary,
        "attractToolPath": attractToolPath,
        "attractBinPath": attractBinPath,
        "attractParFile": attractParFile,
        "basePath":basePath,

        "modes":{
            "in": {"protein": "reduce" },
            "out": {"out": "modes" },
            "numModes":numModes,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "alphabet":{
            "in": {"protein": "reduce" },
            "out": {"out": "alphabet" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "reduce":{
            "in": {"protein": "pdb" },
            "out": {"out": "reduce" },
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "allAtom":{
            "in": {"protein": "pdb" },
            "out": {"out": "allAtom","mapping": "mapping" },
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "heavy":{
            "in": {"protein": "pdb"},
            "out": {"out": "heavy"},
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "grid":{
            "in": {"alphabet": "alphabetPartner", "protein": "reduce"},
            "out": {"out": "grid" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "CA":{
            "in": {"protein":"pdb"},
            "out":{'out':'ca'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
         "cut":{
            "in": {"pdb":"pdb"},
            "out":{'out':'cut','cutlog':'cutlog'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
            "cutoff":5.5
        },
         "secondary":{
            "in": {"pdb":"pdb"},
            "out":{'out':'secondary'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        'saveSettings': {
            'in':{
                },
            'out':{
                'out':'settingsfile'
            },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
            }
    }
}

    return singleSetting

def getDefaultPairSetting(benchmarkName,protein, protType, protTypeRef, numModesRec,numModesLig,basePath, 
        evScale = 1.0,
        modeType =  "hin-99",
        pythonBinary = "python2",
        attractToolPath = "$ATTRACTTOOLS",
        attractBinPath = "$ATTRACTDIR",
        attractParFile = "$ATTRACTDIR/../attract.par",
        attractBinary = "$ATTRACTDIR/attract",
        attractBinaryGPU = "/home/glenn/Documents/Masterarbeit/git/gpuATTRACT_2.0/AttractServer_RELEASE1",
        dry = False,
        verbose = False,
        overwrite = False,
        deviceIds = [0],
        interfaceCutoff = 5,
        numCollectStructures = 50
        ):
    pairSettings = {
    "id":"{}-pair".format(protein),
    "files":
        {
        'receptor':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    "-r.pdb"
            },
        'ligand':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    "-r.pdb"
            },
        'receptorRef':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, "A",protTypeRef),
                "extension":    "-r.pdb"
            },
        'ligandRef':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, "B",protTypeRef),
                "extension":    "-r.pdb"
            },
        'alphabetRec':
            {
                'folder':       "input/grid",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":     "-alphabet.grid"
            },
        'alphabetLig':
            {
                'folder':       "input/grid",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":     "-alphabet.grid"
            },
        'gridRec':
            {
                'folder':       "input/grid",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    "-grid.grid"
            },
        'gridLig':
            {
                'folder':       "input/grid",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    "-grid.grid"
            },
        'modesRec':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, "A",protType),
                "extension":    "-modes-{}-{}.dat".format(numModesRec,modeType)
            },
        'modesLig':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    "-modes-{}-{}.dat".format(numModesLig,modeType)
            },
        'joinedModes':
            {
                'folder':       "input/modes",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-joinedModes-r{}-l{}-{}.dat".format(numModesRec,numModesLig,modeType)
            },
        'translate':
            {
                'folder':       "input/dof",
                "name":         'translate',
                "extension":    ".dat"
            },
        'dof':
            {
                'folder':       "input/dof",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-dof.dat"
            },
        'dockingResult':
            {
                'folder':       "{}/result".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-docking.dat"
            },
        'scoringResult':
            {
                'folder':       "{}/result".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-scoring.dat"
            },
        'sortedResult':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-sorted.dat"
            },
        'deRedundantResult':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-deredundant.dat"
            },
        'demodeResult':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-demode.dat"
            },
        'topResult':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-top.dat"
            },
        'rmsd':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-rmsd.dat"
            },
        'irmsd':{
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-irmsd.dat"
            },
        'fnat':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-fnat.dat"
            },
        'collect':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-collect.pdb"
            },
        'interface':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-interface.json"
            },
        'settingsfile':
            {
                'folder':       "{}/settings".format(benchmarkName),
                "name":         'rec{}{}-{}-mn{}-mt_{}-lig{}{}-{}-mn{}-mt_{}-ev{}'.format(protein, "A",protType,numModesRec,modeType,protein, "B",protType,numModesLig,modeType, evScale),
                "extension":    ".json"
            }
        },
    "settings": 
        {
            "attractToolPath": attractToolPath,
            "attractBinPath": attractBinPath,
            "basePath":basePath,
            "pythonBinary":pythonBinary,

            "dof":{
                "in": {"receptor": "receptor","ligand": "ligand"  },
                "out": {"out": "dof", "translate":"translate" },
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
                "dofbinary": "/home/glenn/Documents/Masterarbeit/git/Attract_benchmark/tools/systsearch"
            },
            "joinModes":{
                "in": {"receptorModes": "modesRec","ligandModes": "modesLig"  },
                "out": {"out": "joinedModes" },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose
            },

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
                "evScale":evScale,
                "dryRun": dry,
                "overwrite": overwrite,
                "GPUattractBinary":attractBinaryGPU,
                "verbose": verbose,
                "attractParFile": attractParFile,
                "GPUdeviceIds": deviceIds
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
                "evScale":evScale,
                "dryRun": dry,
                "overwrite": overwrite,
                "attractBinary":attractBinary,
                "verbose": verbose,
                "attractParFile": attractParFile,

            },
            "sorting":
            {
                 "in": 
                {
                    "dockingResult":"dockingResult",
                    "scoringResult":"scoringResult"
                },
                "out": 
                {
                    "out": "sortedResult"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,

            },
             "top":
            {
                 "in": 
                {
                    "inputDof":"deRedundantResult",
                },
                "out": 
                {
                    "out": "topResult"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numDof": numCollectStructures
,

            },
            'deredundant':{
                "in": 
                {
                    "inputDof":"sortedResult"
                },
                "out": 
                {
                    "out": "deRedundantResult"
                },
                "dryRun": dry, "overwrite": overwrite,"verbose": verbose,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig,
            },
            "rmsd":{
                "in": 
                {
                    "inputDof":"deRedundantResult",
                    "ligand":"ligand",
                    "ligandRef":"ligandRef",
                    "modes":"joinedModes",
                    "receptorRef":"receptorRef"          
                    },
                "out": 
                {
                    "out": "rmsd"
                },
                "numModesRec":numModesRec,
                "numModesLig":numModesLig,
                "dryRun": dry,"overwrite": overwrite,"attractBinary":attractBinary,
                "verbose": verbose,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
                },
            "irmsd":
            {
                 "in": 
                {
                    "inputDof":"deRedundantResult",
                    "receptor":"receptor",
                    "receptorRef":"receptorRef",
                    "ligand":"ligand",
                    "ligandRef":"ligandRef",
                    "modes":"joinedModes",
                },
                "out": 
                {
                    "out": "irmsd"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
            },
            "fnat":{
                  "in": 
                {
                    "inputDof":"deRedundantResult",
                    "receptor":"receptor",
                    "receptorRef":"receptorRef",
                    "ligand":"ligand",
                    "ligandRef":"ligandRef",
                    "modes":"joinedModes"
                },
                "out": 
                {
                    "out": "fnat"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
            },
            "collect":{
                "in": 
                {
                    "inputDof":"topResult",
                    "receptor":"receptor",
                    "ligand":"ligand",
                    "modes":"joinedModes"
                },
                "out": 
                {
                    "out": "collect"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
            },
            'demode':
            {
                "in": 
                {
                    "inputDof":"topResult"
                },
                "out": 
                {
                    "out": "demodeResult"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
            },
            'interface':
            {
                "in": 
                {
                    # "receptor": "receptor",
                    # "ligand": "ligand",
                    'pdb': 'collect'
                },
                "out": 
                {
                    "out": "interface"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "cutoff": interfaceCutoff
            },
            'saveSettings':
            {
                'in':
                    {

                    },
                'out':{
                    'out':'settingsfile'
                },
                 "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
            }
        }
}

    return pairSettings





