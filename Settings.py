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
        attractParFile = "$ATTRACTDIR/..",
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
                'folder':       "pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-r.pdb"
            },
        'allAtom':
            {
                'folder':       "pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-aa.pdb"
            },
        'heavy':
            {
                'folder':       "pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-heavy.pdb"
            },
        'alphabet':
            {
                'folder':       "grid",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":     "-alphabet.grid"
            },
        'alphabetPartner':
            {
                'folder':       "grid",
                "name":         '{}{}-{}'.format(protein, "A" if chain == "B" else "B",protType),
                "extension":     "-alphabet.grid"
            },
        'grid':
            {
                'folder':       "grid",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-grid.grid"
            },
        'ca':
            {
                'folder':       "pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-ca.pdb"
            },
        'cut':
            {
                'folder':       "pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-cut.pdb"
            },
        'modes':
            {
                'folder':       "modes",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,modeType)
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
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        },
        "alphabet":{
            "in": {"protein": "reduce" },
            "out": {"out": "alphabet" },
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        },
        "reduce":{
            "in": {"protein": "pdb" },
            "out": {"out": "reduce" },
            "chain": chain,
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        },
        "allAtom":{
            "in": {"protein": "pdb" },
            "out": {"out": "allAtom","mapping": "mapping" },
            "chain": chain,
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        },
        "heavy":{
            "in": {"protein": "pdb"},
            "out": {"out": "heavy"},
            "chain": chain,
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        },
        "grid":{
            "in": {"alphabet": "alphabetPartner", "protein": "reduce"},
            "out": {"out": "grid" },
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        },
        "CA":{
            "in": {"protein":"pdb"},
            "out":{'out':'ca'},
            "dryRun": dry,
            "overwrite": overwrite,
            "verbose": verbose
        }
    }
}

    return singleSetting

def getDefaultPairSetting(benchmarkName,protein, protType, protTypeRef, numModesRec,numModesLig,basePath, 
        modeType =  "hin-99",
        pythonBinary = "python2",
        attractToolPath = "$ATTRACTTOOLS",
        attractBinPath = "$ATTRACTDIR",
        attractParFile = "$ATTRACTDIR/..",
        attractBinary = "$ATTRACTDIR/attract",
        attractBinaryGPU = "$ATTRACTDIR/AttractServer",
        dry = False,
        verbose = False,
        overwrite = False
        ):
    pairSettings = {
    "id":"{}-pair".format(protein),
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
        'receptorRef':
            {
                'folder':       "",
                "name":         '{}{}-{}'.format(protein, "A",protTypeRef),
                "extension":    ".pdb"
            },
        'ligandRef':
            {
                'folder':       "",
                "name":         '{}{}-{}'.format(protein, "B",protTypeRef),
                "extension":    ".pdb"
            },
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
                "extension":    "-modes-{}-{}.dat".format(numModesRec,modeType)
            },
        'modesLig':
            {
                'folder':       "modes",
                "name":         '{}{}-{}'.format(protein, "B",protType),
                "extension":    "-modes-{}-{}.dat".format(numModesLig,modeType)
            },
        'joinedModes':
            {
                'folder':       "modes",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-joinedModes-r{}-l{}-{}.dat".format(numModesRec,numModesLig,modeType)
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
                "dofbinary": "/home/glenn/Documents/Masterarbeit/git/Attract_benchmark/tools/systseach"
            },
            "joinModes":{
                "in": {"receptorModes": "modesRec","ligandModes": "modesLig"  },
                "out": {"out": "joinedModes" },
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose
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
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,

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
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
                "numDof":1000,

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
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
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
                "dryRun": dry,
                "overwrite": overwrite,
                "attractBinary":attractBinary,
                "verbose": verbose,
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
                    "modes":"joinedModes"
                },
                "out": 
                {
                    "out": "irmsd"
                },
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
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
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
            },
            "collect":{
                "in": 
                {
                    "inputDof":"deRedundantResult",
                    "receptor":"receptor",
                    "ligand":"ligand",
                    "modes":"joinedModes"
                },
                "out": 
                {
                    "out": "collect"
                },
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
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
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
            }
        }
}

    return pairSettings





