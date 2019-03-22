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
        overwrite = False,
        terminiCutoff = 5.5
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
        'refpdb':
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
        'mapping':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-mapping.pdb"
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
        'modes_heavy':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-modes-heavy-{}-{}.dat".format(numModes,modeType)
            },
        'secondary':
            {
                'folder':       "input/secondary",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    ".stride"
            },
        'superimpose':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-super.pdb"
            },
        'settingsfile':
            {
                'folder':       "",
                "name":         '{}{}-{}-mn{}-mt_{}'.format(protein, chain,protType,numModes,modeType),
                "extension":    "-settings.json"
            },
        'mode_evaluation':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}-mn{}-mt_{}'.format(protein, chain,protType,numModes,modeType),
                "extension":    "-eval.json"
            },
        'bound_modes':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-modes-{}-{}.dat".format(1,'bound')
            }
        },
    "settings": {
        "pythonBinary":pythonBinary,
        "attractToolPath": attractToolPath,
        "attractBinPath": attractBinPath,
        "attractParFile": attractParFile,
        "basePath":basePath,

        "modes":{
            'configurator':'modes',
            "in": {"protein": "reduce" },
            "out": {"out": "modes" },
            "numModes":numModes,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "bound_mode":{
            'configurator':'bound_mode',
            "in": {"protein_bound": "partner_bound",'protein_unbound':'reduce' },
            "out": {"out": "bound_modes" },
            "numModes":numModes,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "mode_evaluation":{
            'configurator':'mode_evaluation',
            "in": {"protein_bound": "partner_bound",'protein_unbound':'reduce','mode_file':'modes' },
            "out": {"out": "mode_evaluation" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "mode_manipulate":{
            'configurator':'mode_manipulate',
            "in": {"protein": "reduce", 'mode_file':'modes' },
            "out": {"out": "mode_evaluation" },
            'manipulate':['T'],
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "modes_heavy":{
            'configurator':'modes',
            "in": {"protein": "heavy" },
            "out": {"out": "modes_heavy" },
            "numModes":numModes,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "alphabet":{
            'configurator':'alphabet',
            "in": {"protein": "reduce" },
            "out": {"out": "alphabet" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "reduce":{
            'configurator':'reduce',
            "in": {"protein": "allAtom" },
            "out": {"out": "reduce" },
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "allAtom":{
            'configurator':'allAtom',
            "in": {"protein": "pdb" },
            "out": {"out": "allAtom","mapping": "mapping" },
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "heavy":{
            'configurator':'heavy',
            "in": {"protein": "allAtom"},
            "out": {"out": "heavy"},
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "grid":{
            'configurator':'grid',
            "in": {"alphabet": "alphabetPartner", "protein": "reduce"},
            "out": {"out": "grid" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "CA":{
            'configurator':'CA',
            "in": {"protein":"pdb"},
            "out":{'out':'ca'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "findtermini":{
            'configurator':'findtermini',
            "in": {"pdb":"pdb"},
            "out":{'out':'cutlog'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
            "cutoff": terminiCutoff
        },
        "cut":{
            'configurator':'cut',
            "in": {"pdb":"pdb", 'cutlog':'cutlog'},
            "out":{'out':'cut'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        "secondary":{
            'configurator':'secondary',
            "in": {"pdb":"pdb"},
            "out":{'out':'secondary'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        'superimpose':{
            'configurator':'superimpose',
            "in": {"pdb":"pdb", "refpdb": 'refpdb'},
            "out":{'out':'superimpose'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose
        },
        'saveSettings': {
            'configurator':'saveSettings',
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
        numCollectStructures = 50,
        scoringCutoff = 50
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
        'joinedModes_heavy':
            {
                'folder':       "input/modes",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-joinedModes-heavy-r{}-l{}-{}.dat".format(numModesRec,numModesLig,modeType)
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
        'filled':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-filled.dat"
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
        'rmsd_nomodes':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-rmsd_nomodes.dat"
            },
        'irmsd':{
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-irmsd.dat"
            },
        'irmsd_nomodes':{
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-irmsd-nomodes.dat"
            },
        'fnat':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-fnat.dat"
            },
        'fnat_nomodes':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-fnat-nomodes.dat"
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
            "attractToolPath":  attractToolPath,
            "attractBinPath":   attractBinPath,
            "basePath":         basePath,
            "pythonBinary":     pythonBinary,

            "dof":{
                'configurator':'dof',
                "in": {"receptor": "receptor","ligand": "ligand"  },
                "out": {"out": "dof", "translate":"translate" },
                "dryRun": dry,
                "overwrite": overwrite,
                "verbose": verbose,
                "dofbinary": "/home/glenn/Documents/Masterarbeit/git/Attract_benchmark/tools/systsearch"
            },
            "joinModes":{
                'configurator':'joinModes',
                "in": {"receptorModes": "modesRec","ligandModes": "modesLig"  },
                "out": {"out": "joinedModes" },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose
            },

            "joinModes_heavy":{
                'configurator':'joinModes',
                "in": {"receptorModes": "modesRec_heavy","ligandModes": "modesLig_heavy"  },
                "out": {"out": "joinedModes_heavy" },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose
            },

            "docking":{
                'configurator':'docking',
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
                'configurator':'scoring',
                "in": 
                {
                    "dof":"dockingResult",
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
                "cutoff":scoringCutoff

            },
            "fill_energy":
            {
                'configurator':'fill_energy',
                "in": 
                {
                    "dockingResult":"dockingResult",
                    "scoringResult":"scoringResult"
                },
                "out": 
                {
                    "out": "filled"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,

            },
            "sorting":
            {
                'configurator':'sorting',
                "in": 
                {
                    "input_dof":"filled",
                },
                "out": 
                {
                    "out": "sortedResult"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,

            },
             "top":
            {
                'configurator':'top',
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
                'configurator':'deredundant',
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
                'configurator':'rmsd',
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

                },
            "rmsd_nomodes":{
                'configurator':'rmsd',
                "in": 
                {
                    "inputDof":"demodeResult",
                    "ligand":"ligand",
                    "ligandRef":"ligandRef",
                    "modes":"joinedModes",
                    "receptorRef":"receptorRef"          
                    },
                "out": 
                {
                    "out": "rmsd_nomodes"
                },
                "numModesRec":0,
                "numModesLig":0,
                "dryRun": dry,"overwrite": overwrite,"attractBinary":attractBinary,
                "verbose": verbose,
                },
            "irmsd":
            {
                'configurator':'irmsd',
                "in": 
                {
                    "inputDof":"deRedundantResult",
                    "receptor":"receptor_heavy",
                    "receptorRef":"receptorRef_heavy",
                    "ligand":"ligand_heavy",
                    "ligandRef":"ligandRef_heavy",
                    "modes":"joinedModes_heavy",
                },
                "out": 
                {
                    "out": "irmsd"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
            },
            "irmsd_nomodes":
            {
                'configurator':'irmsd',
                "in": 
                {
                    "inputDof":"demodeResult",
                    "receptor":"receptor_heavy",
                    "receptorRef":"receptorRef_heavy",
                    "ligand":"ligand_heavy",
                    "ligandRef":"ligandRef_heavy",
                    "modes":"joinedModes_heavy",
                },
                "out": 
                {
                    "out": "irmsd_nomodes"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numModesRec":0,
                "numModesLig":0
            },
            "fnat":{
                'configurator':'fnat',
                "in": 
                {
                    "inputDof":"deRedundantResult",
                    "receptor":"receptor_heavy",
                    "receptorRef":"receptorRef_heavy",
                    "ligand":"ligand_heavy",
                    "ligandRef":"ligandRef_heavy",
                    "modes":"joinedModes_heavy"
                },
                "out": 
                {
                    "out": "fnat"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
            },
            "fnat_nomodes":{
                'configurator':'fnat',
                "in": 
                {
                    "inputDof":"demodeResult",
                    "receptor":"receptor_heavy",
                    "receptorRef":"receptorRef_heavy",
                    "ligand":"ligand_heavy",
                    "ligandRef":"ligandRef_heavy",
                    "modes":"joinedModes_heavy"
                },
                "out": 
                {
                    "out": "fnat_nomodes"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
                "numModesRec":0,
                "numModesLig":0
            },
            "collect":{
                'configurator':'collect',
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
                'configurator':'demode',
                "in": 
                {
                    "inputDof":"deRedundantResult"
                },
                "out": 
                {
                    "out": "demodeResult"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,
            },
            'interface':
            {
                'configurator':'joinModes',
                "in": 
                {

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
                'configurator':'saveSettings',
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





