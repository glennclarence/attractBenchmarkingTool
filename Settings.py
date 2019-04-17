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
        terminiCutoff = 5.5,
        checkInput = True
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
        'partner_bound':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,"refe"),
                "extension":    "-r.pdb"
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
        'prune':
            {
                'folder':       "input/pdb",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-pruned.pdb"
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
        'modes_manipulate':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,'manipulate_1')
            },
        'modes_manipulate_heavy':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-modes-{}-{}.dat".format(numModes,'manipulate_heavy_1')
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
                'folder':       "input/settings",
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
            },
        'bound_modes_heavy':
            {
                'folder':       "input/modes",
                "name":         '{}{}-{}'.format(protein, chain,protType),
                "extension":    "-modes-{}-{}.dat".format(1,'bound_heavy')
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
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "bound_mode":{
            'configurator':'bound_mode',
            "in": {"protein_bound": "partner_bound",'protein_unbound':'reduce' },
            "out": {"out": "bound_modes" },
            "numModes":numModes,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "bound_mode_heavy":{
            'configurator':'bound_mode',
            "in": {"protein_bound": "partner_bound_heavy",'protein_unbound':'heavy' },
            "out": {"out": "bound_modes_heavy" },
            "numModes":numModes,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "mode_evaluation":{
            'configurator':'mode_evaluation',
            "in": {"protein_bound": "partner_bound",'protein_unbound':'reduce','mode_file':'modes' },
            "out": {"out": "mode_evaluation" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
            "numModes":numModes
        },
        "mode_manipulate":{
            'configurator':'mode_manipulate',
            "in": {"protein": "reduce", 'mode_file':'modes', 'secondary':'secondary' },
            "out": {"out": "modes_manipulate" },
            'manipulate':['C'],
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "mode_manipulate_heavy":{
            'configurator':'mode_manipulate',
            "in": {"protein": "heavy", 'mode_file':'modes_heavy', 'secondary':'secondary' },
            "out": {"out": "modes_manipulate_heavy" },
            'manipulate':['C'],
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "modes_heavy":{
            'configurator':'modes',
            "in": {"protein": "heavy" },
            "out": {"out": "modes_heavy" },
            "numModes":numModes,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "alphabet":{
            'configurator':'alphabet',
            "in": {"protein": "reduce" },
            "out": {"out": "alphabet" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "reduce":{
            'configurator':'reduce',
            "in": {"protein": "allAtom" },
            "out": {"out": "reduce" },
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "prune":{
            'configurator':'prune',
            "in": {"pdb": "reduce" },
            "out": {"out": "prune" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "allAtom":{
            'configurator':'allAtom',
            "in": {"protein": "pdb" },
            "out": {"out": "allAtom","mapping": "mapping" },
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "heavy":{
            'configurator':'heavy',
            "in": {"protein": "allAtom"},
            "out": {"out": "heavy"},
            "chain": chain,
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "grid":{
            'configurator':'grid',
            "in": {"alphabet": "alphabetPartner", "protein": "reduce"},
            "out": {"out": "grid" },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "CA":{
            'configurator':'CA',
            "in": {"protein":"pdb"},
            "out":{'out':'ca'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "findtermini":{
            'configurator':'findtermini',
            "in": {"pdb":"pdb"},
            "out":{'out':'cutlog'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
            "cutoff": terminiCutoff
        },
        "cut":{
            'configurator':'cut',
            "in": {"pdb":"pdb", 'cutlog':'cutlog'},
            "out":{'out':'cut'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        "secondary":{
            'configurator':'secondary',
            "in": {"pdb":"pdb"},
            "out":{'out':'secondary'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        'superimpose':{
            'configurator':'superimpose',
            "in": {"pdb":"pdb", "refpdb": 'refpdb'},
            "out":{'out':'superimpose'},
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
        },
        'saveSettings': {
            'configurator':'saveSettings',
            'in':{
                },
            'out':{
                'out':'settingsfile'
            },
            "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
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
        scoringCutoff = 50,
        num_dof_eval = 100,
        checkInput = True,
        dofBinary = "/home/glenn/Documents/Masterarbeit/git/Attract_benchmark/tools/systsearch"
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
                "extension":    "-joinedModes-r{}-l{}-{}.dat".format(20,20,modeType)
            },
        'joinedModes_heavy':
            {
                'folder':       "input/modes",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-joinedModes-heavy-r{}-l{}-{}.dat".format(20,20,modeType)
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
        'dof_test':
            {
                'folder':       "input/dof",
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-dof_test.dat"
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
        'rmsd_unsorted':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-rmsd_unsorted.dat"
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
        'collect_low':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-collect_low.pdb"
            },
        'collect_high':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-collect_high.pdb"
            },
        'interface':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-interface.json"
            },
        'interface_low':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-interface_low.json"
            },
        'interface_high':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-interface_high.json"
            },
        'dof_evaluation':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-dof_eval.json"
            },
        'dofs_extractedLow':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-dof_extractedLow.dat"
            },
        'dofs_extractedHigh':
            {
                'folder':       "{}/analysis".format(benchmarkName),
                "name":         '{}-{}'.format(protein,protType),
                "extension":    "-dof_extractedHigh.dat"
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
                "dryRun": dry, "overwrite": overwrite, "verbose": verbose,'checkInput':checkInput,
                "dofbinary": dofBinary
            },
            "dof_test":{
                'configurator':'dof_test',
                "in": { },
                "out": {"out": "dof_test" },
                "dryRun": dry, "overwrite": overwrite, "verbose": verbose,'checkInput':checkInput,
            },
            "joinModes":{
                'configurator':'joinModes',
                "in": {"receptorModes": "modesRec","ligandModes": "modesLig"  },
                "out": {"out": "joinedModes" },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':False,
            },

            "joinModes_heavy":{
                'configurator':'joinModes',
                "in": {"receptorModes": "modesRec_heavy","ligandModes": "modesLig_heavy"  },
                "out": {"out": "joinedModes_heavy" },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':False,
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
                "verbose": verbose,'checkInput':True,
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
                "verbose": verbose,'checkInput':True,
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
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,

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
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,

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
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
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
                "dryRun": dry, "overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
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
                "verbose": verbose,'checkInput':False,

                },
            "rmsd_unsorted":{
                'configurator':'rmsd',
                "in": 
                {
                    "inputDof":"filled",
                    "ligand":"ligand",
                    "ligandRef":"ligandRef",
                    "modes":"joinedModes",
                    "receptorRef":"receptorRef"          
                    },
                "out": 
                {
                    "out": "rmsd_unsorted"
                },
                "numModesRec":numModesRec,
                "numModesLig":numModesLig,
                "dryRun": dry,"overwrite": overwrite,"attractBinary":attractBinary,
                "verbose": verbose,'checkInput':False,

                },
            "rmsd_nomodes":{
                'configurator':'rmsd',
                "in": 
                {
                    "inputDof":"demodeResult",
                    "ligand":"ligand",
                    "ligandRef":"ligandRef",
                    "receptorRef":"receptorRef"          
                    },
                "out": 
                {
                    "out": "rmsd_nomodes"
                },
                "numModesRec":0,
                "numModesLig":0,
                "dryRun": dry,"overwrite": overwrite,"attractBinary":attractBinary,
                "verbose": verbose,'checkInput':checkInput,
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
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':False,
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
                },
                "out": 
                {
                    "out": "irmsd_nomodes"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
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
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':False,
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
                },
                "out": 
                {
                    "out": "fnat_nomodes"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
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
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':False,
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
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
            },
            'interface':
            {
                'configurator':'interface',
                "in": 
                {
                    'receptor': 'receptor',
                    'ligand': 'ligand',
                     'receptorSec': 'receptorSec',
                    'ligandSec': 'ligandSec',
                    'pdb': 'collect'
                },
                "out": 
                {
                    "out": "interface"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
                "cutoff": interfaceCutoff
            },
            'interface_low':
            {
                'configurator':'interface',
                "in": 
                {
                    'receptor': 'receptor',
                    'ligand': 'ligand',
                    'receptorSec': 'receptorSec',
                    'ligandSec': 'ligandSec',
                    'pdb': 'collect_low'
                },
                "out": 
                {
                    "out": "interface_low"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
                "cutoff": interfaceCutoff
            },
            'interface_high':
            {
                'configurator':'interface',
                "in": 
                {
                    'receptor': 'receptor',
                    'ligand': 'ligand',
                     'receptorSec': 'receptorSec',
                    'ligandSec': 'ligandSec',
                    'pdb': 'collect_high'
                },
                "out": 
                {
                    "out": "interface_high"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
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
                 "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
            },
            'dof_evaluation': {
                'configurator':'dof_evaluation',
                'in':{"input_dof":"deRedundantResult",
                'mode_evaluation_rec':'mode_evaluation_rec',
                'mode_evaluation_lig':'mode_evaluation_lig'
                    },
                'out':{
                    'out':'dof_evaluation'
                },
                'numModesRec': numModesRec,
                'numModesLig': numModesLig,
                'num_eval':num_dof_eval,
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
            },
            'dof_extractionLow': {
                'configurator':'dof_extraction',
                'in':{"inpuf_dofs":"deRedundantResult",
                        'rmsd_file':"rmsd"
                    },
                'out':{
                    'out':'dofs_extractedLow'
                },
                'maxDof': 50,
                'threshold': 10,
                'type':'max',
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
            },
            'dof_extractionHigh': {
                'configurator':'dof_extraction',
                'in':{"inpuf_dofs":"deRedundantResult",
                        'rmsd_file':"rmsd"
                    },
                'out':{
                    'out':'dofs_extractedHigh'
                },
                'maxDof': 50,
                'threshold':  25,
                'type':'min',
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':checkInput,
            },
            "collect_Low":{
                'configurator':'collect',
                "in": 
                {
                    "inputDof":"dofs_extractedLow",
                    "receptor":"receptor",
                    "ligand":"ligand",
                    "modes":"joinedModes"
                },
                "out": 
                {
                    "out": "collect_low"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':False,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
            },
            "collect_High":{
                'configurator':'collect',
                "in": 
                {
                    "inputDof":"dofs_extractedHigh",
                    "receptor":"receptor",
                    "ligand":"ligand",
                    "modes":"joinedModes"
                },
                "out": 
                {
                    "out": "collect_high"
                },
                "dryRun": dry,"overwrite": overwrite,"verbose": verbose,'checkInput':False,
                "numModesRec":numModesRec,
                "numModesLig":numModesLig
            },
        }
}

    return pairSettings





