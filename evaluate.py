import utils as utils
import numpy as np
from collections import Counter
import os
import pandas as pd
import json
import glob

def getCapriCriteria(lrmsd, irmsd, fnat):
    if fnat >= 0.5 and( lrmsd <= 1 or irmsd <= 1):
        return 3
    elif (fnat >= 0.5 and lrmsd > 1 and irmsd > 1) or ( fnat >= 0.3 and fnat < 0.5 and lrmsd <= 5 and irmsd <= 2):
        return 2
    elif (fnat >=0.3 and lrmsd  > 5 and irmsd >2)  or (fnat >= 0.1 and fnat < 0.3 and lrmsd <= 10 and irmsd <= 4 ):
        return 1
    elif fnat < 0.1 or (lrmsd > 10 and irmsd > 4):
        return 0
    else:
        return 0

def readRMSD(filename):
    lines = utils.readFileToList(filename)
    if len(lines) == 0:
        return np.nan
    rmsd = np.zeros(len(lines))
    for idx,line in enumerate(lines):
        line_split = line.split()
        if len(line_split) == 2:
            rmsd[idx] = float( line_split[1])
    return rmsd

def readFNAT(filename):
    lines = utils.readFileToList(filename)
    if len(lines) == 0:
        return np.nan
    fnat = np.zeros(len(lines))
    for idx,line in enumerate(lines):
        line_split = line.split()
        if len(line_split) == 1:
            fnat[idx] = float( line_split[0])
    return fnat


def readEnergies( filename ):
    energies = []
    lines = utils.readFileToList(filename)
    if( len(lines ) > 0):
        count = 0
        idx = 1
        for line in lines:
            if "Energy" in line:
                energy = float(line.split()[2])
                energies.append( energy)
    return np.asarray( energies )

def getCAPRIRating(lrmsd_list, irmsd_list, fnat_list):
    ratings = np.zeros(len(lrmsd_list))
    for i,(lrmsd, irmsd, fnat) in enumerate(zip(lrmsd_list, irmsd_list, fnat_list)):
        ratings[i] = getCapriCriteria(lrmsd, irmsd, fnat)
    return ratings

def getCapriOccurances(capri_list):
    return Counter(capri_list)


def getBenchmarkpath(base,protein,benchmark_name):
    return os.path.join(base, os.path.join(protein, benchmark_name))


def evaluateBenchmarks( base_path, protein_list, benchmark_names):
    result = []
    for  bm in benchmark_names:
        for protein in protein_list:
            print("processing protein", protein ,"in benchmark", bm)
            path = getBenchmarkpath(base_path, protein, bm)
            try:
                # rmsd_filename = os.path.join(path,"analysis/{}-{}-rmsd.dat".format(protein, proteinType ))
                # irmsd_filename = os.path.join(path,"analysis/{}-{}-irmsd.dat".format(protein, proteinType ))
                # fnat_filename = os.path.join(path,"analysis/{}-{}-fnat.dat".format(protein, proteinType ))
                # energy_filename = os.path.join(path,"analysis/{}-{}-deredundant.dat".format(protein, proteinType ))
                if len(glob.glob( os.path.join(path,"analysis/*-rmsd.dat"))[0]) > 0:
                    rmsd_filename =glob.glob( os.path.join(path,"analysis/*-rmsd.dat"))[0]
                    irmsd_filename = glob.glob(os.path.join(path,"analysis/*-irmsd.dat"))[0]
                    fnat_filename = glob.glob(os.path.join(path,"analysis/*-fnat.dat"))[0]
                    energy_filename = glob.glob(os.path.join(path,"analysis/*-deredundant.dat"))[0]
                else:
                    continue
                rmsd = readRMSD(rmsd_filename)
                irmsd = readRMSD(irmsd_filename)
                #energies = readEnergies(energy_filename)
                fnat = readFNAT(fnat_filename)
                if rmsd is np.nan or irmsd is np.nan or fnat is np.nan: 
                    print("could not evaluate protein" , protein, " at benchmark " , bm)
                    continue
                capri_list = getCAPRIRating(rmsd, irmsd, fnat)


                ratings_10 = getCapriOccurances(capri_list[:10])
                ratings_50 = getCapriOccurances(capri_list[:50])
                ratings_100 = getCapriOccurances(capri_list[:100])
                ratings_1000 = getCapriOccurances(capri_list[:1000])

                modeType ="none"
                if "boundModes" in bm:
                    modeType = "boundModes"
                if "brownian" in bm:
                    modeType = "brownian"
                if "hin99" in bm:
                    modeType = "hin99"
                bound = "bound" in bm and not "boundModes" in bm

                row = {                "protein":protein,
                    "cut": "cut" in bm,
                    "numModesRec": int(bm.split('_')[2][2:]),
                    "numModesLig": int(bm.split('_')[3][2:]),
                    "scale": float(bm.split('_')[4][1:].replace('p','.')),
                    "bound": bound,
                    "modetype": modeType,
                    "num1_10":ratings_10[1],
                    "num1_50":ratings_50[1],
                    "num1_100":ratings_100[1],
                    "num1_1000":ratings_1000[1],
                    "num2_10":ratings_10[2],
                    "num2_50":ratings_50[2],
                    "num2_100":ratings_100[2],
                    "num2_1000":ratings_1000[2],
                    "num3_10":ratings_10[3],
                    "num3_50":ratings_50[3],
                    "num3_100":ratings_100[3],
                    "num3_1000":ratings_1000[3],
                    'bm':bm
                }
                result.append(row)
            except:
                pass
    return pd.DataFrame(result)



def evaluateSingle( base_path, protein_list, benchmark_names):
    result =[]
    for  bm in benchmark_names:
        for protein in protein_list:
            #try:
                path = getBenchmarkpath(base_path, protein, bm)
                # rmsd_filename = os.path.join(path,"analysis/{}-{}-rmsd.dat".format(protein, proteinType ))
                # irmsd_filename = os.path.join(path,"analysis/{}-{}-irmsd.dat".format(protein, proteinType ))
                # fnat_filename = os.path.join(path,"analysis/{}-{}-fnat.dat".format(protein, proteinType ))
                # energy_filename = os.path.join(path,"analysis/{}-{}-deredundant.dat".format(protein, proteinType ))
                # dof_filename = os.path.join(path,"analysis/{}-{}-dof_eval.json".format(protein, proteinType ))

                if len(glob.glob( os.path.join(path,"analysis/*-rmsd.dat"))[0]) > 0:
                    rmsd_filename =glob.glob( os.path.join(path,"analysis/*-rmsd.dat"))[0]
                    irmsd_filename = glob.glob(os.path.join(path,"analysis/*-irmsd.dat"))[0]
                    fnat_filename = glob.glob(os.path.join(path,"analysis/*-fnat.dat"))[0]
                    energy_filename = glob.glob(os.path.join(path,"analysis/*-deredundant.dat"))[0]
                    dofl = glob.glob(os.path.join(path,"analysis/*-dof_eval.json"))
                    if len(dofl) > 0:
                        dof_filename = dofl[0]

                rmsd = readRMSD(rmsd_filename)
                irmsd = readRMSD(irmsd_filename)
                energies = readEnergies(energy_filename)
                fnat = readFNAT(fnat_filename)
                if rmsd is np.nan or irmsd is np.nan or fnat is np.nan: 
                    print("could not evaluate protein" , protein, " at benchmark " , bm)
                    continue
                capri_list = getCAPRIRating(rmsd, irmsd, fnat)

                
                modeType ="none"
                if "boundModes" in bm:
                    modeType = "boundModes"
                if "brownian" in bm:
                    modeType = "brownian"
                if "hin99" in bm:
                    modeType = "hin99"
                bound = "bound" in bm and not "boundModes" in bm

                row = {                "protein":protein,
                    "cut": "cut" in bm,
                    "numModesRec": int(bm.split('_')[2][2:]),
                    "numModesLig": int(bm.split('_')[3][2:]),
                    "scale": float(bm.split('_')[4][1:].replace('p','.')),
                    "bound": bound,
                    "modetype": modeType,
                    "capri":capri_list[0],
                    "irmsd":irmsd[0],
                    "lrmsd":rmsd[0],
                    "fnat":fnat[0],
                    "energy" : energies[0],

                    'bm':bm
                }
                if len(dofl) > 0:
                    dof_eval = json.load(open(dof_filename))["1"]
                    for modeIdx,value  in dof_eval['rec'].items():
                        row['rec_ratiomode_{}'.format(modeIdx)] = value['ratio']
                    for modeIdx,value  in dof_eval['lig'].items():
                        row['lig_ratiomode_{}'.format(modeIdx)] = value['ratio']
                result.append(row)
            #except:
            #    print("couldn't get protein", protein,bm)
            #    pass
    return pd.DataFrame(result)



def interface( benchmark_names, base_path, protein_list, protein_types):

    result_low_rec_sec = {'C': 0, 'E': 0, 'B': 0, 'T': 0, 'H': 0, 'G': 0, 'b': 0}
    result_high_rec_sec = {'C': 0, 'E': 0, 'B': 0, 'T': 0, 'H': 0, 'G': 0, 'b': 0}

    result_low_lig_sec = {'C': 0, 'E': 0, 'B': 0, 'T': 0, 'H': 0, 'G': 0, 'b': 0}
    result_high_lig_sec = {'C': 0, 'E': 0, 'B': 0, 'T': 0, 'H': 0, 'G': 0, 'b': 0}


    result_low_rec_aa =  {'LYS': 0, 'PRO': 0,'ILE': 0,'TRP': 0,'GLU': 0,'GLN': 0,'GLY': 0,'SER': 0,'PHE': 0,'HIS': 0,'TYR': 0,'LEU': 0,'ASP': 0,'ASN': 0,'ARG': 0,'THR': 0,'ALA': 0,'CYS': 0,'VAL': 0,'MET': 0}
    result_high_rec_aa = {'LYS': 0, 'PRO': 0,'ILE': 0,'TRP': 0,'GLU': 0,'GLN': 0,'GLY': 0,'SER': 0,'PHE': 0,'HIS': 0,'TYR': 0,'LEU': 0,'ASP': 0,'ASN': 0,'ARG': 0,'THR': 0,'ALA': 0,'CYS': 0,'VAL': 0,'MET': 0}
    result_low_lig_aa =  {'LYS': 0, 'PRO': 0,'ILE': 0,'TRP': 0,'GLU': 0,'GLN': 0,'GLY': 0,'SER': 0,'PHE': 0,'HIS': 0,'TYR': 0,'LEU': 0,'ASP': 0,'ASN': 0,'ARG': 0,'THR': 0,'ALA': 0,'CYS': 0,'VAL': 0,'MET': 0}
    result_high_lig_aa = {'LYS': 0, 'PRO': 0,'ILE': 0,'TRP': 0,'GLU': 0,'GLN': 0,'GLY': 0,'SER': 0,'PHE': 0,'HIS': 0,'TYR': 0,'LEU': 0,'ASP': 0,'ASN': 0,'ARG': 0,'THR': 0,'ALA': 0,'CYS': 0,'VAL': 0,'MET': 0}

    area_low_rec = 0
    area_low_lig = 0

    area_high_rec = 0
    area_high_lig = 0


    count_low = 0
    count_high = 0
    for proteinType, bm in zip(protein_types,benchmark_names):
            for protein in protein_list:
                try:
                    path_low = getBenchmarkpath(base_path, protein, bm) + "/analysis/{}-{}-interface_high.json".format( protein, proteinType)
                    path_high = getBenchmarkpath(base_path, protein, bm) + "/analysis/{}-{}-interface_low.json".format( protein, proteinType)
                    data_low = json.load(open(path_low))
                    data_high = json.load(open(path_high))
                    for i in data_low['interfaces']:
                        count_low += 1
                       
                        #rec
                        area_low_rec += i['areaRec']
                        for key, val in i['countSecRec'].items():
                            result_low_rec_sec[key] += val
                        for key, val in i['AARecCount'].items():
                            result_low_rec_aa[key] += val
                        #lig
                        area_low_lig += i['areaLig']
                        for key, val in i['countSecLig'].items():
                            result_low_lig_sec[key] += val
                        for key, val in i['AALigCount'].items():
                            result_low_lig_aa[key] += val




                    for i in data_high['interfaces']:
                        count_high += 1
                        #rec
                        area_high_rec += i['areaRec']
                        for key, val in i['countSecRec'].items():
                            result_high_rec_sec[key] += val
                        for key, val in i['AARecCount'].items():
                            result_high_rec_aa[key] += val
                        #lig
                        area_high_lig += i['areaLig']
                        for key, val in i['countSecLig'].items():
                            result_high_lig_sec[key] += val
                        for key, val in i['AALigCount'].items():
                            result_high_lig_aa[key] += val

                except:
                    print('unable to read', protein)
                    pass

    count_low = float(count_low)
    count_high = float(count_high)
    area_high_lig /= count_high
    area_high_rec /= count_high
    area_low_rec  /= count_low
    area_low_lig  /= count_low


    for key in result_low_rec_sec.keys():
        result_low_rec_sec[key] /= float(count_low)
    for key in result_low_rec_aa.keys():
        result_low_rec_aa[key] /= float(count_low)
    for key in result_low_lig_sec.keys():
        result_low_lig_sec[key] /= float(count_low)
    for key in result_low_lig_aa.keys():
        result_low_lig_aa[key] /= float(count_low)

    for key in result_high_rec_sec.keys():
        result_high_rec_sec[key] /= float(count_high)
    for key in result_high_rec_aa.keys():
        result_high_rec_aa[key] /= float(count_high)
    for key in result_high_lig_sec.keys():
        result_high_lig_sec[key] /= float(count_high)
    for key in result_high_lig_aa.keys():
        result_high_lig_aa[key] /= float(count_high)

    return result_low_rec_sec,result_low_rec_aa,result_low_lig_sec,result_low_lig_aa,result_high_rec_sec,result_high_rec_aa,result_high_lig_sec,result_high_lig_aa,area_high_lig,area_high_rec,area_low_rec ,area_low_lig 



def dofs( benchmark_names, base_path, protein_list, protein_types):
    count_low_rec = 0
    count_high_rec = 0
    count_low_lig = 0
    count_high_lig = 0
    for proteinType, bm in zip(protein_types,benchmark_names):
            for protein in protein_list:
                    path_low = getBenchmarkpath(base_path, protein, bm) + "/analysis/{}-{}-dof_eval_low.json".format( protein, proteinType)
                    path_high = getBenchmarkpath(base_path, protein, bm) + "/analysis/{}-{}-dof_eval_high.json".format( protein, proteinType)
                    evaluation_low  = json.load(open(path_low, 'r'))
                    mean_low_rec = 0
                    mean_high_rec = 0
                    mean_low_lig = 0
                    mean_high_lig = 0
                    res  ={'rec_low':[],'lig_low':[],'rec_high':[],'lig_high':[]}
                    for key,i in evaluation_low.items():
                        for key1,k in i['rec'].items():
                            mean_low_rec += np.float64(k['dof']) / k['mode']
                            res['rec_low'].append(np.float64(k['dof']) / k['mode'])
                            count_low_rec += 1
                        for key1,k in i['lig'].items():
                            mean_low_lig += np.float64(k['dof'])/ k['mode']
                            count_low_lig += 1
                            res['lig_low'].append(np.float64(k['dof']) / k['mode'])


                    evaluation_high  = json.load(open(path_high, 'r'))
                    
                    for key,i in evaluation_high.items():
                        for key1,k in i['rec'].items():
                            mean_high_rec += np.float64(k['dof'])/ k['mode']
                            res['rec_high'].append(np.float64(k['dof']) / k['mode'])

                            count_high_rec += 1
                        for key1,k in i['lig'].items():
                            mean_high_lig += np.float64(k['dof'])/ k['mode']
                            count_high_lig += 1
                            res['lig_high'].append(np.float64(k['dof']) / k['mode'])
                    print(protein, res['rec_low'],res['lig_low'])
    return mean_low_rec/ count_low_rec,mean_low_lig / count_low_lig,mean_high_rec/ count_high_rec,mean_high_lig / count_high_lig

import math

import numpy
import pandas as pd

def modes( benchmark_names, base_path, protein_list, protein_types):
    count_low_rec = 0
    count_high_rec = 0
    count_low_lig = 0
    count_high_lig = 0
    aa = ['LYS', 'PRO','ILE','TRP','GLU','GLN','GLY','SER','PHE','HIS','TYR','LEU','ASP','ASN','ARG','THR','ALA','CYS','VAL','MET']
    sec = ['C', 'E', 'B', 'T', 'H', 'G', 'b']

    Y = []
    X = []
    for proteinType, bm in zip(protein_types,benchmark_names):
            for protein in protein_list:
                    
                    rec_modes = json.load(open(base_path +'/'+ protein+ "/input/modes/{}A-{}-mn20-mt_hin99-eval.json".format( protein, proteinType),'r'))['modes']
                    lig_modes =json.load(open(base_path +'/'+ protein+ "/input/modes/{}B-{}-mn20-mt_hin99-eval.json".format( protein, proteinType),'r'))['modes']
                    rec_prot = json.load(open(base_path +'/'+ protein+ "/input/pdb/{}A-{}-evaluation.json".format( protein, proteinType),'r'))
                    lig_prot =json.load(open(base_path +'/'+ protein+ "/input/pdb/{}B-{}-evaluation.json".format( protein, proteinType),'r'))
                    for key, val in rec_modes.items():
                        Y.append(val['overlap'])
                        dat =  [val['eigenvalue']]
                        dat_d ={}
                        dat_d['y'] = val['overlap']
                        dat_d['eigenvalue'] = val['eigenvalue']
                        dat_d.update(val['area_aa'])
                        dat_d.update(val['area_sec'])

                        dat_d.update(rec_prot['aa_area'])
                        dat_d.update(rec_prot['sec_area'])
                        dat_d.update(rec_prot['aminoAcids'])
                        dat_d.update(rec_prot['secondary'])

                        # for key in aa:
                        #     dat.append(val['area_aa'][key])
                        # for key in sec:
                        #     dat.append(val['area_sec'][key])


                        # for key in aa:
                        #     dat.append(rec_prot['aa_area'][key])
                        # for key in sec:
                        #     dat.append(rec_prot['sec_area'][key])
                        # for key in aa:
                        #     dat.append(rec_prot['aminoAcids'][key])
                        # for key in sec:
                        #     dat.append(rec_prot['secondary'][key])
                        #for a in dat:
                        # for i in range(len(dat)):
                        #     k = dat[i]
                        #     if k == math.nan or k == math.inf or k == np.inf or k== np.nan:
                        #         dat[i] = 0    
                        X.append(dat_d)
                        

                    for key, val in lig_modes.items():
                        # Y.append(val['overlap'])
                        # dat =  [val['eigenvalue']]
                        # for key in aa:
                        #     dat.append(val['area_aa'][key])
                        # for key in sec:
                        #     dat.append(val['area_sec'][key])
                        # for key in aa:
                        #     dat.append(lig_prot['aa_area'][key])
                        # for key in sec:
                        #     dat.append(lig_prot['sec_area'][key])
                        # for key in aa:
                        #     dat.append(lig_prot['aminoAcids'][key])
                        # for key in sec:
                        #     dat.append(lig_prot['secondary'][key])
                        # for i in range(len(dat)):
                        #     k = dat[i]
                        #     if k == math.nan or k == math.inf or k == np.inf or k== np.nan:
                        #         dat[i] = 0
                        
                        dat_d = {}
                        dat_d['y'] = val['overlap']
                        dat_d['eigenvalue'] = val['eigenvalue']
                        dat_d.update(val['area_aa'])
                        dat_d.update(val['area_sec'])

                        dat_d.update(lig_prot['aa_area'])
                        dat_d.update(lig_prot['sec_area'])
                        dat_d.update(lig_prot['aminoAcids'])
                        dat_d.update(lig_prot['secondary'])
                        X.append(dat_d)
                        #X.append(np.asarray(dat))
    #X = np.asarray(X) 
    df =pd.DataFrame(X)
    #Y = np.asarray(Y) 
    return X,Y



def isConnected(base_path, protein_list):

    distance_receptor = 0
    distance_ligand = 0
    distance_receptor_list = []
    distance_ligand_list = []
    count_receptor = 0
    count_ligand = 0
    max_rec = 0
    max_lig = 0
    for protein in protein_list:
        receptor = os.path.join(base_path , protein) + "/{}A-unbound.pdb".format(protein)
        ligand = os.path.join(base_path , protein) + "/{}B-unbound.pdb".format(protein)
        receptor_capos = utils.getCoordinatesFromPDBlines(utils.getCAOnlyFromPDBLines(utils.readFileToList(receptor)))
        ligand_capos = utils.getCoordinatesFromPDBlines(utils.getCAOnlyFromPDBLines(utils.readFileToList(ligand)))
        rec_list = []
        lig_list = []
        for i in range(len(receptor_capos)-1):
            dpos = receptor_capos[i] - receptor_capos[i+1]
            distance_receptor +=np.sqrt( dpos.dot(dpos))
            #distance_receptor_list.append(np.sqrt( dpos.dot(dpos))) 
            rec_list.append(np.sqrt( dpos.dot(dpos))) 
            count_receptor += 1
        for i in range(len(ligand_capos)-1):
            dpos = ligand_capos[i] - ligand_capos[i+1]
            distance_ligand += np.sqrt(dpos.dot(dpos))
            #distance_ligand_list.append(np.sqrt( dpos.dot(dpos))) 
            lig_list.append(np.sqrt( dpos.dot(dpos))) 
            count_ligand += 1
        distance_receptor_list.append(np.asarray(rec_list))
        distance_ligand_list.append(np.asarray(lig_list))

    return distance_receptor/ count_receptor, distance_ligand / count_ligand, np.asarray(distance_receptor_list),np.asarray(distance_ligand_list)



def sortbyOverlap(base_path, protein_list):
    receptor_overlap = []
    ligand_overlap = []
    receptor_dict = {}
    ligand_dict = {}
    for protein in protein_list:
        receptor = os.path.join(base_path , protein) + "/input/modes/{}A-unbound-mn20-mt_hin99-eval.json".format(protein)
        ligand = os.path.join(base_path , protein) + "/input/modes/{}B-unbound-mn20-mt_hin99-eval.json".format(protein)
        receptor_json = json.load(open(receptor))
        ligand_json = json.load(open(ligand))
        ligand_overlap.append(ligand_json['modes']['1']['overlap'])
        receptor_overlap.append(receptor_json['modes']['1']['overlap'])
        ligand_dict[protein] = ligand_json['modes']['1']['overlap']
        receptor_dict[protein] = receptor_json['modes']['1']['overlap']
    receptor_overlap = np.asarray(receptor_overlap)
    ligand_overlap = np.asarray(ligand_overlap)
    protein_list = np.asarray(protein_list)
    rec_sorted = np.argsort(receptor_overlap)
    lig_sorted = np.argsort(ligand_overlap)
    return rec_sorted, lig_sorted, receptor_overlap[rec_sorted], ligand_overlap[lig_sorted], protein_list[lig_sorted], protein_list[rec_sorted], receptor_dict, ligand_dict


#evaluateBenchmarks(base_path, protein_list,[ 'unbound'], ['bm_dG_mr1_ml0_s1p000000_sO_c50_mr1_ml0_s1p000000_hin99'])  
bms = ['bm_dG_mr0_ml0_s1p000000_sO_c50_mr0_ml0_s1p000000_cut',
'bm_dG_mr0_ml0_s1p000000_sO_c50_mr0_ml0_s1p000000_hin99',
'bm_dG_mr0_ml0_s1p000000_sO_c50_mr0_ml0_s1p000000_hin99_singleDof',
'bm_dG_mr0_ml1_s1p000000_sO_c50_mr0_ml1_s1p000000_cut',
'bm_dG_mr0_ml1_s1p000000_sO_c50_mr0_ml1_s1p000000_hin99_singleDof',
'bm_dG_mr1_ml0_s1p000000_sO_c50_mr1_ml0_s1p000000_cut',
'bm_dG_mr1_ml0_s1p000000_sO_c50_mr1_ml0_s1p000000_hin99',
'bm_dG_mr1_ml0_s1p000000_sO_c50_mr1_ml0_s1p000000_hin99_singleDof',
'bm_dG_mr1_ml1_s0p100000_sO_c50_mr1_ml1_s0p100000_hin99_singleDof',
'bm_dG_mr1_ml1_s1p000000_sO_c50_mr1_ml1_s1p000000_cut']

types = [
"unbound-cut",
"unbound",
"unbound",
"unbound-cut",
"unbound",
"unbound-cut",
"unbound",
"unbound",
"unbound",
"unbound-cut"]

import prody as dy

def calcPCA(filenameEnsemble):
    ubi = dy.parsePDB(filenameEnsemble, subset='ca')
    ensemble = dy.Ensemble('ensemble')
    ensemble.setCoords( ubi.getCoords() )
    ensemble.addCoordset( ubi.getCoordsets() )
    ensemble.iterpose()
    pca = dy.PCA('Ubiquitin')
    pca.buildCovariance(ensemble)
    pca.calcModes()
    return pca


from subprocess import Popen, PIPE

# p = Popen("/path/to/env.sh", stdout=PIPE, stdin=PIPE, stderr=STDOUT) 

# p = Popen(["/home/glenn/Documents/Masterarbeit/concoord/bin/dist","-p",receptor,"-op",dist_rec], stdin=PIPE) #, shell=True
# p.communicate(input=b'1\n1\n')
# #p.communicate(input=b"1")
# #p.communicate(input=b"1")
# p.stdin.write(b'1\n')
# p.stdin.close()
# p.wait()
# p.stdin.write(b'1\n')
# p.stdin.close()
# p.wait()
def createPCAMOdes(base_path, protein_list):
    for protein in protein_list:
        receptor = os.path.join(base_path , protein) + "/{}A-unbound.pdb".format(protein)
        ligand = os.path.join(base_path , protein) + "/{}B-unbound.pdb".format(protein)
        pca_rec_folder = "{}/{}/input/pca/concoord/receptor".format(base_path, protein)
        pca_lig_folder = "{}/{}/input/pca/concoord/ligand".format(base_path, protein)


        dist_rec = "{}/{}A-dist".format(pca_rec_folder, protein)
        dist_lig = "{}/{}B-dist".format(pca_lig_folder, protein)

        disco_rec = "{}/{}A-disco.pdb".format(pca_rec_folder, protein)
        disco_lig = "{}/{}B-disco.pdb".format(pca_lig_folder, protein)

        nmdfile_rec = "{}/{}A-nmd".format(pca_rec_folder, protein)
        nmdfile_lig = "{}/{}B-nmd".format(pca_lig_folder, protein)

        os.system("mkdir -p {}".format(pca_rec_folder))
        os.system("mkdir -p {}".format(pca_lig_folder))
        #pwd = os.getcwd()

        os.chdir(pca_rec_folder)
        p = Popen(["/home/glenn/Documents/Masterarbeit/concoord/bin/dist","-p",receptor], stdin=PIPE) #, shell=True #,"-op",dist_rec
        p.communicate(input=b'1\n1\n')        
        os.system("/home/glenn/Documents/Masterarbeit/concoord/bin/disco -on {} -n 200 -i 1000 -viol 1. -bump ".format( disco_rec ))

        os.chdir(pca_lig_folder)
        p = Popen(["/home/glenn/Documents/Masterarbeit/concoord/bin/dist","-p",ligand ], stdin=PIPE) #, shell=True
        p.communicate(input=b'1\n1\n')
        os.system("/home/glenn/Documents/Masterarbeit/concoord/bin/disco -on {} -n 200 -i 1000 -viol 1. -bump  ".format( disco_lig ))

        try:
            pca_rec = calcPCA(disco_rec)
            atoms_rec = dy.parsePDB(receptor, subset='ca')
            dy.writeNMD(nmdfile_rec,pca_rec,atoms_rec)
        except:
            pass


        try:
            pca_lig = calcPCA(disco_lig)
            atoms_lig = dy.parsePDB(ligand, subset='ca')
            dy.writeNMD(nmdfile_lig,pca_lig,atoms_lig)
        except:
            pass


                    # os.system("rm {}/BONDS.DAT".format(pwd))
        # os.system("rm {}/MARGINS.DAT".format(pwd))
        # os.system("rm {}/ATOMS.DAT".format(pwd))

        #os.system("cd {}".format(pca_lig_folder))
        #os.system("/home/glenn/Documents/Masterarbeit/concoord/bin/dist -p {} -op {} ".format(ligand,dist_lig))



def transposeFromCAToFull(mode_dict, residueNamesFull):
    index = 0
    
    new_dict ={}
    size = len(residueNamesFull)
    for key, val in mode_dict.items():
        new_dict[key] = {'eval': val['eval'], 'evec':np.zeros((size,3))}
    currRes = residueNamesFull[0]
    for i,res in enumerate(residueNamesFull):
        for key in new_dict.keys():
            new_dict[key]['evec'][i] = mode_dict[key]['evec'][index] 
        if currRes != res:
            currRes = res
            index += 1
    return new_dict


def brownianTonormalModes(base_path, prot_list, mode_extension = "-unbound-cut-cut_modes.dat", pdb_extension  = "-unbound-cut-r.pdb", out_modes_extension = "-unbound-cut-modes-20-brow.dat"):
    for prot in prot_list:
        rec_mode_name = "{}/{}/input/webnma/{}A{}".format(base_path, prot, prot,mode_extension)
        lig_mode_name = "{}/{}/input/webnma/{}B{}".format(base_path, prot, prot,mode_extension)

        rec_mode_full_name = "{}/{}/input/modes/{}A{}".format(base_path, prot, prot,out_modes_extension)
        lig_mode_full_name = "{}/{}/input/modes/{}B{}".format(base_path, prot, prot,out_modes_extension)

        rec_protein_name = "{}/{}/input/pdb/{}A{}".format(base_path, prot, prot,pdb_extension)
        lig_protein_name = "{}/{}/input/pdb/{}B{}".format(base_path, prot, prot,pdb_extension)

        rec_modes = utils.readHinsenModeFile(rec_mode_name)
        lig_modes = utils.readHinsenModeFile(lig_mode_name)

        rec_mode_full = transposeFromCAToFull(rec_modes, utils.getResidueFromPDBlines(utils.readFileToList(rec_protein_name)))
        lig_mode_full = transposeFromCAToFull(lig_modes, utils.getResidueFromPDBlines(utils.readFileToList(lig_protein_name)))
        
        rec_truncated_20 ={}
        lig_truncated_20 ={}
        for key in range(1,21):
            rec_truncated_20[key] = rec_mode_full[key]
            lig_truncated_20[key] = lig_mode_full[key]

        utils.writeModeFileFromDict(rec_truncated_20, rec_mode_full_name)
        utils.writeModeFileFromDict(lig_truncated_20, lig_mode_full_name)


def createModeFileFromNMD(base_path, prot_list):
    pcaType = "modeller"
    pcaType = "concoord"
    modeType = "pca_concoord"
    for prot in prot_list:
       # try:
            modeFileName_rec = "modes_{}_A_200_strt0_noT.nmd".format(prot)
            modeFileName_lig = "modes_{}_B_200_strt0_noT.nmd".format(prot)
            modeFileName_rec_concoord = "{}A-nmd.nmd".format(prot)
            modeFileName_lig_concoord = "{}B-nmd.nmd".format(prot)

            pdb_rec_lines = utils.readFileToList( "{}/{}/input/pdb/{}A-unbound-r.pdb".format(base_path, prot, prot))
            pdb_lig_lines = utils.readFileToList("{}/{}/input/pdb/{}B-unbound-r.pdb".format(base_path, prot, prot))

            nmd_rec_file = "{}/{}/input/pca/{}/receptor/{}".format(base_path, prot, pcaType, modeFileName_rec_concoord)
            nmd_lig_file = "{}/{}/input/pca/{}/ligand/{}".format(base_path, prot, pcaType, modeFileName_lig_concoord)

            modes_rec = "{}/{}/input/modes/{}A-unbound-modes-20-{}.dat".format(base_path, prot, prot , modeType )
            modes_lig = "{}/{}/input/modes/{}B-unbound-modes-20-{}.dat".format(base_path, prot, prot , modeType )


            nmd_rec = utils.parseNMDTomodeDict(nmd_rec_file)
            nmd_lig = utils.parseNMDTomodeDict(nmd_lig_file)

            norm_rec = np.sqrt(len(nmd_rec[1]['evec']))
            norm_lig = np.sqrt(len(nmd_lig[1]['evec']))

            print(norm_rec, norm_lig)
            modes_rec_transpose = transposeFromCAToFull(nmd_rec,utils.getResidueFromPDBlines(pdb_rec_lines))
            modes_lig_transpose = transposeFromCAToFull(nmd_lig,utils.getResidueFromPDBlines(pdb_lig_lines))

            for key in modes_rec_transpose.keys():
                print(modes_rec_transpose[key]['eval']  )
               # modes_rec_transpose[key]['eval'] = modes_rec_transpose[key]['eval']/ norm_rec
            # for key, val in modes_lig_transpose.items():
            #     print(modes_lig_transpose[key]['eval'] )

            #     modes_lig_transpose[key]['eval'] /= norm_lig

            utils.writeModeFileFromDict(modes_rec_transpose,modes_rec)
            utils.writeModeFileFromDict(modes_lig_transpose,modes_lig)

#createModeFileFromNMD(base_path, prot_list)

        # except:
        #     pass


def getbms(base_path):
    return [os.path.basename(os.path.normpath(d))  for d in glob.glob(base_path + '/1ACB/**/') if 'bm' in d]




df = pd.read_csv('dataframe_2842019') 
df = dfmul
df['sum'] = 0
df['sum'] = df['num1_50'] + df['num2_50'] + df['num3_50']
keys = df[df['sum'] > 0 ].groupby('bm').size().sort_values()
print(keys)
evalcols = ['num1_50', 'num2_50', 'num3_50'] 
#evalcols = ['sum'] 

plt.clf() 
colors = ['b','r','g'] 
res_dict = {} 
df_group = df.groupby('bm')

xticks = []
for i,name in enumerate(keys.index): 
    d = df_group.get_group(name)
    lastval = 0 
    xticks.append(name)
    for idx,col in enumerate(evalcols):  
        d2 = d[d[col] > 0]  
        plt.bar(i,d2.shape[0], bottom = lastval, color = colors[idx]) 
        lastval += d2.shape[0] 
    res_dict[name] = lastval 

plt.xticks(np.arange(0, len(xticks)),xticks,rotation=90)
plt.show()     



#plotsingle
df = dfsin
plt.clf() 
colors = ['b','r','g'] 
res_dict = {} 
df_group = df.groupby('bm')
keys = df.groupby('bm')['lrmsd'].mean().sort_values()
print(keys)

xticks = []
for i,name in enumerate(keys.index): 
    d = df_group.get_group(name)
    xticks.append(name)   
    plt.bar(i,d['lrmsd'].mean(),  color = colors[0]) 
    res_dict[name] = lastval 

plt.xticks(np.arange(0, len(xticks)),xticks,rotation=90)
plt.show()     

#print benchmarks
print("ModesRec ModesLig scale    modetype     bound    cut pruned rigidStart manipulate")
for name,d in dfmul.groupby('bm'):
    row = d.iloc[0]
    print("{:2d}      {:2d}        {:4.3f}    {:10s} {:6} {:6} {:6} {:6} {:6} ".format(  row.numModesRec, row.numModesLig, row.scale, row.modetype, row.bound, row.cut, 'pruned' in row.bm, 'rigidStart' in row.bm, 'manipulate' in row.bm))


benchmark_names = ["bm_dG_mr1_ml1_s1p000000_sO_c50_mr1_ml1_s1p000000_hin99",  "bm_dG_mr0_ml0_s1p000000_sO_c50_mr0_ml0_s1p000000_hin99"]
for  bm in benchmark_names:
    for protein in prot_list:
        print("processing protein", protein ,"in benchmark", bm)
        path = getBenchmarkpath(base_path, protein, bm)
        rmsd_filename =glob.glob( os.path.join(path,"analysis/*-rmsd.dat"))[0]
        energy_filename = glob.glob(os.path.join(path,"analysis/*-deredundant.dat"))[0]




#plot best or worst

protein_list = worstList
benchmark_names = ["bm_dG_mr1_ml1_s1p000000_sO_c50_mr1_ml1_s1p000000_hin99",  "bm_dG_mr0_ml0_s1p000000_sO_c50_mr0_ml0_s1p000000_hin99"]
for protein in protein_list:
    plt.clf()
    for  bm in benchmark_names:
        print("processing protein", protein ,"in benchmark", bm)
        path = getBenchmarkpath(base_path, protein, bm)
        rmsd_filename =glob.glob( os.path.join(path,"analysis/*-rmsd.dat"))[0]
        energy_filename = glob.glob(os.path.join(path,"analysis/*-deredundant.dat"))[0]
        energies = readEnergies(energy_filename)
        rmsd = readRMSD(rmsd_filename)
        plt.scatter(rmsd, energies,label = bm)
    plt.xlabel("RMSD (A)")
    plt.title(protein)
    plt.legend()
    plt.ylabel("Energy (kcal/mol)")
    #plt.show()
    #plt.savefig("/home/glenn/Documents/Masterarbeit/evaluation/rmsdPlot_{}_best.png".format(protein))
    fig = plt.gcf()
    fig.set_figheight(7)
    fig.set_figwidth(8)
    plt.savefig("/home/glenn/Documents/Masterarbeit/evaluation/rmsdPlot_{}_worst.png".format(protein))
    #break


#overlapplot
# rec_sorted, lig_sorted, receptor_overlap, ligand_overlap, protein_list, protein_list, receptor_dict, ligand_dict = sortbyOverlap(base_path, prot_list)
    
# sorted_diff = diff.sort_values().keys().tolist()
# lig_ol = []
# rec_ol = []
# for prot in sorted_diff:
#     lig_ol.append(ligand_dict[prot])
#     rec_ol.append(receptor_dict[prot])

# plt.plot( savgol_filter(lig_ol, 3,1), label = "ligand")
# plt.plot( savgol_filter(rec_ol, 3,1), label = "receptor")
# plt.xticks(np.arange(0, len(sorted_diff)),sorted_diff, rotation = 90)
# plt.legend()
# plt.ylabel('overlap')
# plt.xlabel('protein')
# plt.show()



#copy best collect files to evaluation 
best = diff.sort_values().nlargest(10).keys().tolist()
benchmark_modes = "bm_dG_mr1_ml1_s0p100000_sO_c50_mr1_ml1_s0p100000_hin99"
for b in best:
    collectfile = "{}/{}/{}/analysis/{}-unbound-collect_high.pdb".format(base_path, b, benchmark_modes, b)
    print(collectfile)
    os.system("cp {} /home/glenn/Documents/Masterarbeit/evaluation/best_collect_high/".format(collectfile))
    

best = diff.sort_values().nlargest(10).keys().tolist()
benchmark_modes = "bm_dG_mr1_ml1_s0p100000_sO_c50_mr1_ml1_s0p100000_hin99"
for b in best:
    dofeval = "{}/{}/{}/analysis/{}-unbound-dof_eval_high.json".format(base_path, b, benchmark_modes, b)
    evaluation = json.load(open(dofeval, 'r'))
    










#plot sorted diff sum 50 
best = diff.sort_values().nlargest(10).keys().tolist()
benchmark_modes = "bm_dG_mr1_ml1_s1p000000_sO_c50_mr1_ml1_s1p000000_hin99"
dof_rec = []
dof_lig = []
best = diff.sort_values().keys().tolist()
for b in best:
    dofeval = "{}/{}/{}/analysis/{}-unbound-dof_eval_low.json".format(base_path, b, benchmark_modes, b)
    evaluation = json.load(open(dofeval, 'r'))
    for key, val in evaluation.items():
        print(b,key, val['rec']['1']['dof'],val['rec']['1']['mode'], val['lig']['1']['dof'],val['lig']['1']['mode'], val['rec']['1']['dof']/val['rec']['1']['mode'], val['lig']['1']['dof']/val['lig']['1']['mode'])
    dof_rec.append(val['rec']['1']['dof']/val['rec']['1']['mode'])
    dof_lig.append(val['lig']['1']['dof']/val['lig']['1']['mode'])
plt.plot(dof_rec, label = 'rec')
plt.plot(dof_lig, label = 'lig')
plt.plot(np.ones(len(dof_lig)))
plt.legend()
plt.xticks(np.arange(0, len(dof_lig)), best, rotation = 90)
plt.xlabel('protein')
plt.ylabel("mode dof result / mode dof optimal")
plt.ylim([-10,10])
plt.show()