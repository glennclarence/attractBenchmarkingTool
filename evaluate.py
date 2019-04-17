import utils as utils
import numpy as np
from collections import Counter
import os
import pandas as pd
import json

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


def evaluateBenchmarks( base_path, protein_list,proteinTypes, benchmark_names):
    result = []
    for proteinType, bm in zip(proteinTypes,benchmark_names):
        for protein in protein_list:
            print("processing protein", protein ,"in benchmark", bm)
            path = getBenchmarkpath(base_path, protein, bm)
            rmsd_filename = os.path.join(path,"analysis/{}-{}-rmsd.dat".format(protein, proteinType ))
            irmsd_filename = os.path.join(path,"analysis/{}-{}-irmsd.dat".format(protein, proteinType ))
            fnat_filename = os.path.join(path,"analysis/{}-{}-fnat.dat".format(protein, proteinType ))
            energy_filename = os.path.join(path,"analysis/{}-{}-deredundant.dat".format(protein, proteinType ))
            
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
    return pd.DataFrame(result)



def evaluateSingle( base_path, protein_list,proteinTypes, benchmark_names):
    result =[]
    for proteinType, bm in zip(proteinTypes,benchmark_names):
        for protein in protein_list:
            try:
                path = getBenchmarkpath(base_path, protein, bm)
                rmsd_filename = os.path.join(path,"analysis/{}-{}-rmsd.dat".format(protein, proteinType ))
                irmsd_filename = os.path.join(path,"analysis/{}-{}-irmsd.dat".format(protein, proteinType ))
                fnat_filename = os.path.join(path,"analysis/{}-{}-fnat.dat".format(protein, proteinType ))
                energy_filename = os.path.join(path,"analysis/{}-{}-deredundant.dat".format(protein, proteinType ))
                dof_filename = os.path.join(path,"analysis/{}-{}-dof_eval.json".format(protein, proteinType ))

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

                dof_eval = json.load(open(dof_filename))["1"]
                for modeIdx,value  in dof_eval['rec'].items():
                    row['rec_ratiomode_{}'.format(modeIdx)] = value['ratio']
                for modeIdx,value  in dof_eval['lig'].items():
                    row['lig_ratiomode_{}'.format(modeIdx)] = value['ratio']
                result.append(row)
            except:
                print("couldn't get protein", protein,bm)
                pass
    return pd.DataFrame(result)



def interface( benchmark_names, base_path, protein_list, protein_types):

    result = {'C': 0, 'E': 0, 'B': 0, 'T': 0, 'H': 0, 'G': 0, 'b': 0}
    result_high = {'C': 0, 'E': 0, 'B': 0, 'T': 0, 'H': 0, 'G': 0, 'b': 0}

    count = 0
    count_high = 0
    for proteinType, bm in zip(protein_types,benchmark_names):
            for protein in protein_list:
                try:
                    path_low = getBenchmarkpath(base_path, protein, bm) + "/analysis/{}-{}-interface_high.json".format( protein, proteinType)
                    path_high = getBenchmarkpath(base_path, protein, bm) + "/analysis/{}-{}-interface_low.json".format( protein, proteinType)
                    data_low = json.load(open(path_low))
                    data_high = json.load(open(path_high))
                    for i in data_low['interfaces']:
                        count += 1
                        for key, val in i['countSecRec'].items():
                            result[key] += val
                    for i in data_high['interfaces']:
                        count_high += 1
                        for key, val in i['countSecRec'].items():
                            result_high[key] += val
                except:
                    print('unable to read', protein)
                    pass
    for key, val in result.items():
        result[key] /= float(count)
    for key, val in result_high.items():
        result_high[key] /= float(count_high)
    return result, result_high
#evaluateBenchmarks(base_path, protein_list,[ 'unbound'], ['bm_dG_mr1_ml0_s1p000000_sO_c50_mr1_ml0_s1p000000_hin99'])  
bms = [
'bm_dG_mr0_ml0_s1p000000_sO_c50_mr0_ml0_s1p000000_cut',
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




