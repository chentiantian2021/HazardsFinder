# -*- coding: utf-8 -*-
# Created by: ctt@2022
# chentiantian@dicp.ac.cn



from threading import Thread
import pandas as pd
import os
import numpy as np
from pyteomics import mzml
import get_area
import get_isotopic
import get_sim
from collections import defaultdict
from datetime import datetime

def generate_standard_file():
    file_name = pd.DataFrame()
    file_name['Sample_Name'] = ''
    file_name['Name'] = ''
    file_name['PubChemCID'] = ''
    file_name['CAS'] = ''
    file_name['Formula'] = ''
    file_name['SMILES'] = ''
    file_name['InChIKey'] = ''
    file_name['m/ztheor'] = ''
    file_name['m/zexper'] = ''
    file_name['m/z-inten'] = ''
    file_name['tR/mintheor'] = ''
    file_name['tR/min-MS'] = ''
    file_name['tR/min-MSMS'] = ''
    file_name['EIC'] = ''
    file_name['SIM'] = ''
    file_name['MS/MStheor'] = ''
    file_name['MS/MSexper'] = ''
    file_name['mass error matching'] = ''
    file_name['isotopic deviation matching'] = ''
    file_name['tR deviation matching'] = ''
    file_name['MS/MS matched matching'] = ''
    file_name['Confidence level'] = ''
    return file_name

def suspect1_sceening(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc):
    output_standard=generate_standard_file()
    row=-1

    for mzML_file in file_list:
        print(mzML_file)
        print('  ')
        os.chdir(sample_dir)
        file = mzml.read(mzML_file)
        for spe in file:
            if spe['ms level'] == 2:
                precursor = spe['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                spectral_rt = spe['scanList']['scan'][0]['scan start time']
                values_to_find = [round(precursor) - 1, round(precursor), round(precursor) + 1]
                standard_indices = [index for index, value in enumerate(standard_mz_round) if value in values_to_find]
                for i in standard_indices:
                    if 10 ** 6 * abs(precursor - standard_precursorMass[i]) / precursor < ms1_error and abs(
                            spectral_rt - standard_rt[i]) < rt_error:
                        spectral = spe['m/z array']
                        spectral_inten = spe['intensity array'] / max(spe['intensity array'])
                        msms_standard = standard_msms[i]
                        similarity, msms_exper = get_sim.sim_cal_2(precursor, spectral, spectral_inten, msms_standard)
                        if similarity > SIM_sc:
                            file_raw = mzML_file.split('.')[0]
                            row = row + 1
                            output_standard.loc[row, 'Sample_Name'] = file_raw
                            output_standard.loc[row, 'Name'] = standard_name[i]
                            output_standard.loc[row, 'PubChemCID'] = standard_cid[i]
                            output_standard.loc[row, 'CAS'] = standard_cas[i]
                            output_standard.loc[row, 'SMILES'] = standard_smiles[i]
                            output_standard.loc[row, 'Formula'] = standard_formula[i]
                            output_standard.loc[row, 'InChIKey'] = standard_InchiKey[i]
                            output_standard.loc[row, 'm/ztheor'] = round(standard_precursorMass[i], 4)
                            output_standard.loc[row, 'm/zexper'] = round(precursor, 4)
                            output_standard.loc[row, 'tR/min-MSMS'] = round(spectral_rt, 2)
                            output_standard.loc[row, 'tR/mintheor'] = round(standard_rt[i], 2)
                            output_standard.loc[row, 'SIM'] = round(similarity, 2)
                            output_standard.loc[row, 'MS/MStheor'] = msms_standard
                            output_standard.loc[row, 'MS/MSexper'] = msms_exper
                            output_standard.loc[row, 'mass error matching'] = 'Success'
                            output_standard.loc[row, 'MS/MS matched matching'] = 'Success'
                            output_standard.loc[row, 'tR deviation matching'] = 'Success'
                            output_standard.loc[row, 'Confidence level'] = "CL1"

    return output_standard

def suspect1_threading1(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc):
    global output_standard1
    output_standard1=suspect1_sceening(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc)

def suspect1_threading2(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc):
    global output_standard2
    output_standard2=suspect1_sceening(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc)

def suspect1_threading3(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc):
    global output_standard3
    output_standard3=suspect1_sceening(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc)

def suspect1_threading4(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc):
    global output_standard4
    output_standard4=suspect1_sceening(file_list,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc)
def check_standardisotopes_etc(output_standard,sample_dir,current_path,isotope_sc,rt_error,ms1_error):
    print('......正在进行同位素等核查......')
    output_standard = output_standard.sort_values(by=['PubChemCID','SIM'], ascending=[False, False])
    output_standard = output_standard.reset_index(drop=True)
    pubchemcid = list( output_standard['PubChemCID'])
    pubchemcid_statistics = defaultdict(list)
    for n, v in enumerate(pubchemcid):
        pubchemcid_statistics[v].append(n)
    for e in pubchemcid_statistics:
        pubchemcid_eve_num = len(pubchemcid_statistics[e])
        for n in range(pubchemcid_eve_num):
            index = pubchemcid_statistics[e][n]
            rt =  output_standard.loc[index, 'tR/min-MSMS']
            mz =  output_standard.loc[index, 'm/zexper']
            mzML_file =  output_standard.loc[index, 'Sample_Name'] + '.mzML'
            os.chdir(sample_dir)
            max_height, max_height_time, scan_times, chromatogram = get_area.get_area(mz, rt, mzML_file)
            output_standard.loc[index, 'm/z-inten'] = round(max_height)
            output_standard.loc[index, 'tR/min-MS'] = round(max_height_time, 2)
            eic = ''
            for x in range(len(scan_times)):
                eic = eic + str(round(scan_times[x], 3)) + ' ' + str(round(chromatogram[x], 1)) + ';'
            eic = eic.strip(';')
            output_standard.loc[index, 'EIC'] = eic

            formula_eve =  output_standard.loc[index, 'Formula']
            ratio_ther, M_1, M_2 = get_isotopic.isotopic(formula_eve, mz)
            os.chdir(sample_dir)
            file = mzml.read(mzML_file)
            isotopic_num = 0
            for spe in file:
                p_1_exp = 0
                p_2_exp = 0
                if spe['ms level'] == 1:
                    mzml_rt = float(spe['scanList']['scan'][0]['scan start time'])
                    if mzml_rt - rt < rt_error:
                        mzml_ms = spe['m/z array']
                        mzml_ms_inten = spe['intensity array']
                        keep = np.where(mzml_ms_inten > 0)
                        mzml_ms = mzml_ms[keep]
                        mzml_ms_inten = mzml_ms_inten[keep]
                        if abs(mzml_rt - rt) < rt_error:
                            for g in range(len(mzml_ms)):
                                if (10 ** 6) * abs(M_1 - mzml_ms[g]) / mzml_ms[g] < 2*ms1_error:
                                    p_1_exp = mzml_ms_inten[g]
                                if (10 ** 6) * abs(M_2 - mzml_ms[g]) / mzml_ms[g] < 2*ms1_error:
                                    p_2_exp = mzml_ms_inten[g]
                    if mzml_rt - rt > rt_error:
                        break
                if p_1_exp > 0 and p_2_exp > 0:
                    ratio_exp = p_1_exp/p_2_exp
                    if abs(ratio_ther - ratio_exp) / ratio_exp <= isotope_sc:
                        isotopic_num = isotopic_num + 1
                        if isotopic_num >= 1:
                            break
            if isotopic_num >= 1:
                output_standard.loc[index, 'isotopic deviation matching'] = 'Success'
                os.chdir(sample_dir)
                break

    retain_list = []
    for i in range(output_standard.shape[0]):
        if output_standard.loc[i, 'isotopic deviation matching'] == 'Success':
            retain_list.append(i)
    output_standard = output_standard.loc[retain_list, :]
    output_standard.reset_index(drop=True, inplace=True)
    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name=str(current_date)+'_standard.xlsx'
    output_standard.to_excel(file_name,index=False)

def screening(sample_dir,current_path,standard_database,ms1_error,rt_error,isotope_sc,SIM_sc):
    standard_database = standard_database.sort_values(by='m/z')
    standard_database = standard_database.reset_index(drop=True)
    standard_precursorMass = standard_database['m/z']
    standard_mz_round = [round(number) for number in list(standard_precursorMass)]
    standard_rt=standard_database['tR/min']
    standard_name=standard_database['Name']
    standard_cid=standard_database['PubChemCID']
    standard_cas=standard_database['CAS']
    standard_formula=standard_database['Formula']
    standard_smiles=standard_database['SMILES']
    standard_InchiKey=standard_database['InChIKey']
    standard_msms=standard_database['MS/MS']
    files = [f for f in os.listdir(sample_dir) if f.endswith('.mzML')]
    file_list1 = []
    file_list2 = []
    file_list3 = []
    file_list4 = []
    for x in range(0, len(files), 4):
        file_list1.append(files[x])
    for x in range(1, len(files), 4):
        file_list2.append(files[x])
    for x in range(2, len(files), 4):
        file_list3.append(files[x])
    for x in range(3, len(files), 4):
        file_list4.append(files[x])
    for file in files:
        if file not in file_list1 and file not in file_list2 and file not in file_list3 and file not in file_list4:
            file_list4.append(file)

    t1 = Thread(target=suspect1_threading1, args=(file_list1,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc,))
    t2 = Thread(target=suspect1_threading2, args=(file_list2,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc,))
    t3 = Thread(target=suspect1_threading3, args=(file_list3,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc,))
    t4 = Thread(target=suspect1_threading4, args=(file_list4,sample_dir,standard_mz_round,standard_precursorMass,standard_rt,standard_msms,standard_name,standard_cid,standard_cas,standard_smiles,standard_formula,standard_InchiKey,ms1_error,rt_error,SIM_sc,))
    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t1.join()
    t2.join()
    t3.join()
    t4.join()
    output_standard = pd.concat([output_standard1, output_standard2, output_standard3, output_standard4], axis=0,
                                ignore_index=True)

    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name = str(current_date) + '_original_standard.xlsx'
    output_standard.to_excel(file_name, index=False)
    check = Thread(target=check_standardisotopes_etc, args=(output_standard,sample_dir,current_path,isotope_sc,rt_error,ms1_error,))
    check.start()
    check.join()





