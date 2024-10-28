# -*- coding: utf-8 -*-
# Created by: ctt@2022
# chentiantian@dicp.ac.cn


from threading import Thread
import pandas as pd
import os
import numpy as np
from pyteomics import mzml
import get_area
import get_sim
import get_isotopic
from datetime import datetime
from collections import defaultdict

def remove_extra(result_df):
    remove_list=[]
    for i in range(result_df.shape[0]):
        for g in range(i,result_df.shape[0]):
            if i != g:
                if result_df.loc[i, 'PubChemCID'] == result_df.loc[g, 'PubChemCID'] and abs(result_df.loc[i, 'tR/min-MSMS'] -result_df.loc[g, 'tR/min-MSMS']) <= 0.21:
                    msms_1 = result_df.loc[i, 'MS/MSexper']
                    msms_2 = result_df.loc[g, 'MS/MSexper']
                    similarity = get_sim.sim_cal_1(msms_1, msms_2)
                    if similarity > 0.6:
                        remove_list.append(g)
    retain_list=[]
    for i in range(result_df.shape[0]):
        if i not in remove_list:
            retain_list.append(i)
    result_df=result_df.loc[retain_list,:]
    result_df.reset_index(drop=True, inplace=True)
    return result_df



def remove_extra2(result_df):
    sample_list = list(result_df['Sample_Name'])
    sample_list = list(set(sample_list))
    new_df = result_df[['Sample_Name', 'PubChemCID','tR/min-MSMS']]
    specific_new_df = new_df.drop_duplicates()
    specific_new_df.reset_index(drop=True, inplace=True)

    remove_pubchemcid = []
    for i in range(len(sample_list)):
        pubchemcid=[]
        for l in range(specific_new_df.shape[0]):
            if sample_list[i]==specific_new_df.loc[l,'Sample_Name']:
                pubchemcid.append(specific_new_df.loc[l,'PubChemCID'])
        pubchemcid_statistics = defaultdict(list)
        for n, v in enumerate(pubchemcid):
            pubchemcid_statistics[v].append(n)
        for e in pubchemcid_statistics:
            if len(pubchemcid_statistics[e])>10:
                remove_pubchemcid.append(e)

    remove_list=[]
    for e in remove_pubchemcid:
        for i in range(result_df.shape[0]):
            if str(e) ==str(result_df.loc[i,'PubChemCID']):
                remove_list.append(i)

    retain_list=[]
    for i in range(result_df.shape[0]):
        if i not in remove_list:
            retain_list.append(i)
    result_df=result_df.loc[retain_list,:]
    result_df.reset_index(drop=True, inplace=True)
    return result_df

def generate_standard_online_file():
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


def generate_literature_file():
    output_literature = pd.DataFrame()
    output_literature['Sample_Name'] = ''
    output_literature['Name'] = ''
    output_literature['PubChemCID'] = ''
    output_literature['CAS'] = ''
    output_literature['Formula'] = ''
    output_literature['SMILES'] = ''
    output_literature['InChIKey'] = ''
    output_literature['m/ztheor'] = ''
    output_literature['m/zexper'] = ''
    output_literature['tR/mintheor'] = ''
    output_literature['m/z-inten'] = ''
    output_literature['tR/min-MS'] = ''
    output_literature['tR/min-MSMS'] = ''
    output_literature['EIC'] = ''
    output_literature['SIM-num'] = ''
    output_literature['SIM-Frg'] = ''
    output_literature['MS/MStheor'] = ''
    output_literature['MS/MSexper'] = ''
    output_literature['mass error matching'] = ''
    output_literature['isotopic deviation matching'] = ''
    output_literature['tR deviation matching'] = ''
    output_literature['MS/MS matched matching'] = ''
    output_literature['Confidence level'] = ''

    return output_literature

def suspect_sceening(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature):
    output_online=generate_standard_online_file()
    row_online = -1
    output_literature=generate_literature_file()
    row_literature=-1

    for mzML_file in file_list:
        print(mzML_file)
        print('  ')
        os.chdir(sample_dir)
        file = mzml.read(mzML_file)
        for spe in file:
            if spe['ms level'] == 2:
                precursor = spe['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                spectral_rt = spe['scanList']['scan'][0]['scan start time']
                values_to_find = [round(precursor)-1,round(precursor),round(precursor)+1]
                online_indices = [index for index, value in enumerate(online_mz_round) if value in values_to_find]
                for i in online_indices:
                    if 10 ** 6 * abs(precursor - online_precursorMass[i]) / precursor < ms1_error:
                        spectral = spe['m/z array']
                        spectral_inten = spe['intensity array'] / max(spe['intensity array'])
                        msms_online = online_msms[i]
                        similarity, msms_exper = get_sim.sim_cal_2(precursor, spectral, spectral_inten, msms_online)
                        if similarity > SIM_online and abs(spectral_rt - online_rt[i]) < rt_error:
                            file_raw = mzML_file.split('.')[0]
                            row_online = row_online + 1
                            output_online.loc[row_online, 'Sample_Name'] = file_raw
                            output_online.loc[row_online, 'Name'] = online_name[i]
                            output_online.loc[row_online, 'PubChemCID'] = online_cid[i]
                            output_online.loc[row_online, 'CAS'] = online_cas[i]
                            output_online.loc[row_online, 'SMILES'] = online_smiles[i]
                            output_online.loc[row_online, 'Formula'] = online_formula[i]
                            output_online.loc[row_online, 'InChIKey'] = online_InchiKey[i]
                            output_online.loc[row_online, 'm/ztheor'] = round(online_precursorMass[i], 4)
                            output_online.loc[row_online, 'm/zexper'] = round(precursor, 4)
                            output_online.loc[row_online, 'tR/min-MSMS'] = round(spectral_rt, 2)
                            output_online.loc[row_online, 'tR/mintheor'] = round(online_rt[i], 2)
                            output_online.loc[row_online, 'SIM'] = round(similarity, 2)
                            output_online.loc[row_online, 'MS/MStheor'] = online_msms[i]
                            output_online.loc[row_online, 'MS/MSexper'] = msms_exper
                            output_online.loc[row_online, 'mass error matching'] = 'Success'
                            output_online.loc[row_online, 'MS/MS matched matching'] = 'Success'
                            output_online.loc[row_online, 'tR deviation matching'] = 'Success'
                            output_online.loc[row_online, 'Confidence level'] = "CL2a"
                        if similarity > SIM_online and abs(spectral_rt - online_rt[i]) >= rt_error:
                            file_raw = mzML_file.split('.')[0]
                            row_online = row_online + 1
                            output_online.loc[row_online, 'Sample_Name'] = file_raw
                            output_online.loc[row_online, 'Name'] = online_name[i]
                            output_online.loc[row_online, 'PubChemCID'] = online_cid[i]
                            output_online.loc[row_online, 'CAS'] = online_cas[i]
                            output_online.loc[row_online, 'SMILES'] = online_smiles[i]
                            output_online.loc[row_online, 'Formula'] = online_formula[i]
                            output_online.loc[row_online, 'InChIKey'] = online_InchiKey[i]
                            output_online.loc[row_online, 'm/ztheor'] = round(online_precursorMass[i], 4)
                            output_online.loc[row_online, 'm/zexper'] = round(precursor, 4)
                            output_online.loc[row_online, 'tR/min-MSMS'] = round(spectral_rt, 2)
                            output_online.loc[row_online, 'tR/mintheor'] = round(online_rt[i], 2)
                            output_online.loc[row_online, 'SIM'] = round(similarity, 2)
                            output_online.loc[row_online, 'MS/MStheor'] = online_msms[i]
                            output_online.loc[row_online, 'MS/MSexper'] = msms_exper
                            output_online.loc[row_online, 'mass error matching'] = 'Success'
                            output_online.loc[row_online, 'MS/MS matched matching'] = 'Success'
                            output_online.loc[row_online, 'Confidence level'] = "CL3"

                literature_indices = [index for index, value in enumerate(literature_mz_round) if value in values_to_find]
                for i in literature_indices:
                    if 10 ** 6 * abs(precursor - literature_precursorMass[i]) / precursor < ms1_error:
                        num = 0
                        spectral = spe['m/z array']
                        spectral_inten = spe['intensity array'] / max(spe['intensity array'])
                        spectral_inten = spectral_inten / max(spectral_inten)
                        keep1 = np.where(spectral_inten >= 0.01)
                        spectral = spectral[keep1]
                        spectral_inten=spectral_inten[keep1]
                        SIM_Frg = ''
                        msms_reference = literature_msms[i]
                        try:
                            msms_reference = msms_reference.strip(',')
                            msms_reference = msms_reference.split(',')
                            length = len(msms_reference)
                            for m in range(len(spectral)):
                                for n in range(len(msms_reference)):
                                    if 10 ** 6 * (
                                            abs(spectral[m] - float(msms_reference[n])) / float(msms_reference[n])) < 2*ms1_error:
                                        num = num + 1
                                        SIM_Frg = SIM_Frg + str(round(float(msms_reference[n]),4))+ ','
                        except:
                            length = 1
                            for m in range(len(spectral)):
                                if 10 ** 6 * (abs(spectral[m] - float(msms_reference)) / float(msms_reference)) < 2*ms1_error:
                                    num = num + 1
                                    SIM_Frg = str(round(float(msms_reference),4))
                        SIM_Frg = SIM_Frg.strip(',')
                        similarity = num / length
                        if similarity >= SIM_online or num >= frgnum_literature:
                            msms_exper = ''
                            for m in range(len(spectral)):
                                msms_exper = msms_exper + str(round(spectral[m], 4)) + ' ' + str(
                                    round(spectral_inten[m], 2)) + ';'
                            msms_exper = msms_exper.strip(';')
                            if abs(spectral_rt - literature_rt[i]) < rt_error:
                                file_raw = mzML_file.split('.')[0]
                                row_literature = row_literature + 1
                                output_literature.loc[row_literature, 'Sample_Name'] = file_raw
                                output_literature.loc[row_literature, 'Name'] = literature_name[i]
                                output_literature.loc[row_literature, 'PubChemCID'] = literature_cid[i]
                                output_literature.loc[row_literature, 'CAS'] = literature_cas[i]
                                output_literature.loc[row_literature, 'SMILES'] = literature_smiles[i]
                                output_literature.loc[row_literature, 'Formula'] = literature_formula[i]
                                output_literature.loc[row_literature, 'InChIKey'] = literature_InchiKey[i]
                                output_literature.loc[row_literature, 'm/ztheor'] = round(literature_precursorMass[i], 4)
                                output_literature.loc[row_literature, 'm/zexper'] = round(precursor, 4)
                                output_literature.loc[row_literature, 'tR/min-MSMS'] = round(spectral_rt, 2)
                                output_literature.loc[row_literature, 'tR/mintheor'] = round(literature_rt[i], 2)
                                output_literature.loc[row_literature, 'SIM-num'] = int(num)
                                output_literature.loc[row_literature, 'SIM-Frg'] = SIM_Frg
                                output_literature.loc[row_literature, 'MS/MStheor'] = literature_msms[i]
                                output_literature.loc[row_literature, 'MS/MSexper'] = msms_exper
                                output_literature.loc[row_literature, 'mass error matching'] = 'Success'
                                output_literature.loc[row_literature, 'MS/MS matched matching'] = 'Success'
                                output_literature.loc[row_literature, 'tR deviation matching'] = 'Success'
                                output_literature.loc[row_literature, 'Confidence level'] = "CL2a"
                            if abs(spectral_rt - literature_rt[i]) >= rt_error:
                                file_raw = mzML_file.split('.')[0]
                                row_literature = row_literature + 1
                                output_literature.loc[row_literature, 'Sample_Name'] = file_raw
                                output_literature.loc[row_literature, 'Name'] = literature_name[i]
                                output_literature.loc[row_literature, 'PubChemCID'] = literature_cid[i]
                                output_literature.loc[row_literature, 'CAS'] = literature_cas[i]
                                output_literature.loc[row_literature, 'SMILES'] = literature_smiles[i]
                                output_literature.loc[row_literature, 'Formula'] = literature_formula[i]
                                output_literature.loc[row_literature, 'InChIKey'] = literature_InchiKey[i]
                                output_literature.loc[row_literature, 'm/ztheor'] = round(literature_precursorMass[i], 4)
                                output_literature.loc[row_literature, 'm/zexper'] = round(precursor, 4)
                                output_literature.loc[row_literature, 'tR/min-MSMS'] = round(spectral_rt, 2)
                                output_literature.loc[row_literature, 'tR/mintheor'] = round(literature_rt[i], 2)
                                output_literature.loc[row_literature, 'SIM-num'] = int(num)
                                output_literature.loc[row_literature, 'SIM-Frg'] = SIM_Frg
                                output_literature.loc[row_literature, 'MS/MStheor'] = literature_msms[i]
                                output_literature.loc[row_literature, 'MS/MSexper'] = msms_exper
                                output_literature.loc[row_literature, 'mass error matching'] = 'Success'
                                output_literature.loc[row_literature, 'MS/MS matched matching'] = 'Success'
                                output_literature.loc[row_literature, 'Confidence level'] = "CL3"
    return output_online,output_literature


def suspect_threading1(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature):
    global output_online1,output_literature1
    output_online1,output_literature1=suspect_sceening(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature)

def suspect_threading2(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature):
    global output_online2, output_literature2
    output_online2,output_literature2=suspect_sceening(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature)

def suspect_threading3(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature):
    global output_online3, output_literature3
    output_online3,output_literature3=suspect_sceening(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature)

def suspect_threading4(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature):
    global output_online4, output_literature4
    output_online4,output_literature4=suspect_sceening(file_list,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature)

def check_onlineisotopes_etc(output_online,sample_dir,current_path,isotope_sc,rt_error,ms1_error):
    print('......正在进行同位素等核查......')
    output_online = output_online.sort_values(by=['Confidence level', 'SIM'], ascending=[True, False])
    output_online = output_online.reset_index(drop=True)
    output_online = remove_extra2(output_online)
    output_online = remove_extra(output_online)

    for index in range(output_online.shape[0]):
        rt =  output_online.loc[index, 'tR/min-MSMS']
        mz =  output_online.loc[index, 'm/zexper']
        mzML_file =  output_online.loc[index, 'Sample_Name'] + '.mzML'

        os.chdir(sample_dir)
        max_height, max_height_time, scan_times, chromatogram = get_area.get_area(mz, rt,mzML_file)
        output_online.loc[index, 'm/z-inten'] = round(max_height)
        output_online.loc[index, 'tR/min-MS'] = round(max_height_time, 2)
        eic = ''
        for x in range(len(scan_times)):
            eic = eic + str(round(scan_times[x], 3)) + ' ' + str(round(chromatogram[x], 1)) + ';'
        eic = eic.strip(';')
        output_online.loc[index, 'EIC'] = eic
        formula_eve =  output_online.loc[index, 'Formula']
        formula_eve = formula_eve.strip('+').strip('-')

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
                if abs(ratio_ther - ratio_exp) / ratio_exp < isotope_sc:
                    isotopic_num = isotopic_num + 1
                    if isotopic_num >= 1:
                        break
        if isotopic_num >= 1:
            output_online.loc[index, 'isotopic deviation matching'] = 'Success'



    # retain_list = []
    # for i in range(output_online.shape[0]):
    #     if output_online.loc[i, 'isotopic deviation matching'] == 'Success':
    #         retain_list.append(i)
    # output_online = output_online.loc[retain_list, :]
    # output_online.reset_index(drop=True, inplace=True)

    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name = str(current_date) + '_online.xlsx'
    output_online.to_excel(file_name, index=False)

def check_literatureisotopes_etc(output_literature,sample_dir,current_path,isotope_sc,rt_error,ms1_error):
    output_literature = output_literature.sort_values(by=['Confidence level', 'SIM-num'], ascending=[True, False])
    output_literature = output_literature.reset_index(drop=True)

    output_literature = remove_extra2(output_literature)
    output_literature = remove_extra(output_literature)


    for index in range(output_literature.shape[0]):
        rt =  output_literature.loc[index, 'tR/min-MSMS']
        mz =  output_literature.loc[index, 'm/zexper']
        mzML_file =  output_literature.loc[index, 'Sample_Name'] + '.mzML'
        os.chdir(sample_dir)
        max_height, max_height_time, scan_times, chromatogram = get_area.get_area(mz, rt, mzML_file)
        output_literature.loc[index, 'm/z-inten'] = round(max_height)
        output_literature.loc[index, 'tR/min-MS'] = round(max_height_time, 2)
        eic = ''
        for x in range(len(scan_times)):
            eic = eic + str(round(scan_times[x], 3)) + ' ' + str(round(chromatogram[x], 1)) + ';'
        eic = eic.strip(';')
        output_literature.loc[index, 'EIC'] = eic

        formula_eve =  output_literature.loc[index, 'Formula'].strip()
        formula_eve=formula_eve.strip('+').strip('-')

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
                if abs(ratio_ther - ratio_exp) / ratio_exp < isotope_sc:
                    isotopic_num = isotopic_num + 1
                    if isotopic_num >= 1:
                        break
        if isotopic_num >= 1:
            output_literature.loc[index, 'isotopic deviation matching'] = 'Success'

    # retain_list = []
    # for i in range(output_literature.shape[0]):
    #     if output_literature.loc[i, 'isotopic deviation matching'] == 'Success':
    #         retain_list.append(i)
    # output_literature = output_literature.loc[retain_list, :]
    # output_literature.reset_index(drop=True, inplace=True)
    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name = str(current_date) + '_literature.xlsx'
    output_literature.to_excel(file_name, index=False)

def screening(sample_dir,current_path,online_database,literature_database,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature):
    online_database = online_database.sort_values(by='m/z')
    online_database = online_database.reset_index(drop=True)
    online_precursorMass = online_database['m/z']
    online_mz_round = [round(number) for number in list(online_precursorMass)]

    online_rt=online_database['Predicted tR/min']
    online_name=online_database['Name']
    online_cid=online_database['PubChemCID']
    online_cas=online_database['CAS']
    online_formula=online_database['Formula']
    online_smiles=online_database['SMILES']
    online_InchiKey=online_database['InChIKey']
    online_msms=online_database['MS/MS']

    literature_database = literature_database.sort_values(by='m/z')
    literature_database = literature_database.reset_index(drop=True)

    literature_precursorMass = literature_database['m/z']
    literature_mz_round = [round(number) for number in list(literature_precursorMass)]

    literature_rt=literature_database['Predicted tR/min']
    literature_name=literature_database['Name']

    literature_cid=literature_database['PubChemCID']
    literature_cas=literature_database['CAS']
    literature_formula=literature_database['Formula']
    literature_smiles=literature_database['SMILES']
    literature_InchiKey=literature_database['InChIKey']
    literature_msms=literature_database['MS/MS']
    files = [f for f in os.listdir(sample_dir) if f.endswith('.mzML')]

    file_list1 = []
    file_list2 = []
    file_list3 = []
    file_list4 = []
    for x in range(0,len(files),4):
        file_list1.append(files[x])
    for x in range(1,len(files),4):
        file_list2.append(files[x])
    for x in range(2,len(files),4):
        file_list3.append(files[x])
    for x in range(3,len(files),4):
        file_list4.append(files[x])
    for file in files:
        if file not in file_list1 and file not in file_list2 and file not in file_list3 and file not in file_list4 :
            file_list4.append(file)


    t1=Thread(target=suspect_threading1,args=(file_list1,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature,))
    t2=Thread(target=suspect_threading2,args=(file_list2,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature,))
    t3=Thread(target=suspect_threading3,args=(file_list3,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature,))
    t4=Thread(target=suspect_threading4,args=(file_list4,sample_dir,online_mz_round, online_precursorMass, online_msms, online_rt,
                         online_name, online_cid, online_cas, online_smiles, online_formula, online_InchiKey,
                         literature_mz_round, literature_precursorMass, literature_msms, literature_rt, literature_name,
                         literature_cid, literature_cas, literature_smiles, literature_formula, literature_InchiKey,ms1_error, rt_error,isotope_sc, SIM_online,frgnum_literature,))
    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t1.join()
    t2.join()
    t3.join()
    t4.join()


    output_online=pd.concat([output_online1,output_online2,output_online3,output_online4],axis=0,ignore_index=True)
    output_literature=pd.concat([output_literature1,output_literature2,output_literature3,output_literature4],axis=0,ignore_index=True)


    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name1 = str(current_date) + '_original_online.xlsx'
    output_online.to_excel(file_name1, index=False)
    file_name2 = str(current_date) + '_original_literature.xlsx'
    output_literature.to_excel(file_name2, index=False)


    check_online=Thread(target=check_onlineisotopes_etc,args=(output_online,sample_dir,current_path,isotope_sc,rt_error,ms1_error,))
    check_literature=Thread(target=check_literatureisotopes_etc,args=(output_literature,sample_dir,current_path,isotope_sc, rt_error,ms1_error,))

    check_online.start()
    check_literature.start()

    check_online.join()
    check_literature.join()


