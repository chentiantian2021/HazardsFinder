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
from collections import defaultdict
from datetime import datetime



def remove_extra(result_df):
    remove_list=[]
    for i in range(result_df.shape[0]):
        for g in range(i,result_df.shape[0]):
            if i != g:
                if result_df.loc[i, 'InChIKey'] == result_df.loc[g, 'InChIKey'] and abs(result_df.loc[i, 'tR/min-MSMS'] -result_df.loc[g, 'tR/min-MSMS']) <= 0.21:
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
    new_df = result_df[['Sample_Name', 'InChIKey','tR/min-MSMS']]
    specific_new_df = new_df.drop_duplicates()
    specific_new_df.reset_index(drop=True, inplace=True)

    remove_pubchemcid = []
    for i in range(len(sample_list)):
        pubchemcid=[]
        for l in range(specific_new_df.shape[0]):
            if sample_list[i]==specific_new_df.loc[l,'Sample_Name']:
                pubchemcid.append(specific_new_df.loc[l,'InChIKey'])
        pubchemcid_statistics = defaultdict(list)
        for n, v in enumerate(pubchemcid):
            pubchemcid_statistics[v].append(n)
        for e in pubchemcid_statistics:
            if len(pubchemcid_statistics[e])>=8:
                remove_pubchemcid.append(e)

    remove_list=[]
    for e in remove_pubchemcid:
        for i in range(result_df.shape[0]):
            if str(e) ==str(result_df.loc[i,'InChIKey']):
                remove_list.append(i)

    retain_list=[]
    for i in range(result_df.shape[0]):
        if i not in remove_list:
            retain_list.append(i)
    result_df=result_df.loc[retain_list,:]
    result_df.reset_index(drop=True, inplace=True)
    return result_df


def generate_tps_file():
    output_tps = pd.DataFrame()
    output_tps['Sample_Name'] = ''
    output_tps['Name'] = ''
    output_tps['PubChemCID'] = ''
    output_tps['CAS'] = ''
    output_tps['Formula'] = ''
    output_tps['SMILES'] = ''
    output_tps['InChIKey'] = ''
    output_tps['Predecessor_PubchemCID'] = ''
    output_tps['Predecessor_Name'] = ''
    output_tps['Predecessor_SMILES'] = ''

    output_tps['m/ztheor'] = ''
    output_tps['m/zexper'] = ''
    output_tps['m/z-inten'] = ''
    output_tps['tR/mintheor'] = ''
    output_tps['tR/min-MS'] = ''
    output_tps['tR/min-MSMS'] = ''
    output_tps['EIC'] = ''
    output_tps['SIM-num'] = ''
    output_tps['SIM-Frg'] = ''
    output_tps['MS/MStheor'] = ''
    output_tps['MS/MSexper'] = ''
    output_tps['mass error matching'] = ''
    output_tps['isotopic deviation matching'] = ''
    output_tps['tR deviation matching'] = ''
    output_tps['MS/MS matched matching'] = ''
    output_tps['Confidence level'] = ''

    return output_tps

def nontaget_sceening(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid):
    output_tps=generate_tps_file()
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
                tps_indices = [index for index, value in enumerate(tps_mz_round) if value in values_to_find]
                for i in tps_indices:
                    if 10 ** 6 * abs(precursor - tps_precursorMass[i]) / precursor < ms1_error:
                        num = 0
                        spectral = spe['m/z array']
                        spectral_inten = spe['intensity array'] / max(spe['intensity array'])
                        spectral_inten = spectral_inten / max(spectral_inten)
                        keep1 = np.where(spectral_inten >= 0.01)
                        spectral = spectral[keep1]
                        spectral_inten=spectral_inten[keep1]
                        SIM_Frg = ''
                        try:
                            msms_reference = tps_msms[i].split(';')
                            for m in range(len(spectral)):
                                for n in range(len(msms_reference)):
                                    if 10 ** 6 * abs(spectral[m] - float(msms_reference[n].split(' ')[0])) / spectral[m] < ms1_error:
                                        num = num + 1
                                        SIM_Frg = SIM_Frg + str(round(float(msms_reference[n].split(' ')[0]), 4)) + ','
                        except:
                            pass
                        if num >= frgnum_cfmid:
                            msms_exper = ''
                            SIM_Frg = SIM_Frg.strip(',')
                            for m in range(len(spectral)):
                                msms_exper = msms_exper + str(round(spectral[m], 4)) + ' ' + str(round(spectral_inten[m], 3)) + ';'
                            msms_exper = msms_exper.strip(';')
                            if abs(spectral_rt - tps_rt[i]) < rt_error:
                                file_raw = mzML_file.split('.')[0]
                                row = row + 1
                                output_tps.loc[row, 'Sample_Name'] = file_raw
                                output_tps.loc[row, 'Name'] = tps_name[i]
                                output_tps.loc[row, 'PubChemCID'] = tps_cid[i]
                                output_tps.loc[row, 'CAS'] = tps_cas[i]
                                output_tps.loc[row, 'SMILES'] = tps_smiles[i]
                                output_tps.loc[row, 'Formula'] = tps_formula[i]
                                output_tps.loc[row, 'InChIKey'] = tps_InchiKey[i]
                                output_tps.loc[row, 'm/ztheor'] = round(tps_precursorMass[i], 4)
                                output_tps.loc[row, 'm/zexper'] = round(precursor, 4)
                                output_tps.loc[row, 'tR/min-MSMS'] = round(spectral_rt, 2)
                                output_tps.loc[row, 'tR/mintheor'] = round(tps_rt[i], 2)
                                output_tps.loc[row, 'SIM-num'] = int(num)
                                output_tps.loc[row, 'SIM-Frg'] = SIM_Frg
                                output_tps.loc[row, 'MS/MStheor'] = tps_msms[i]
                                output_tps.loc[row, 'MS/MSexper'] = msms_exper
                                output_tps.loc[row, 'mass error matching'] = 'Success'
                                output_tps.loc[row, 'MS/MS matched matching'] = 'Success'
                                output_tps.loc[row, 'tR deviation matching'] = 'Success'
                                output_tps.loc[row, 'Confidence level'] = "CL3"
                                output_tps.loc[row, 'Predecessor_PubchemCID'] = tps_Predecessor_PubchemCID[i]
                                output_tps.loc[row, 'Predecessor_Name'] = tps_Predecessor_Name[i]
                                output_tps.loc[row, 'Predecessor_SMILES'] = tps_Predecessor_SMILES[i]
                            if abs(spectral_rt - tps_rt[i]) >= rt_error:
                                file_raw = mzML_file.split('.')[0]
                                row = row + 1
                                output_tps.loc[row, 'Sample_Name'] = file_raw
                                output_tps.loc[row, 'Name'] = tps_name[i]
                                output_tps.loc[row, 'PubChemCID'] = tps_cid[i]
                                output_tps.loc[row, 'CAS'] = tps_cas[i]
                                output_tps.loc[row, 'SMILES'] = tps_smiles[i]
                                output_tps.loc[row, 'Formula'] = tps_formula[i]
                                output_tps.loc[row, 'InChIKey'] = tps_InchiKey[i]
                                output_tps.loc[row, 'm/ztheor'] = round(tps_precursorMass[i], 4)
                                output_tps.loc[row, 'm/zexper'] = round(precursor, 4)
                                output_tps.loc[row, 'tR/min-MSMS'] = round(spectral_rt, 2)
                                output_tps.loc[row, 'tR/mintheor'] = round(tps_rt[i], 2)
                                output_tps.loc[row, 'SIM-num'] = int(num)
                                output_tps.loc[row, 'SIM-Frg'] = SIM_Frg
                                output_tps.loc[row, 'MS/MStheor'] = tps_msms[i]
                                output_tps.loc[row, 'MS/MSexper'] = msms_exper
                                output_tps.loc[row, 'mass error matching'] = 'Success'
                                output_tps.loc[row, 'MS/MS matched matching'] = 'Success'
                                output_tps.loc[row, 'tR deviation matching'] = 'Success'
                                output_tps.loc[row, 'Confidence level'] = "CL4"
                                output_tps.loc[row, 'Predecessor_PubchemCID'] = tps_Predecessor_PubchemCID[i]
                                output_tps.loc[row, 'Predecessor_Name'] = tps_Predecessor_Name[i]
                                output_tps.loc[row, 'Predecessor_SMILES'] = tps_Predecessor_SMILES[i]
    return output_tps

def screening_threading1(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid):
    global output_tps1
    output_tps1=nontaget_sceening(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid)

def screening_threading2(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid):
    global output_tps2
    output_tps2=nontaget_sceening(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid)

def screening_threading3(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid):
    global output_tps3
    output_tps3=nontaget_sceening(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid)

def screening_threading4(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid):
    global output_tps4
    output_tps4=nontaget_sceening(file_list,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid)

def check_tpsisotopes_etc(output_tps,sample_dir,current_path,isotope_sc,rt_error,ms1_error):
    print('......正在进行同位素等核查......')
    output_tps = output_tps.sort_values(by=['Confidence level', 'SIM-num'], ascending=[True, False])
    output_tps = output_tps.reset_index(drop=True)
    output_tps = remove_extra2(output_tps)
    output_tps = remove_extra(output_tps)
    output_tps.reset_index(drop=True, inplace=True)
    for index in range(output_tps.shape[0]):
        rt =  output_tps.loc[index, 'tR/min-MSMS']
        mz =  output_tps.loc[index, 'm/zexper']
        mzML_file =  output_tps.loc[index, 'Sample_Name'] + '.mzML'
        os.chdir(sample_dir)
        max_height, max_height_time, scan_times, chromatogram = get_area.get_area(mz, rt, mzML_file)
        output_tps.loc[index, 'tR/min-MS'] = round(max_height_time, 2)
        output_tps.loc[index, 'm/z-inten'] = round(max_height)
        eic = ''
        for x in range(len(scan_times)):
            eic = eic + str(round(scan_times[x], 3)) + ' ' + str(round(chromatogram[x], 1)) + ';'
        eic = eic.strip(';')
        output_tps.loc[index, 'EIC'] = eic
        formula_eve =  output_tps.loc[index, 'Formula']
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
                ratio_exp = p_1_exp / p_2_exp
                if abs(ratio_ther - ratio_exp) / ratio_exp < isotope_sc:
                    isotopic_num = isotopic_num + 1
                    if isotopic_num >= 1:
                        break
        if isotopic_num >= 1:
            output_tps.loc[index, 'isotopic deviation matching'] = 'Success'

    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name = str(current_date) + '_tps.xlsx'
    output_tps.to_excel(file_name, index=False)

def screening(sample_dir,current_path,tps_dadabase,peaklistpath,ms1_error, rt_error,isotope_sc,frgnum_cfmid):
    peak_df=pd.read_excel(peaklistpath)
    tps_precursorMass = tps_dadabase['m/z']
    peak_mz=peak_df['m/z']
    retain_list = []
    for i in range(len(tps_precursorMass)):
        for l in range(len(peak_mz)):
            if 10 ** 6 * abs(peak_mz[l] - tps_precursorMass[i]) / tps_precursorMass[i] < 20:
                retain_list.append(i)
                break
    tps_dadabase = tps_dadabase.loc[retain_list, :]
    tps_dadabase.reset_index(drop=True, inplace=True)


    tps_dadabase = tps_dadabase.sort_values(by='m/z')
    tps_dadabase = tps_dadabase.reset_index(drop=True)

    tps_precursorMass = tps_dadabase['m/z']

    tps_mz_round = [round(number) for number in list(tps_precursorMass)]

    tps_rt=tps_dadabase['Predicted tR/min']
    tps_msms=tps_dadabase['Predicted MS/MS']
    tps_name=tps_dadabase['Name']
    tps_cid=tps_dadabase['PubChemCID']
    tps_cas=tps_dadabase['CAS']
    tps_formula=tps_dadabase['Formula']
    tps_smiles=tps_dadabase['SMILES']
    tps_InchiKey=tps_dadabase['InChIKey']
    tps_Predecessor_PubchemCID=tps_dadabase['Predecessor_PubchemCID']
    tps_Predecessor_Name=tps_dadabase['Predecessor_Name']
    tps_Predecessor_SMILES=tps_dadabase['Predecessor_SMILES']

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


    t1=Thread(target=screening_threading1,args=(file_list1,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid,))
    t2=Thread(target=screening_threading2,args=(file_list2,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid,))
    t3=Thread(target=screening_threading3,args=(file_list3,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid,))
    t4=Thread(target=screening_threading4,args=(file_list4,sample_dir,tps_mz_round,tps_precursorMass,tps_msms,tps_rt,tps_name,tps_cid,tps_cas,tps_smiles,tps_formula,tps_InchiKey,tps_Predecessor_PubchemCID,tps_Predecessor_Name,tps_Predecessor_SMILES,ms1_error, rt_error,isotope_sc,frgnum_cfmid,))
    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t1.join()
    t2.join()
    t3.join()
    t4.join()

    output_tps=pd.concat([output_tps1,output_tps2,output_tps3,output_tps4],axis=0,ignore_index=True)
    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name = str(current_date) + '_original_tps.xlsx'
    output_tps.to_excel(file_name, index=False)

    check=Thread(target=check_tpsisotopes_etc,args=(output_tps,sample_dir,current_path,isotope_sc,rt_error,ms1_error,))

    check.start()
    check.join()


