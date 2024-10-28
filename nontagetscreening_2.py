# -*- coding: utf-8 -*-
# Created by: ctt@2022
# chentiantian@dicp.ac.cn

from threading import Thread
import pandas as pd
import os
import numpy as np
from pyteomics import mzml
from collections import defaultdict
from datetime import datetime

def generate_cfi_file(CFI_list,CNL_list):
    output_cfi = pd.DataFrame()
    output_cfi['Sample_Name'] = ''
    output_cfi['Scan'] = ''
    output_cfi['EIC'] = ''
    output_cfi['tR/min'] = ''
    output_cfi['m/zexper'] = ''
    output_cfi['MS/MSexper'] = ''
    output_cfi['CFI_NUM']=''

    if len(CFI_list)>0:
        for i in range(len(CFI_list)):
            output_cfi[('CFI-'+str(CFI_list[i]))]=''


    if len(CNL_list)>0:
        for i in range(len(CNL_list)):
            output_cfi[('CNL-'+str(CNL_list[i]))]=''


    return output_cfi

def nontaget_sceening(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list):
    output_cfi=generate_cfi_file(CFI_list,CNL_list)
    CFI_num=len(CFI_list)+len(CNL_list)
    row=-1

    for mzML_file in file_list:
        print(mzML_file)
        print('  ')
        os.chdir(sample_dir)
        file = mzml.read(mzML_file)
        for spe in file:
            if spe['ms level'] == 2:
                spectral_rt = spe['scanList']['scan'][0]['scan start time']
                precursor = spe['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                scan=int(spe['id'].split('scan=')[1])
                spectral = spe['m/z array']
                spectral_inten = spe['intensity array'] / max(spe['intensity array'])
                spectral_inten = spectral_inten / max(spectral_inten)
                keep1 = np.where(spectral_inten >= inten_filter)
                spectral = spectral[keep1]
                spectral_inten = spectral_inten[keep1]
                zero_length_list=[0] * CFI_num
                for l in range(len(spectral)):
                    for i in range(len(CFI_list)):
                        if 10 ** 6 * abs(CFI_list[i] - spectral[l]) / 300 < ms2_error:
                            zero_length_list[i]=1
                    for n in range(len(CNL_list)):
                        if 10 ** 6 * abs(CNL_list[n] -precursor + spectral[l])/300 < ms2_error:
                            zero_length_list[n+len(CFI_list)] = 1
                if CFI_num>4 and sum(zero_length_list)>=2:
                    msms_exper = ''
                    for m in range(len(spectral)):
                        msms_exper = msms_exper + str(round(spectral[m], 4)) + ' ' + str(round(spectral_inten[m], 3)) + ';'
                    msms_exper = msms_exper.strip(';')
                    row=row+1
                    output_cfi.loc[row,'Sample_Name']=(mzML_file.split('.')[0])
                    output_cfi.loc[row, 'Scan'] = scan
                    output_cfi.loc[row, 'tR/min'] = round(spectral_rt, 2)
                    output_cfi.loc[row, 'm/zexper'] = round(precursor, 4)
                    output_cfi.loc[row, 'MS/MSexper'] = msms_exper
                    output_cfi.loc[row, 'CFI_NUM'] = sum(zero_length_list)
                    for w in range(len(zero_length_list)):
                        if zero_length_list[w]==1:
                            output_cfi.iloc[row, 7+w] = 1
                if CFI_num<=4 and sum(zero_length_list)>=1:
                    msms_exper = ''
                    for m in range(len(spectral)):
                        msms_exper = msms_exper + str(round(spectral[m], 4)) + ' ' + str(round(spectral_inten[m], 3)) + ';'
                    msms_exper = msms_exper.strip(';')
                    row=row+1
                    output_cfi.loc[row,'Sample_Name']=(mzML_file.split('.')[0])
                    output_cfi.loc[row, 'Scan'] = scan
                    output_cfi.loc[row, 'tR/min'] = round(spectral_rt, 2)
                    output_cfi.loc[row, 'm/zexper'] = round(precursor, 4)
                    output_cfi.loc[row, 'MS/MSexper'] = msms_exper
                    output_cfi.loc[row, 'CFI_NUM'] = sum(zero_length_list)
                    for w in range(len(zero_length_list)):
                        if zero_length_list[w]==1:
                            output_cfi.iloc[row, 7+w] = 1
    return output_cfi




def screening_threading1(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list):
    global output_cfi1
    output_cfi1=nontaget_sceening(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list)

def screening_threading2(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list):
    global output_cfi2
    output_cfi2=nontaget_sceening(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list)

def screening_threading3(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list):
    global output_cfi3
    output_cfi3=nontaget_sceening(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list)

def screening_threading4(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list):
    global output_cfi4
    output_cfi4=nontaget_sceening(file_list,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list)

def check_removal_etc(output_cfi,sample_dir,current_path,cls,cl_name,peak_df):
    print('......正在进行去重......')
    l_list = []
    for i in range(output_cfi.shape[0]):
        for l in range(i, output_cfi.shape[0]):
            if i != l and l not in l_list:
                if 10 ** 6 * abs(output_cfi.loc[i, 'm/zexper'] - output_cfi.loc[l, 'm/zexper']) / 300 < 20 and abs(output_cfi.loc[i, 'tR/min'] - output_cfi.loc[l, 'tR/min']) < 0.2:
                    mz = output_cfi.loc[i, 'm/zexper']
                    rt = output_cfi.loc[i, 'tR/min']
                    output_cfi.loc[l, 'm/zexper'] = mz
                    output_cfi.loc[l, 'tR/min'] = rt
                    l_list.append(l)
    output_cfi = output_cfi.drop_duplicates(subset=cl_name)
    output_cfi = output_cfi.reset_index(drop=True)

    for i in range(output_cfi.shape[0]):
        for l in range(i, output_cfi.shape[0]):
            if i != l:
                if 10 ** 6 * abs(output_cfi.loc[i, 'm/zexper'] - output_cfi.loc[l, 'm/zexper']) / 300 < 20:
                    mz = output_cfi.loc[i, 'm/zexper']
                    output_cfi.loc[l, 'm/zexper'] = mz

    remove_mz_list = []
    mz_statistics = defaultdict(list)
    for n, v in enumerate(output_cfi['m/zexper']):
        mz_statistics[v].append(n)
    for e in mz_statistics:
        if len(mz_statistics[e]) >= 8:
            remove_mz_list.extend(mz_statistics[e])

    retain_list = []
    for i in range(output_cfi.shape[0]):
        if i not in remove_mz_list:
            retain_list.append(i)
    output_cfi = output_cfi.loc[retain_list, :]
    output_cfi.reset_index(drop=True, inplace=True)


    retain_list = []
    peak_mz = peak_df['m/z']
    peak_rt = peak_df['tR/min']
    for i in range(len(output_cfi['m/zexper'])):
        for l in range(len(peak_mz)):
            if 10 ** 6 * abs(peak_mz[l] - output_cfi.loc[i, 'm/zexper']) / peak_mz[l] < 20 and abs(peak_rt[l] - output_cfi.loc[i, 'tR/min']) < 0.2:
                retain_list.append(i)
                break

    output_cfi = output_cfi.loc[retain_list, :]
    output_cfi.reset_index(drop=True, inplace=True)

    for w in range(output_cfi.shape[0]):
        target_mz=output_cfi.loc[w, 'm/zexper']
        rt = output_cfi.loc[w, 'tR/min']
        file_name=output_cfi.loc[w, 'Sample_Name']
        os.chdir(sample_dir)
        file = mzml.read(file_name+'.mzML')
        eic = ''
        for spe in file:
            if spe['ms level'] == 1:
                scan_time = spe['scanList']['scan'][0]['scan start time']
                if scan_time - rt < 0.2:
                    if abs(scan_time - rt) < 0.2:
                        mz_array = spe['m/z array']
                        intensity_array = spe['intensity array']
                        keep = np.where(intensity_array > 0)
                        mz_array = mz_array[keep]
                        intensity_array = intensity_array[keep]
                        for m in range(len(mz_array)):
                            if 10 ** 6 * abs(target_mz - mz_array[m]) / target_mz < 10:
                                eic = eic + str(round(scan_time*60, 3)) + ' ' + str(round(intensity_array[m], 1)) + ';'
                                break
                if scan_time - rt > 0.2:
                    break
        eic = eic.strip(';')
        output_cfi.loc[w, 'EIC']=eic

    os.chdir(current_path)
    os.chdir('results')
    current_date = datetime.now().date()
    file_name = str(current_date) + '_cfi_'+str(cls)+'.xlsx'
    output_cfi.to_excel(file_name, index=False)

def screening(sample_dir,current_path,cfi_dadabase,peaklistpath,ms2_error,inten_filter):

    peak_df=pd.read_excel(peaklistpath)

    cfi_category=cfi_dadabase['Category']
    cfi_CFI=cfi_dadabase['CFI']
    cfi_CNL=cfi_dadabase['CNL']

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
    for c in range(len(cfi_category)):
        CFI_list=[]
        CNL_list = []
        cls=cfi_category[c]
        print('......' + str(cls) + '......')
        CFI=cfi_CFI[c]
        if pd.isnull(CFI) == False:
            try:
                CFI=CFI.strip(',').split(',')
                for cfi_e in CFI:
                    CFI_list.append(float(cfi_e))
            except:
                CFI = CFI
                CFI_list.append(float(CFI))

        CNL = cfi_CNL[c]
        if pd.isnull(CNL) == False:
            try:
                CNL = CNL.strip(',').split(',')
                for cnl_e in CNL:
                    CNL_list.append(float(cnl_e))
            except:
                CNL = CNL
                CNL_list.append(float(CNL))
        cl_name = ['tR/min', 'm/zexper']

        if len(CFI_list) > 0:
            for i in range(len(CFI_list)):
                cl_name.append(('CFI-' + str(CFI_list[i])))


        if len(CNL_list) > 0:
            for i in range(len(CNL_list)):
                cl_name.append(('CNL-' + str(CNL_list[i])))



        t1=Thread(target=screening_threading1,args=(file_list1,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list,))
        t2=Thread(target=screening_threading2,args=(file_list2,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list,))
        t3=Thread(target=screening_threading3,args=(file_list3,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list,))
        t4=Thread(target=screening_threading4,args=(file_list4,sample_dir,ms2_error,inten_filter,CFI_list,CNL_list,))
        t1.start()
        t2.start()
        t3.start()
        t4.start()
        t1.join()
        t2.join()
        t3.join()
        t4.join()

        output_cfi=pd.concat([output_cfi1,output_cfi2,output_cfi3,output_cfi4],axis=0,ignore_index=True)
        os.chdir(current_path)
        os.chdir('results')
        current_date = datetime.now().date()
        file_name = str(current_date) +'_original_cfi_'+str(cls)+'.xlsx'
        output_cfi.to_excel(file_name, index=False)

        check=Thread(target=check_removal_etc,args=(output_cfi,sample_dir,current_path,cls,cl_name,peak_df,))

        check.start()
        check.join()


