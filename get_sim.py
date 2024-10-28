# -*- coding: utf-8 -*-
# Created by: ctt@2022
# chentiantian@dicp.ac.cn
import numpy as np
import spectral_entropy
import tools



def sim_cal_1(msms_1,msms_2):
    try:
        msms_1s=msms_1.split(';')
        msms_2s=msms_2.split(';')
    except:
        msms_1s = msms_1
        msms_2s = msms_2
    spec_1 = []
    for n in range(len(msms_1s)):
        spec0 = []
        spec0.append(msms_1s[n].split(' ')[0])
        spec0.append(msms_1s[n].split(' ')[1])
        spec_1.append(spec0)
    spec_1 = np.array(spec_1, dtype=np.float32)
    spec_2 = []
    for n in range(len(msms_2s)):
        spec0 = []
        spec0.append(msms_2s[n].split(' ')[0])
        spec0.append(msms_2s[n].split(' ')[1])
        spec_2.append(spec0)
    spec_2 = np.array(spec_2, dtype=np.float32)
    spec_2_clean = tools.clean_spectrum(spec_2, noise_removal=0.01, ms2_ppm=10)
    spec_1_clean = tools.clean_spectrum(spec_1, noise_removal=0.01, ms2_ppm=10)
    similarity = spectral_entropy.calculate_entropy_similarity(spec_1_clean, spec_2_clean, ms2_da=0.01)
    return similarity

def sim_cal_2(precursor,spectral,spectral_inten, msms_standard):
    spec_query = []
    msms_exper = ''
    for m in range(len(spectral)):
        spec0 = []
        spec0.append(spectral[m])
        spec0.append(spectral_inten[m])
        spec_query.append(spec0)
        if spectral_inten[m] >= 0.01:
            msms_exper = msms_exper + str(round(spectral[m], 4)) + ' ' + str(round(spectral_inten[m], 2)) + ';'
    msms_exper = msms_exper.strip(';')
    spec_query = np.array(spec_query, dtype=np.float32)
    try:
        msms_reference = msms_standard.split(';')
        spec_reference = []
        for n in range(len(msms_reference)):
            spec0 = []
            spec0.append(float(msms_reference[n].split(' ')[0]))
            spec0.append(float(msms_reference[n].split(' ')[1]))
            spec_reference.append(spec0)
    except:
        try:
            msms_reference = msms_standard.split(';')
            spec_reference = []
            for n in range(len(msms_reference)-1):
                spec0 = []
                spec0.append(float(msms_reference[n].split(' ')[0]))
                spec0.append(float(msms_reference[n].split(' ')[1]))
                spec_reference.append(spec0)
        except:
            msms_reference = msms_standard
            spec_reference = []
            spec0 = []
            spec0.append(float(msms_reference.split(' ')[0]))
            spec0.append(float(msms_reference.split(' ')[1]))
            spec_reference.append(spec0)
    spec_reference = np.array(spec_reference, dtype=np.float32)
    spec_query_clean = tools.clean_spectrum(spec_query,max_mz=precursor+1, noise_removal=0.01,
                                            ms2_ppm=10)
    spec_reference_clean = tools.clean_spectrum(spec_reference, max_mz=precursor+1,
                                                noise_removal=0.01, ms2_ppm=10)
    similarity = spectral_entropy.calculate_entropy_similarity(spec_query_clean, spec_reference_clean,
                                                               ms2_da=0.01)

    return similarity,msms_exper