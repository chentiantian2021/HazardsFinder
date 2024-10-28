# -*- coding: utf-8 -*-
# Created by: ctt@2022
# chentiantian@dicp.ac.cn


from pyteomics import mzml
import numpy as np

def get_area(target_mz,rt,mzML_file):
    max_height=0
    max_height_time=0
    intensities = []
    times = []
    mzml_scr = mzml.read(mzML_file)
    for spectrum in mzml_scr:
        if spectrum['ms level'] == 1:
            match_num=0
            scan_time = spectrum['scanList']['scan'][0]['scan start time']
            if scan_time - rt < 0.5:
                if abs(scan_time - rt) < 0.5:
                    mz_array = spectrum['m/z array']
                    intensity_array = spectrum['intensity array']
                    keep = np.where(intensity_array > 0)
                    mz_array = mz_array[keep]
                    intensity_array = intensity_array[keep]
                    for w in range(len(mz_array)):
                        if 10 ** 6 * abs(target_mz - mz_array[w]) / target_mz < 10:
                            intensities.append(intensity_array[w])
                            times.append(scan_time*60)
                            match_num=match_num+1
                    if match_num == 0:
                        intensities.append(0)
                        times.append(scan_time*60)
            if scan_time-rt >=0.5:
                break



    # 计算面积和最大高度
    if len(intensities)>0:
        max_index=intensities.index(max(intensities))
        max_height = intensities[max_index]
        max_height_time=times[max_index]/60
    return max_height,max_height_time,times,intensities

