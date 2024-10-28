# -*- coding: utf-8 -*-
# Created by: ctt@2022
# chentiantian@dicp.ac.cn
import formula

def isotopic(formula_eve,mz):
    f = formula.Formula(formula_eve)
    M, p = f.get_isotopic_envelope(4)
    M = list(M)
    p = list(p)
    max_index = p.index(max(p[1:]))
    M_2 = M[max_index] + mz - M[0]
    p_2 = p[max_index]
    M_1 = M[0] + mz - M[0]
    p_1 = p[0]
    ratio_ther = p_1/p_2
    return ratio_ther,M_1,M_2
