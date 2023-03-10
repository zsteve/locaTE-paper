#####################################################
import autograd.numpy as np
# This file is created automatically
def Model(Y,t,pars):
    # Parameters
    a_g1_g4 = pars[0]
    a_g1_g4_g5 = pars[1]
    a_g1_g5 = pars[2]
    a_g2_g1 = pars[3]
    a_g3_g2 = pars[4]
    a_g4_g3 = pars[5]
    a_g4_g3_g5 = pars[6]
    a_g4_g4 = pars[7]
    a_g4_g4_g3 = pars[8]
    a_g4_g4_g3_g5 = pars[9]
    a_g4_g4_g5 = pars[10]
    a_g4_g5 = pars[11]
    a_g5_g3 = pars[12]
    a_g5_g4 = pars[13]
    a_g5_g4_g3 = pars[14]
    a_g5_g4_g5 = pars[15]
    a_g5_g4_g5_g3 = pars[16]
    a_g5_g5 = pars[17]
    a_g5_g5_g3 = pars[18]
    a_g6_g4 = pars[19]
    a_g6_g4_g5 = pars[20]
    a_g6_g4_g6 = pars[21]
    a_g6_g4_g6_g5 = pars[22]
    a_g6_g4_g6_g9 = pars[23]
    a_g6_g4_g6_g9_g5 = pars[24]
    a_g6_g4_g9 = pars[25]
    a_g6_g4_g9_g5 = pars[26]
    a_g6_g5 = pars[27]
    a_g6_g6 = pars[28]
    a_g6_g6_g5 = pars[29]
    a_g6_g6_g9 = pars[30]
    a_g6_g6_g9_g5 = pars[31]
    a_g6_g7 = pars[32]
    a_g6_g7_g4 = pars[33]
    a_g6_g7_g4_g5 = pars[34]
    a_g6_g7_g4_g6 = pars[35]
    a_g6_g7_g4_g6_g5 = pars[36]
    a_g6_g7_g4_g6_g9 = pars[37]
    a_g6_g7_g4_g6_g9_g5 = pars[38]
    a_g6_g7_g4_g9 = pars[39]
    a_g6_g7_g4_g9_g5 = pars[40]
    a_g6_g7_g5 = pars[41]
    a_g6_g7_g6 = pars[42]
    a_g6_g7_g6_g5 = pars[43]
    a_g6_g7_g6_g9 = pars[44]
    a_g6_g7_g6_g9_g5 = pars[45]
    a_g6_g7_g9 = pars[46]
    a_g6_g7_g9_g5 = pars[47]
    a_g6_g9 = pars[48]
    a_g6_g9_g5 = pars[49]
    a_g7_g4 = pars[50]
    a_g7_g4_g5 = pars[51]
    a_g7_g4_g5_g6 = pars[52]
    a_g7_g4_g5_g8 = pars[53]
    a_g7_g4_g5_g8_g6 = pars[54]
    a_g7_g4_g6 = pars[55]
    a_g7_g4_g8 = pars[56]
    a_g7_g4_g8_g6 = pars[57]
    a_g7_g5 = pars[58]
    a_g7_g5_g6 = pars[59]
    a_g7_g5_g8 = pars[60]
    a_g7_g5_g8_g6 = pars[61]
    a_g7_g6 = pars[62]
    a_g7_g8 = pars[63]
    a_g7_g8_g6 = pars[64]
    a_g8_g4 = pars[65]
    a_g8_g5 = pars[66]
    a_g8_g5_g4 = pars[67]
    a_g8_g5_g9 = pars[68]
    a_g8_g5_g9_g4 = pars[69]
    a_g8_g7 = pars[70]
    a_g8_g7_g4 = pars[71]
    a_g8_g7_g5 = pars[72]
    a_g8_g7_g5_g4 = pars[73]
    a_g8_g7_g5_g9 = pars[74]
    a_g8_g7_g5_g9_g4 = pars[75]
    a_g8_g7_g9 = pars[76]
    a_g8_g7_g9_g4 = pars[77]
    a_g8_g9 = pars[78]
    a_g8_g9_g4 = pars[79]
    a_g9_g4 = pars[80]
    a_g9_g4_g5 = pars[81]
    a_g9_g4_g6 = pars[82]
    a_g9_g4_g6_g5 = pars[83]
    a_g9_g4_g6_g8 = pars[84]
    a_g9_g4_g6_g8_g5 = pars[85]
    a_g9_g4_g6_g9 = pars[86]
    a_g9_g4_g6_g9_g5 = pars[87]
    a_g9_g4_g6_g9_g8 = pars[88]
    a_g9_g4_g6_g9_g8_g5 = pars[89]
    a_g9_g4_g8 = pars[90]
    a_g9_g4_g8_g5 = pars[91]
    a_g9_g4_g9 = pars[92]
    a_g9_g4_g9_g5 = pars[93]
    a_g9_g4_g9_g8 = pars[94]
    a_g9_g4_g9_g8_g5 = pars[95]
    a_g9_g5 = pars[96]
    a_g9_g6 = pars[97]
    a_g9_g6_g5 = pars[98]
    a_g9_g6_g8 = pars[99]
    a_g9_g6_g8_g5 = pars[100]
    a_g9_g6_g9 = pars[101]
    a_g9_g6_g9_g5 = pars[102]
    a_g9_g6_g9_g8 = pars[103]
    a_g9_g6_g9_g8_g5 = pars[104]
    a_g9_g8 = pars[105]
    a_g9_g8_g5 = pars[106]
    a_g9_g9 = pars[107]
    a_g9_g9_g5 = pars[108]
    a_g9_g9_g8 = pars[109]
    a_g9_g9_g8_g5 = pars[110]
    alpha_g1 = pars[111]
    alpha_g2 = pars[112]
    alpha_g3 = pars[113]
    alpha_g4 = pars[114]
    alpha_g5 = pars[115]
    alpha_g6 = pars[116]
    alpha_g7 = pars[117]
    alpha_g8 = pars[118]
    alpha_g9 = pars[119]
    k_g1 = pars[120]
    k_g2 = pars[121]
    k_g3 = pars[122]
    k_g4 = pars[123]
    k_g5 = pars[124]
    k_g6 = pars[125]
    k_g7 = pars[126]
    k_g8 = pars[127]
    k_g9 = pars[128]
    l_p_g1 = pars[129]
    l_p_g2 = pars[130]
    l_p_g3 = pars[131]
    l_p_g4 = pars[132]
    l_p_g5 = pars[133]
    l_p_g6 = pars[134]
    l_p_g7 = pars[135]
    l_p_g8 = pars[136]
    l_p_g9 = pars[137]
    l_x_g1 = pars[138]
    l_x_g2 = pars[139]
    l_x_g3 = pars[140]
    l_x_g4 = pars[141]
    l_x_g5 = pars[142]
    l_x_g6 = pars[143]
    l_x_g7 = pars[144]
    l_x_g8 = pars[145]
    l_x_g9 = pars[146]
    m_g1 = pars[147]
    m_g2 = pars[148]
    m_g3 = pars[149]
    m_g4 = pars[150]
    m_g5 = pars[151]
    m_g6 = pars[152]
    m_g7 = pars[153]
    m_g8 = pars[154]
    m_g9 = pars[155]
    n_g1 = pars[156]
    n_g2 = pars[157]
    n_g3 = pars[158]
    n_g4 = pars[159]
    n_g5 = pars[160]
    n_g6 = pars[161]
    n_g7 = pars[162]
    n_g8 = pars[163]
    n_g9 = pars[164]
    r_g1 = pars[165]
    r_g2 = pars[166]
    r_g3 = pars[167]
    r_g4 = pars[168]
    r_g5 = pars[169]
    r_g6 = pars[170]
    r_g7 = pars[171]
    r_g8 = pars[172]
    r_g9 = pars[173]
    sigmaH_g1 = pars[174]
    sigmaH_g2 = pars[175]
    sigmaH_g3 = pars[176]
    sigmaH_g4 = pars[177]
    sigmaH_g5 = pars[178]
    sigmaH_g6 = pars[179]
    sigmaH_g7 = pars[180]
    sigmaH_g8 = pars[181]
    sigmaH_g9 = pars[182]
    # Variables
    x_g1 = Y[0]
    p_g1 = Y[1]
    x_g2 = Y[2]
    p_g2 = Y[3]
    x_g3 = Y[4]
    p_g3 = Y[5]
    x_g4 = Y[6]
    p_g4 = Y[7]
    x_g5 = Y[8]
    p_g5 = Y[9]
    x_g6 = Y[10]
    p_g6 = Y[11]
    x_g7 = Y[12]
    p_g7 = Y[13]
    x_g8 = Y[14]
    p_g8 = Y[15]
    x_g9 = Y[16]
    p_g9 = Y[17]
    dx_g1 = m_g1*(( alpha_g1 + a_g1_g4*(p_g4/k_g4)**n_g4 + a_g1_g5*(p_g5/k_g5)**n_g5 + a_g1_g4_g5*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 ))-l_x_g1*x_g1
    dp_g1 = r_g1*x_g1- l_p_g1*p_g1
    dx_g2 = m_g2*(( alpha_g2 + a_g2_g1*(p_g1/k_g1)**n_g1 )/( 1 +(p_g1/k_g1)**n_g1 ))-l_x_g2*x_g2
    dp_g2 = r_g2*x_g2- l_p_g2*p_g2
    dx_g3 = m_g3*(( alpha_g3 + a_g3_g2*(p_g2/k_g2)**n_g2 )/( 1 +(p_g2/k_g2)**n_g2 ))-l_x_g3*x_g3
    dp_g3 = r_g3*x_g3- l_p_g3*p_g3
    dx_g4 = m_g4*(( alpha_g4 + a_g4_g4*(p_g4/k_g4)**n_g4 + a_g4_g3*(p_g3/k_g3)**n_g3 + a_g4_g5*(p_g5/k_g5)**n_g5 + a_g4_g4_g3*(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 + a_g4_g4_g5*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 + a_g4_g3_g5*(p_g3/k_g3)**n_g3*(p_g5/k_g5)**n_g5 + a_g4_g4_g3_g5*(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3*(p_g5/k_g5)**n_g5 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g3/k_g3)**n_g3 +(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 +(p_g3/k_g3)**n_g3*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3*(p_g5/k_g5)**n_g5 ))-l_x_g4*x_g4
    dp_g4 = r_g4*x_g4- l_p_g4*p_g4
    dx_g5 = m_g5*(( alpha_g5 + a_g5_g4*(p_g4/k_g4)**n_g4 + a_g5_g5*(p_g5/k_g5)**n_g5 + a_g5_g3*(p_g3/k_g3)**n_g3 + a_g5_g4_g5*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 + a_g5_g4_g3*(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 + a_g5_g5_g3*(p_g5/k_g5)**n_g5*(p_g3/k_g3)**n_g3 + a_g5_g4_g5_g3*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g3/k_g3)**n_g3 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g5/k_g5)**n_g5 +(p_g3/k_g3)**n_g3 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 +(p_g5/k_g5)**n_g5*(p_g3/k_g3)**n_g3 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g3/k_g3)**n_g3 ))-l_x_g5*x_g5
    dp_g5 = r_g5*x_g5- l_p_g5*p_g5
    dx_g6 = m_g6*(( alpha_g6 + a_g6_g7*(p_g7/k_g7)**n_g7 + a_g6_g4*(p_g4/k_g4)**n_g4 + a_g6_g6*(p_g6/k_g6)**n_g6 + a_g6_g9*(p_g9/k_g9)**n_g9 + a_g6_g5*(p_g5/k_g5)**n_g5 + a_g6_g7_g4*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4 + a_g6_g7_g6*(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6 + a_g6_g7_g9*(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9 + a_g6_g7_g5*(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5 + a_g6_g4_g6*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 + a_g6_g4_g9*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9 + a_g6_g4_g5*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 + a_g6_g6_g9*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 + a_g6_g6_g5*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 + a_g6_g9_g5*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g6_g7_g4_g6*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 + a_g6_g7_g4_g9*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9 + a_g6_g7_g4_g5*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 + a_g6_g7_g6_g9*(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 + a_g6_g7_g6_g5*(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 + a_g6_g7_g9_g5*(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g6_g4_g6_g9*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 + a_g6_g4_g6_g5*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 + a_g6_g4_g9_g5*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g6_g6_g9_g5*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g6_g7_g4_g6_g9*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 + a_g6_g7_g4_g6_g5*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 + a_g6_g7_g4_g9_g5*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g6_g7_g6_g9_g5*(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g6_g4_g6_g9_g5*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g6_g7_g4_g6_g9_g5*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 )/( 1 +(p_g7/k_g7)**n_g7 +(p_g4/k_g4)**n_g4 +(p_g6/k_g6)**n_g6 +(p_g9/k_g9)**n_g9 +(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4 +(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6 +(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9 +(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 +(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 +(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 +(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 +(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 ))-l_x_g6*x_g6
    dp_g6 = r_g6*x_g6- l_p_g6*p_g6
    dx_g7 = m_g7*(( alpha_g7 + a_g7_g4*(p_g4/k_g4)**n_g4 + a_g7_g5*(p_g5/k_g5)**n_g5 + a_g7_g8*(p_g8/k_g8)**n_g8 + a_g7_g6*(p_g6/k_g6)**n_g6 + a_g7_g4_g5*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 + a_g7_g4_g8*(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8 + a_g7_g4_g6*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 + a_g7_g5_g8*(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8 + a_g7_g5_g6*(p_g5/k_g5)**n_g5*(p_g6/k_g6)**n_g6 + a_g7_g8_g6*(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 + a_g7_g4_g5_g8*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8 + a_g7_g4_g5_g6*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g6/k_g6)**n_g6 + a_g7_g4_g8_g6*(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 + a_g7_g5_g8_g6*(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 + a_g7_g4_g5_g8_g6*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g5/k_g5)**n_g5 +(p_g8/k_g8)**n_g8 +(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 +(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8 +(p_g5/k_g5)**n_g5*(p_g6/k_g6)**n_g6 +(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 +(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5*(p_g8/k_g8)**n_g8*(p_g6/k_g6)**n_g6 ))-l_x_g7*x_g7
    dp_g7 = r_g7*x_g7- l_p_g7*p_g7
    dx_g8 = m_g8*(( alpha_g8 + a_g8_g7*(p_g7/k_g7)**n_g7 + a_g8_g5*(p_g5/k_g5)**n_g5 + a_g8_g9*(p_g9/k_g9)**n_g9 + a_g8_g4*(p_g4/k_g4)**n_g4 + a_g8_g7_g5*(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5 + a_g8_g7_g9*(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9 + a_g8_g7_g4*(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4 + a_g8_g5_g9*(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9 + a_g8_g5_g4*(p_g5/k_g5)**n_g5*(p_g4/k_g4)**n_g4 + a_g8_g9_g4*(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 + a_g8_g7_g5_g9*(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9 + a_g8_g7_g5_g4*(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5*(p_g4/k_g4)**n_g4 + a_g8_g7_g9_g4*(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 + a_g8_g5_g9_g4*(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 + a_g8_g7_g5_g9_g4*(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 )/( 1 +(p_g7/k_g7)**n_g7 +(p_g5/k_g5)**n_g5 +(p_g9/k_g9)**n_g9 +(p_g4/k_g4)**n_g4 +(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5 +(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9 +(p_g7/k_g7)**n_g7*(p_g4/k_g4)**n_g4 +(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9 +(p_g5/k_g5)**n_g5*(p_g4/k_g4)**n_g4 +(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 +(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9 +(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5*(p_g4/k_g4)**n_g4 +(p_g7/k_g7)**n_g7*(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 +(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 +(p_g7/k_g7)**n_g7*(p_g5/k_g5)**n_g5*(p_g9/k_g9)**n_g9*(p_g4/k_g4)**n_g4 ))-l_x_g8*x_g8
    dp_g8 = r_g8*x_g8- l_p_g8*p_g8
    dx_g9 = m_g9*(( alpha_g9 + a_g9_g4*(p_g4/k_g4)**n_g4 + a_g9_g6*(p_g6/k_g6)**n_g6 + a_g9_g9*(p_g9/k_g9)**n_g9 + a_g9_g8*(p_g8/k_g8)**n_g8 + a_g9_g5*(p_g5/k_g5)**n_g5 + a_g9_g4_g6*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 + a_g9_g4_g9*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9 + a_g9_g4_g8*(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8 + a_g9_g4_g5*(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 + a_g9_g6_g9*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 + a_g9_g6_g8*(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8 + a_g9_g6_g5*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 + a_g9_g9_g8*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 + a_g9_g9_g5*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g9_g8_g5*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 + a_g9_g4_g6_g9*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 + a_g9_g4_g6_g8*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8 + a_g9_g4_g6_g5*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 + a_g9_g4_g9_g8*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 + a_g9_g4_g9_g5*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g9_g4_g8_g5*(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 + a_g9_g6_g9_g8*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 + a_g9_g6_g9_g5*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g9_g6_g8_g5*(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 + a_g9_g9_g8_g5*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 + a_g9_g4_g6_g9_g8*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 + a_g9_g4_g6_g9_g5*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 + a_g9_g4_g6_g8_g5*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 + a_g9_g4_g9_g8_g5*(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 + a_g9_g6_g9_g8_g5*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 + a_g9_g4_g6_g9_g8_g5*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g6/k_g6)**n_g6 +(p_g9/k_g9)**n_g9 +(p_g8/k_g8)**n_g8 +(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9 +(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8 +(p_g4/k_g4)**n_g4*(p_g5/k_g5)**n_g5 +(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 +(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8 +(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 +(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 +(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 +(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 +(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 +(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 +(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 +(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6*(p_g9/k_g9)**n_g9*(p_g8/k_g8)**n_g8*(p_g5/k_g5)**n_g5 ))-l_x_g9*x_g9
    dp_g9 = r_g9*x_g9- l_p_g9*p_g9
    dY = np.array([dx_g1,dp_g1,dx_g2,dp_g2,dx_g3,dp_g3,dx_g4,dp_g4,dx_g5,dp_g5,dx_g6,dp_g6,dx_g7,dp_g7,dx_g8,dp_g8,dx_g9,dp_g9,])
    return(dY)
#####################################################