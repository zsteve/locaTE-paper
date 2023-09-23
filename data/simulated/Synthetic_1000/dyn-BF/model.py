#####################################################
import autograd.numpy as np
# This file is created automatically
def Model(Y,t,pars):
    # Parameters
    a_g1_g4 = pars[0]
    a_g1_g4_g6 = pars[1]
    a_g1_g6 = pars[2]
    a_g2_g1 = pars[3]
    a_g3_g2 = pars[4]
    a_g4_g3 = pars[5]
    a_g4_g3_g6 = pars[6]
    a_g4_g4 = pars[7]
    a_g4_g4_g3 = pars[8]
    a_g4_g4_g3_g6 = pars[9]
    a_g4_g4_g6 = pars[10]
    a_g4_g6 = pars[11]
    a_g6_g3 = pars[12]
    a_g6_g3_g6 = pars[13]
    a_g6_g4 = pars[14]
    a_g6_g4_g3 = pars[15]
    a_g6_g4_g3_g6 = pars[16]
    a_g6_g4_g6 = pars[17]
    a_g6_g6 = pars[18]
    a_g7_g4 = pars[19]
    a_g8_g4 = pars[20]
    alpha_g1 = pars[21]
    alpha_g2 = pars[22]
    alpha_g3 = pars[23]
    alpha_g4 = pars[24]
    alpha_g6 = pars[25]
    alpha_g7 = pars[26]
    alpha_g8 = pars[27]
    k_g1 = pars[28]
    k_g2 = pars[29]
    k_g3 = pars[30]
    k_g4 = pars[31]
    k_g6 = pars[32]
    k_g7 = pars[33]
    k_g8 = pars[34]
    l_p_g1 = pars[35]
    l_p_g2 = pars[36]
    l_p_g3 = pars[37]
    l_p_g4 = pars[38]
    l_p_g6 = pars[39]
    l_p_g7 = pars[40]
    l_p_g8 = pars[41]
    l_x_g1 = pars[42]
    l_x_g2 = pars[43]
    l_x_g3 = pars[44]
    l_x_g4 = pars[45]
    l_x_g6 = pars[46]
    l_x_g7 = pars[47]
    l_x_g8 = pars[48]
    m_g1 = pars[49]
    m_g2 = pars[50]
    m_g3 = pars[51]
    m_g4 = pars[52]
    m_g6 = pars[53]
    m_g7 = pars[54]
    m_g8 = pars[55]
    n_g1 = pars[56]
    n_g2 = pars[57]
    n_g3 = pars[58]
    n_g4 = pars[59]
    n_g6 = pars[60]
    n_g7 = pars[61]
    n_g8 = pars[62]
    r_g1 = pars[63]
    r_g2 = pars[64]
    r_g3 = pars[65]
    r_g4 = pars[66]
    r_g6 = pars[67]
    r_g7 = pars[68]
    r_g8 = pars[69]
    sigmaH_g1 = pars[70]
    sigmaH_g2 = pars[71]
    sigmaH_g3 = pars[72]
    sigmaH_g4 = pars[73]
    sigmaH_g6 = pars[74]
    sigmaH_g7 = pars[75]
    sigmaH_g8 = pars[76]
    # Variables
    x_g1 = Y[0]
    p_g1 = Y[1]
    x_g2 = Y[2]
    p_g2 = Y[3]
    x_g3 = Y[4]
    p_g3 = Y[5]
    x_g4 = Y[6]
    p_g4 = Y[7]
    x_g6 = Y[8]
    p_g6 = Y[9]
    x_g7 = Y[10]
    p_g7 = Y[11]
    x_g8 = Y[12]
    p_g8 = Y[13]
    dx_g1 = m_g1*(( alpha_g1 + a_g1_g4*(p_g4/k_g4)**n_g4 + a_g1_g6*(p_g6/k_g6)**n_g6 + a_g1_g4_g6*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 ))-l_x_g1*x_g1
    dp_g1 = r_g1*x_g1- l_p_g1*p_g1
    dx_g2 = m_g2*(( alpha_g2 + a_g2_g1*(p_g1/k_g1)**n_g1 )/( 1 +(p_g1/k_g1)**n_g1 ))-l_x_g2*x_g2
    dp_g2 = r_g2*x_g2- l_p_g2*p_g2
    dx_g3 = m_g3*(( alpha_g3 + a_g3_g2*(p_g2/k_g2)**n_g2 )/( 1 +(p_g2/k_g2)**n_g2 ))-l_x_g3*x_g3
    dp_g3 = r_g3*x_g3- l_p_g3*p_g3
    dx_g4 = m_g4*(( alpha_g4 + a_g4_g4*(p_g4/k_g4)**n_g4 + a_g4_g3*(p_g3/k_g3)**n_g3 + a_g4_g6*(p_g6/k_g6)**n_g6 + a_g4_g4_g3*(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 + a_g4_g4_g6*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 + a_g4_g3_g6*(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 + a_g4_g4_g3_g6*(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g3/k_g3)**n_g3 +(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 +(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 ))-l_x_g4*x_g4
    dp_g4 = r_g4*x_g4- l_p_g4*p_g4
    dx_g6 = m_g6*(( alpha_g6 + a_g6_g4*(p_g4/k_g4)**n_g4 + a_g6_g3*(p_g3/k_g3)**n_g3 + a_g6_g6*(p_g6/k_g6)**n_g6 + a_g6_g4_g3*(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 + a_g6_g4_g6*(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 + a_g6_g3_g6*(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 + a_g6_g4_g3_g6*(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 )/( 1 +(p_g4/k_g4)**n_g4 +(p_g3/k_g3)**n_g3 +(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3 +(p_g4/k_g4)**n_g4*(p_g6/k_g6)**n_g6 +(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 +(p_g4/k_g4)**n_g4*(p_g3/k_g3)**n_g3*(p_g6/k_g6)**n_g6 ))-l_x_g6*x_g6
    dp_g6 = r_g6*x_g6- l_p_g6*p_g6
    dx_g7 = m_g7*(( alpha_g7 + a_g7_g4*(p_g4/k_g4)**n_g4 )/( 1 +(p_g4/k_g4)**n_g4 ))-l_x_g7*x_g7
    dp_g7 = r_g7*x_g7- l_p_g7*p_g7
    dx_g8 = m_g8*(( alpha_g8 + a_g8_g4*(p_g4/k_g4)**n_g4 )/( 1 +(p_g4/k_g4)**n_g4 ))-l_x_g8*x_g8
    dp_g8 = r_g8*x_g8- l_p_g8*p_g8
    dY = np.array([dx_g1,dp_g1,dx_g2,dp_g2,dx_g3,dp_g3,dx_g4,dp_g4,dx_g6,dp_g6,dx_g7,dp_g7,dx_g8,dp_g8,])
    return(dY)
#####################################################