#####################################################
import autograd.numpy as np
# This file is created automatically
def Model(Y,t,pars):
    # Parameters
    a_g1_g7 = pars[0]
    a_g2_g1 = pars[1]
    a_g3_g2 = pars[2]
    a_g4_g3 = pars[3]
    a_g5_g4 = pars[4]
    a_g6_g5 = pars[5]
    a_g7_g6 = pars[6]
    a_g7_g7 = pars[7]
    a_g7_g7_g6 = pars[8]
    alpha_g1 = pars[9]
    alpha_g2 = pars[10]
    alpha_g3 = pars[11]
    alpha_g4 = pars[12]
    alpha_g5 = pars[13]
    alpha_g6 = pars[14]
    alpha_g7 = pars[15]
    k_g1 = pars[16]
    k_g2 = pars[17]
    k_g3 = pars[18]
    k_g4 = pars[19]
    k_g5 = pars[20]
    k_g6 = pars[21]
    k_g7 = pars[22]
    l_p_g1 = pars[23]
    l_p_g2 = pars[24]
    l_p_g3 = pars[25]
    l_p_g4 = pars[26]
    l_p_g5 = pars[27]
    l_p_g6 = pars[28]
    l_p_g7 = pars[29]
    l_x_g1 = pars[30]
    l_x_g2 = pars[31]
    l_x_g3 = pars[32]
    l_x_g4 = pars[33]
    l_x_g5 = pars[34]
    l_x_g6 = pars[35]
    l_x_g7 = pars[36]
    m_g1 = pars[37]
    m_g2 = pars[38]
    m_g3 = pars[39]
    m_g4 = pars[40]
    m_g5 = pars[41]
    m_g6 = pars[42]
    m_g7 = pars[43]
    n_g1 = pars[44]
    n_g2 = pars[45]
    n_g3 = pars[46]
    n_g4 = pars[47]
    n_g5 = pars[48]
    n_g6 = pars[49]
    n_g7 = pars[50]
    r_g1 = pars[51]
    r_g2 = pars[52]
    r_g3 = pars[53]
    r_g4 = pars[54]
    r_g5 = pars[55]
    r_g6 = pars[56]
    r_g7 = pars[57]
    sigmaH_g1 = pars[58]
    sigmaH_g2 = pars[59]
    sigmaH_g3 = pars[60]
    sigmaH_g4 = pars[61]
    sigmaH_g5 = pars[62]
    sigmaH_g6 = pars[63]
    sigmaH_g7 = pars[64]
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
    dx_g1 = m_g1*(( alpha_g1 + a_g1_g7*(p_g7/k_g7)**n_g7 )/( 1 +(p_g7/k_g7)**n_g7 ))-l_x_g1*x_g1
    dp_g1 = r_g1*x_g1- l_p_g1*p_g1
    dx_g2 = m_g2*(( alpha_g2 + a_g2_g1*(p_g1/k_g1)**n_g1 )/( 1 +(p_g1/k_g1)**n_g1 ))-l_x_g2*x_g2
    dp_g2 = r_g2*x_g2- l_p_g2*p_g2
    dx_g3 = m_g3*(( alpha_g3 + a_g3_g2*(p_g2/k_g2)**n_g2 )/( 1 +(p_g2/k_g2)**n_g2 ))-l_x_g3*x_g3
    dp_g3 = r_g3*x_g3- l_p_g3*p_g3
    dx_g4 = m_g4*(( alpha_g4 + a_g4_g3*(p_g3/k_g3)**n_g3 )/( 1 +(p_g3/k_g3)**n_g3 ))-l_x_g4*x_g4
    dp_g4 = r_g4*x_g4- l_p_g4*p_g4
    dx_g5 = m_g5*(( alpha_g5 + a_g5_g4*(p_g4/k_g4)**n_g4 )/( 1 +(p_g4/k_g4)**n_g4 ))-l_x_g5*x_g5
    dp_g5 = r_g5*x_g5- l_p_g5*p_g5
    dx_g6 = m_g6*(( alpha_g6 + a_g6_g5*(p_g5/k_g5)**n_g5 )/( 1 +(p_g5/k_g5)**n_g5 ))-l_x_g6*x_g6
    dp_g6 = r_g6*x_g6- l_p_g6*p_g6
    dx_g7 = m_g7*(( alpha_g7 + a_g7_g7*(p_g7/k_g7)**n_g7 + a_g7_g6*(p_g6/k_g6)**n_g6 + a_g7_g7_g6*(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6 )/( 1 +(p_g7/k_g7)**n_g7 +(p_g6/k_g6)**n_g6 +(p_g7/k_g7)**n_g7*(p_g6/k_g6)**n_g6 ))-l_x_g7*x_g7
    dp_g7 = r_g7*x_g7- l_p_g7*p_g7
    dY = np.array([dx_g1,dp_g1,dx_g2,dp_g2,dx_g3,dp_g3,dx_g4,dp_g4,dx_g5,dp_g5,dx_g6,dp_g6,dx_g7,dp_g7,])
    return(dY)
#####################################################