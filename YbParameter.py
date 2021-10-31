# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 12:52:02 2021

@author: Riaan
"""
"Parameter list"

import numpy as np

def Ybdata():

    #  parameter list
    c       = 2.99792458*10**8
    e0      = 8.854*10**-12
    hplank  = 6.62606957*10**-34
    hbar    = hplank/2*np.pi
    sig_Pe  = 2.4*10**-24
    sig_Pa  = 2.4254*10**-24     # /1m**2
    sig_Sa  = 0.1*10**-24
    sig_Se  = 0.7*10**-24
    tau     = 840*10**-6             # Spont. life time tau2 (s)
    Nt      = 6.3*10**25         # Doping concentration
    L       = 3            # Fiber length (m)
    CE      = 1              # Coupling efficiency
    P_in0   = 75             # input pump power (W)
    S_in0   = 1
    P_in0   = CE*P_in0
    S_in0   = CE*S_in0
    InputPump = CE*75
    InputSig  = CE*2
    #h       = L*0.01             # step size (m)
    h       = 0.01                 # 10 mm step size
    PLoss   = (5.6)   # dB/m
    SLoss   = (3.6/1000)     # dB/m
    eta_p   = PLoss*h
    eta_s   = SLoss*h                # Fiber signal and pump losses
    r_c     = 7.5*10**-6             # Fiber core radius (m)
    r_clad  = 65*10**-6
    A_c     = np.pi*r_c**2
    A_clad  = np.pi*r_clad**2
    
    R1p     = 0.0001
    R2p     = 0.0001                  # Reflectivities for pump and signal
    R1s     = 0.0336870
    R2s     = 0.0336870
    
    z = np.arange(0,L,h)                          # z vectors to track evolution of Pump and Signal along fiber
    #z =
    gam_p    = A_c/A_clad             # Pump and signal overlap factors with core
    gam_s    = 0.97
    
    # setting up abs cross section data for wavelengths
    Data = np.genfromtxt('C:\\Users\\Riaan\\Dropbox\\Programming\\Fiber laser model\\abs_cross5.csv',delimiter=',')
    Wavelengths1 = Data[:,0]*10**-9
    Cross_sections_Abs_Data = Data[:,1]*10**-25
    # setting up ems cross section data for wavelengths
    Data = np.genfromtxt('C:\\Users\\Riaan\\Dropbox\\Programming\\Fiber laser model\\em_cross1.csv',delimiter=',')
    Wavelengths2 = Data[:,0]*10**-9
    Cross_sections_Ems_Data = Data[:,1]*10**-25
    delL = 0.1*10**-9             #  Wavelength spacing/channels
    WLengthspan = np.arange(Wavelengths1[0],Wavelengths2[-1],delL)
    sig_Abs = np.interp(WLengthspan,Wavelengths1,Cross_sections_Abs_Data)
    sig_Ems = np.interp(WLengthspan,Wavelengths2,Cross_sections_Ems_Data)
    

    return delL,WLengthspan,sig_Abs,sig_Ems,c,e0,hplank,hbar,sig_Pe,sig_Pa,sig_Sa,sig_Se,tau,Nt,L,CE,P_in0,S_in0,InputPump,InputSig,h,PLoss,SLoss,eta_p,eta_s,r_c,r_clad,A_c,A_clad,R1p,R2p,R1s,R2s,z,gam_p,gam_s
            