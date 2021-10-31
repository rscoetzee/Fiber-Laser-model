# -*- coding: utf-8 -*-
"""
Created on Mon May 16 19:56:34 2016

@author: Turiya
"""

# RHS of Differential equations.
#import numpy as np

def RHS(h,PumpP,SigP,z,gam_p,gam_s,eta_p,eta_s,Nt,N2,N1,hplank,c,A_c,lam_p,lam_s,
        k,tau,Cp,Cs,index,delL,sig_Abs,sig_Ems): 
    
    kp = index*h*(gam_p*(sig_Ems*N2-sig_Abs*N1)*PumpP - eta_p*PumpP) 
    ks = index*h*(gam_s*(sig_Ems*N2-sig_Abs*N1)*SigP  - eta_s*SigP  + gam_s*sig_Ems*N2*((delL)*2*hplank*(c**2)/(lam_s**3)))
    
    return kp,ks