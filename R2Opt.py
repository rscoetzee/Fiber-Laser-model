# -*- coding: utf-8 -*-
"""
Created on Thu May 26 13:10:17 2016

@author: Riaan
"""
# script to calculate optimum outcoupling (R2s) for Signal

import numpy as np
from RK4 import Rk4

def R2sOpt(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
           index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,PIN,SIN,
           P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad,FindPump,CE):

           InputPump = 50 #CE*np.linspace(0.1,50,20)
           R2s = np.linspace(0,1,30)
           TotPout = np.zeros(len(R2s))
           TotSout = np.zeros(len(R2s))

           for j in range(0,len(R2s)):

                P_in0 = np.zeros(len(WLengthspan))
                S_in0 = np.zeros(len(WLengthspan))
                P_in0[FindPump[0][0]] = InputPump

                SIN = S_in0
                PIN = P_in0

                kp = np.zeros((len(WLengthspan),4))   # each row reprents prop of signal wavelength channel along fiber.
                ks = np.zeros((len(WLengthspan),4))   # a single coloum contains all wavelengths at one position in the fiber.

                P   = np.zeros(len(z))
                S   = np.zeros(len(z))
                N2  = np.zeros(len(z))
                N1  = np.zeros(len(z))

                PumpP = np.zeros((len(WLengthspan),int(L/h)))
                SigP  = np.zeros((len(WLengthspan),int(L/h)))

                for i in range(0,RT):

                    # Forward propogation
                    P,S,N2,N1,SpecS,SpecP = Rk4(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
                                                          index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s[j],R1s,Pout,Sout,PIN,SIN,
                                                          P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad)

                    Pout[i] = (1-R2p)*P[-1]
                    Sout[i] = (1-R2s[j])*S[-1]
                    P_in0    = R2p*SpecP
                    S_in0    = R2s[j]*SpecS

                    index = 1

                    # Backward propogation
                    P,S,N2,N1,SpecS,SpecP = Rk4(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
                                                          index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s[j],R1s,Pout,Sout,PIN,SIN,
                                                          P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad)

                    P_in0 = SpecP*(R1p) + (1-R1p)*PIN
                    S_in0 = SpecS*(R1s) + (1-R1s)*SIN

                    index = 1

                completed1 = float(j/len(R2s))
                completed2 = 'Progress = {0} %'.format(int(100*(completed1)))
                print(completed2)

                TotPout[j] = Pout[-1]
                TotSout[j] = Sout[-1]

           return TotPout,TotSout,R2s
