# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 13:56:27 2021

@author: Riaan
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 13:25:27 2021

@author: Riaan
"""
"Amplifier model, single pass, co-pump"

import numpy as np
from RK4counter import Rk4Counter
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()

def YbAmpBackward(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
                 index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,
                 Cs,Cp,A_clad,SinglePass,InputPump,FindPump,TotPout,TotSout):
    
        P_in0 = np.zeros(len(WLengthspan))
        S_in0 = np.zeros(len(WLengthspan))
        P_in0[FindPump[0][0]] = InputPump[-1]

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

        P,S,N2,N1,SpecS,SpecP,ASEpower = Rk4Counter(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
                                                    index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,PIN,SIN,
                                                    P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad)

        plt.figure(2,figsize=(15,6))
        plt.subplot(3,3,1)
        plt.plot(z,P,'b',label='Pump')
        plt.plot(z,S,'r',label='Signal')
        plt.plot(z,ASEpower,'k',label='ASE')
        plt.legend()
        plt.grid()
        plt.xlabel('Fiber Length (m)')
        plt.ylabel('Power (W)')

        plt.subplot(3,3,2)
        plt.plot(z,N2,'b',label='N2')
        plt.plot(z,N1,'g',label='N1')
        plt.legend()
        plt.grid()
        plt.xlabel('Fiber Length (m)')
        plt.ylabel('Population Number')

        # plt.subplot(3,3,3)
        # plt.plot((1*10**9)*WLengthspan,SpecS,'r',label='Signal')
        # plt.plot((1*10**9)*WLengthspan,SpecP,'b',label='Pump')
        # plt.legend()
        # plt.grid()
        # plt.xlabel('Wavelength (nm)')
        # plt.ylabel('Powers')
        
        # plt.subplot(3,3,4)
        # plt.plot(z,ASEpower,'k',label='ASE')
        # plt.xlabel('fiber length (nm)')
        # plt.ylabel('Power (W)')
        
        # print(ASEpower)
        

        print(datetime.now() - startTime)
        
        return