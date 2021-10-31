# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 13:26:25 2021

@author: Riaan
"""
"Yblaser"
import numpy as np
from RK4 import Rk4
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()

def yBlaser(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
            index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,
            Cs,Cp,A_clad,SinglePass,InputPump,FindPump,TotPout,TotSout):

    if SinglePass == 0:
        for j in range(0,len(InputPump)):

            P_in0 = np.zeros(len(WLengthspan))
            S_in0 = np.zeros(len(WLengthspan))
            P_in0[FindPump[0][0]] = InputPump[j]
    
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
                                                      index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,PIN[j],SIN,
                                                      P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad)
    
                Pout[i] = (1-R2p)*P[-1]
                Sout[i] = (1-R2s)*S[-1]
                P_in0    = R2p*SpecP
                S_in0    = R2s*SpecS
                SpecSOut = (1-R2s)*SpecS
                SpecPOut = (1-R2p)*SpecP
    
                index = 1
    
                # Backward propogation
                P,S,N2,N1,SpecS,SpecP = Rk4(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
                                                      index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,PIN[j],SIN,
                                                      P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad)
    
                P_in0 = SpecP*(R1p) + (1-R1p)*PIN
                S_in0 = SpecS*(R1s) + (1-R1s)*SIN
    
                index = 1

            completed1 = float(j/len(InputPump))
            completed2 = 'Progress = {0} %'.format(int(100*(completed1)))
            print(completed2)
    
            TotPout[j] = Pout[-1]
            TotSout[j] = Sout[-1]


        plt.figure(2,figsize=(25,20))
        plt.subplot(2, 2, 1)
        plt.plot(Pout,'b',linewidth=2,label='Pump')
        plt.plot(Sout,'r',linewidth=2,label='Signal')
        plt.legend(loc='upper left')
        plt.grid()
        plt.xlabel('Round trips')
        plt.ylabel('Power (W)')
    
        #plt.figure(3)
        plt.subplot(2, 2, 2)
        plt.plot(InputPump,TotPout,'bo',linewidth=2,label='Pump')
        plt.plot(InputPump,TotSout,'ro',linewidth=2,label='Signal')
        plt.grid()
        plt.xlabel('Input Pump Power (W)')
        plt.ylabel('Power (W)')
        ## slope fitting stuff
        b  = np.nonzero(TotSout > 0.002*InputPump[-1])[0][0]
        bb = np.nonzero(TotSout > 0.002*InputPump[-1])[0][-1]
        xx = np.polyfit(InputPump[b:bb+1],TotSout[b:bb+1],deg=1)[0]*InputPump[b:bb+1] + np.polyfit(InputPump[b:bb+1],TotSout[b:bb+1],deg=1)[1]
        plt.plot(InputPump[b:bb+1],xx,'g',label='Slope = {0}'.format(float(round(np.polyfit(InputPump[b:bb+1],TotSout[b:bb+1],deg=1)[0],3))))
    
        ##
        plt.legend(loc='upper left')
    
        #plt.figure(4)
        plt.subplot(2, 2, 3)
        plt.plot((1*10**+9)*WLengthspan,SpecSOut,'b',linewidth=1,label= 'Pump')
        plt.plot((1*10**+9)*WLengthspan,SpecPOut,'r',linewidth=1,label='Signal')
        plt.legend(loc='upper left')
        plt.grid()
        plt.xlabel('Wavelength')
        plt.ylabel('Power (W)')
    
        #plt.figure(5)
        plt.subplot(2, 2, 4)
        plt.plot(z,N2,linewidth=2,label='N2')
        plt.plot(z,N1,linewidth=2,label='N1')
        plt.legend(loc='upper left')
        plt.grid()
        plt.xlabel('Fiber length (m)')
        plt.ylabel('Population')
    
        print(datetime.now() - startTime)
        
        return