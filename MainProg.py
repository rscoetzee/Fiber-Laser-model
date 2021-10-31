# -*- coding: utf-8 -*-

# Fiber laser model, based on ASE channel method.
# Solve with RK4 Cash-Karp scheme. CW model for now. Yb fiber laser.
# Solved in steady-state (dN2/dt = 0)
# RS Coetzee.

from datetime import datetime
startTime = datetime.now()

import numpy as np
import matplotlib.pyplot as plt
from YbParameter import Ybdata
from RK4 import Rk4

plt.close('all')

delL,WLengthspan,sig_Abs,sig_Ems,c,e0,hplank,hbar,sig_Pe,sig_Pa,sig_Sa,sig_Se,tau,Nt,L,CE,P_in0,S_in0,InputPump,InputSig,h,PLoss,SLoss,eta_p,eta_s,r_c,r_clad,A_c,A_clad,R1p,R2p,R1s,R2s,z,gam_p,gam_s = Ybdata()


plt.figure(1,figsize=(10,10))
plt.plot((1*10**9)*WLengthspan,sig_Ems,'r', linewidth = 2,label='emission', antialiased = True)
plt.plot((1*10**9)*WLengthspan,sig_Abs,'b', linewidth = 2,label='absorption',antialiased = True)
plt.legend()
plt.grid()
plt.xlabel('Wavelength (nm)')
plt.ylabel('Cross section (cm^-1)')
#plt.axis((1*10**9)*Wavelengths[0],(1*10**9)*Wavelengths[-1],0,np.max(sig_Ems))

Xabs = np.zeros((len(WLengthspan),len(z)))
Xems = np.zeros((len(WLengthspan),len(z)))
Xabs[:,0] = sig_Abs
Xems[:,0] = sig_Ems

X = np.zeros(len(WLengthspan))
lam_p = WLengthspan
lam_s = WLengthspan

P_in0 = np.zeros(len(WLengthspan))
FindPump = np.nonzero(WLengthspan>975.99*1*10**-9)
P_in0[FindPump[0][0]] = InputPump

S_in0 = np.zeros(len(WLengthspan))
FindSeed = np.nonzero(WLengthspan>1063.99*1*10**-9)
S_in0[FindSeed[0][0]] = InputSig

SinglePass = 0
CounterPump = 0
OptCoupling = 1
VarySeedSpec = 0

RT = 10
Pout = np.zeros(RT)
Sout = np.zeros(RT)
index = 1

InputPump = CE*np.linspace(0.1,100,20)
TotPout = np.zeros(len(InputPump))
TotSout = np.zeros(len(InputPump))

Cp = (gam_p/(hplank*c*A_c))
Cs = (gam_s/(hplank*c*A_c))

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
            P,S,N2,N1,SpecS,SpecP,ASEpower = Rk4(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
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
            P,S,N2,N1,SpecS,SpecP,ASEpower = Rk4(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
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

if SinglePass == 1:
    if CounterPump == 0:
        # from ampForward import YbAmpForward
        # YbAmpForward(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
        #              index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,
        #              Cs,Cp,A_clad,SinglePass,InputPump,InputSig,FindPump,FindSeed,TotPout,TotSout)
        
        P_in0 = np.zeros(len(WLengthspan))
        S_in0 = np.zeros(len(WLengthspan))
        P_in0[FindPump[0][0]] = InputPump[-1]
        S_in0[FindSeed[0][0]] = InputSig

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

        P,S,N2,N1,SpecS,SpecP,ASEpower = Rk4(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
                                             index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,PIN,SIN,
                                             P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad)
        
        plt.figure(2,figsize=(15,6))
        plt.subplot(1,3,1)
        plt.plot(z,P,'b',label='Pump')
        plt.plot(z,S,'r',label='Signal')
        # plt.plot(z,ASEpower,'k',label='ASE')
       # plt.legend()
        plt.grid()
        plt.xlabel('Fiber Length (m)')
        plt.ylabel('Power (W)')
        plt.legend()

        plt.subplot(1,3,2)
        plt.plot(z,N2,'b',label='N2')
        plt.plot(z,N1,'g',label='N1')
       # plt.legend()
        plt.grid()
        plt.xlabel('Fiber Length (m)')
        plt.ylabel('Population Number')
        plt.legend()

        plt.subplot(1,3,3)
        plt.plot((1*10**9)*WLengthspan,SpecS,'r',label='Signal')
        plt.plot((1*10**9)*WLengthspan,SpecP,'b',label='Pump')
       # plt.legend()
        plt.grid()
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Powers')
        plt.legend()
        
        # plt.subplot(3,3,4)
        # plt.plot(z,ASEpower,'k',label='ASE')
        # plt.xlabel('fiber length (nm)')
        # plt.ylabel('Power (W)')
        # print(datetime.now() - startTime)
        
        # return SpecS,ASEpower,P,S,N2,N1
        
    #     if VarySeedSpec == 1:
            
    #         #plt.close('all')
            
    #         from ampForward import YbAmpForward
    #         InputSigPowers = CE*np.linspace(0.00001,5,10)
                
    #         for k in range(0,len(InputSigPowers)-1):
    #             print(k)
    #             InputSig = InputSigPowers[k]
    #             SpecS,ASEpower,P,S,N2,N1 =  YbAmpForward(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
    #                                         index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,
    #                                         Cs,Cp,A_clad,SinglePass,InputPump,InputSig,FindPump,FindSeed,TotPout,TotSout)
                
    #             plt.figure(4,figsize=(15,10))
    #             plt.plot(WLengthspan*10**9,10*np.log10(SpecS))
         
        
    # elif CounterPump == 1:
    #     from ampBackward import YbAmpBackward
    #     YbAmpBackward(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
    #                  index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,
    #                  Cs,Cp,A_clad,SinglePass,InputPump,FindPump,TotPout,TotSout)
              

# if OptCoupling == 0:

#     from R2Opt import R2sOpt

#     TotPout,TotSout,R2s = R2sOpt(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
#                                  index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,PIN,SIN,
#                                  P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad,FindPump,CE)

#     plt.figure(3,figsize=(15,10))
# #    plt.subplot(2, 2, 1)
#     plt.plot(R2s,TotPout,'bo',linewidth=2,label='Pump')
#     plt.plot(R2s,TotSout,'ro',linewidth=2,label='Signal')
#     plt.legend(loc='upper left')
#     plt.grid()
#     plt.xlabel('R2s Reflectivity')
#     plt.ylabel('Power (W)')
#     plt.title('Calculated Output coupling for constant Pump.')

    print(datetime.now() - startTime)

# need to normalize input pump according to a given pump bandwidth.
#DelP = 3*10**-9
#DelP = c/DelP                              # this was for broadband pump, but keep single freq for now
#PumpCenter = 1030*10**-9
#PumpCenter = c/PumpCenter
#start = c/(1029.999*10**-9)
#stop  = c/(1030.001*10**-9)
#Pumprange = np.linspace(start,stop,1000)
#print(np.shape(Pumprange))
#G = (2/DelP)*np.sqrt(np.log(2)/np.pi)*np.exp(-4*((Pumprange-PumpCenter)**2)*(np.log(2)/(np.sqrt(DelP)**2)))
#plt.figure(2)
#plt.plot(Pumprange,G)
