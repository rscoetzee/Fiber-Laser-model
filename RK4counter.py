# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 14:05:18 2021

@author: Riaan
"""
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 19:49:36 2016

@author: Turiya
"""
# Rk4 function for Laser model

import numpy as np #,flipud

def Rk4Counter(h,z,gam_p,gam_s,eta_p,eta_s,P_in0,S_in0,Nt,hplank,c,A_c,lam_p,lam_s,tau,
        index,delL,sig_Ems,sig_Abs,L,WLengthspan,RT,R2p,R1p,R2s,R1s,Pout,Sout,PIN,SIN,
        P,S,N2,N1,PumpP,SigP,kp,ks,Cs,Cp,A_clad):

        P[0] = np.sum(P_in0)
        S[0] = np.sum(S_in0)

        PumpP[:,0] = P_in0
        SigP[:,0]  = S_in0


        N2num   = (Cp)*np.sum(P_in0*lam_p*sig_Abs)*Nt        + (Cs)*np.sum(S_in0*lam_s*sig_Abs)*Nt
        N2denom = (Cp)*np.sum(P_in0*lam_p*(sig_Abs+sig_Ems)) + (Cs)*np.sum(S_in0*lam_s*(sig_Abs+sig_Ems)) + (1/tau)

        ASEpower  = np.zeros(len(z))
        
        N2[0] = N2num/N2denom
        N1[0] = Nt-N2[0]
        
        SpecS = 0
        SpecP = 0
        
        print(P[-1])
        
        while np.absolute(np.sum(P_in0)-P[-1]) > 0.00001:
            
            corr = np.sign(np.sum(P_in0)-P[-1])*10**(np.floor(np.log10(np.absolute(np.sum(P_in0)-P[-1]))))
            PumpP[:,0] = P_in0[12504] + corr
            print(corr)
                
            for k in range(0,len(z)-1):
    
                ASE = (gam_s)*sig_Ems*N2[k]*((delL)*2*hplank*(c**2)/(lam_s**3))
    
                kp[:,0] = -index*h*(gam_p*(sig_Ems*N2[k]-sig_Abs*N1[k])*PumpP[:,k] - eta_p*PumpP[:,k])
                ks[:,0] = index*h*(gam_s*(sig_Ems*N2[k]-sig_Abs*N1[k])*SigP[:,k]  - eta_s*SigP[:,k] + ASE)
    
                kp[:,1] = -index*h*(gam_p*(sig_Ems*N2[k]-sig_Abs*N1[k])*(PumpP[:,k]+kp[:,0]) - eta_p*(PumpP[:,k]+kp[:,0]))
                ks[:,1] = index*h*(gam_s*(sig_Ems*N2[k]-sig_Abs*N1[k])*(SigP[:,k]+ks[:,0])  - eta_s*(SigP[:,k]+ks[:,0]) + ASE)
    
                kp[:,2] = -index*h*(gam_p*(sig_Ems*N2[k]-sig_Abs*N1[k])*(PumpP[:,k]+kp[:,1]) - eta_p*(PumpP[:,k]+kp[:,1]))
                ks[:,2] = index*h*(gam_s*(sig_Ems*N2[k]-sig_Abs*N1[k])*(SigP[:,k]+ks[:,1])  - eta_s*(SigP[:,k]+ks[:,1]) + ASE)
    
                kp[:,3] = -index*h*(gam_p*(sig_Ems*N2[k]-sig_Abs*N1[k])*(PumpP[:,k]+kp[:,2]) - eta_p*(PumpP[:,k]+kp[:,2]))
                ks[:,3] = index*h*(gam_s*(sig_Ems*N2[k]-sig_Abs*N1[k])*(SigP[:,k]+ks[:,2])  - eta_s*(SigP[:,k]+ks[:,2]) + ASE)
    
                PumpP[:,k+1] = PumpP[:,k] + (kp[:,0] + 2*kp[:,1] + 2*kp[:,2] + kp[:,3])/6
                SigP[:,k+1]  = SigP[:,k]  + (ks[:,0] + 2*ks[:,1] + 2*ks[:,2] + ks[:,3])/6            
                
                P[k+1] = np.sum(PumpP[:,k+1])
                S[k+1] = np.sum(SigP[:,k+1])
    
                N2num   = (Cp)*(np.sum(PumpP[:,k+1]*lam_p*sig_Abs)*Nt)        + (Cs)*(np.sum(SigP[:,k+1]*lam_s*sig_Abs)*Nt)
                N2denom = (Cp)*(np.sum(PumpP[:,k+1]*lam_p*(sig_Abs+sig_Ems))) + (Cs)*(np.sum(SigP[:,k+1]*lam_s*(sig_Abs+sig_Ems))) + (1/tau)
    
                N2[k+1] = N2num/N2denom
                N1[k+1] = Nt - N2[k+1]
                
                ASEpower[k] = np.sum(ASE)
    
            SpecS = SigP[:,k+1]
            SpecP = PumpP[:,k+1]
            
            print(P[-1])
             
        return P,S,N2,N1,SpecS,SpecP,ASEpower
    
#                kp[:,0],ks[:,0] = RHS(h,PumpP[:,k],SigP[:,k],z[k],gam_p,gam_s,eta_p,eta_s,                          # k1p, k1s
#                                  Nt,N2[k],N1[k],hplank,c,A_c,lam_p,lam_s,k,tau,Cp,Cs,index,delL,sig_Abs,sig_Ems)
#
#                kp[:,1],ks[:,1] = RHS(h,PumpP[:,k]+kp[:,0]/2,SigP[:,k]+ks[:,0]/2,z[k]+h/2,gam_p,gam_s,eta_p,eta_s,  # k2p, k2s
#                                  Nt,N2[k],N1[k],hplank,c,A_c,lam_p,lam_s,k,tau,Cp,Cs,index,delL,sig_Abs,sig_Ems)
#
#                kp[:,2],ks[:,2] = RHS(h,PumpP[:,k]+kp[:,1]/2,SigP[:,k]+ks[:,1]/2,z[k]+h/2,gam_p,gam_s,eta_p,eta_s,  # k3p, k3s
#                                  Nt,N2[k],N1[k],hplank,c,A_c,lam_p,lam_s,k,tau,Cp,Cs,index,delL,sig_Abs,sig_Ems)
#
#                kp[:,3],ks[:,3] = RHS(h,PumpP[:,k]+kp[:,2],SigP[:,k]+ks[:,2],z[k]+h,gam_p,gam_s,eta_p,eta_s,        # k4p, k4s
#                                  Nt,N2[k],N1[k],hplank,c,A_c,lam_p,lam_s,k,tau,Cp,Cs,index,delL,sig_Abs,sig_Ems)

#    kp = index*h*(gam_p*(sig_Ems*N2-sig_Abs*N1)*PumpP - eta_p*PumpP)
#    ks = index*h*(gam_s*(sig_Ems*N2-sig_Abs*N1)*SigP  - eta_s*SigP  + gam_s*sig_Ems*N2*((delL)*2*hplank*(c**2)/(lam_s**3)))
