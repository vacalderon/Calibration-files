# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:25:35 2022

@author: VACALDER
"""

import matplotlib.pyplot as plt
import numpy as np

iCL = [0]

#plotting commands
linestyle_str = ['solid', 'dotted','dashed','dashdot','dashed',(0, (3, 5, 1, 5, 1, 5))] 
plt.rcParams.update({'font.family':'serif'})
plt.rcParams.update({'font.serif':'Times New Roman'})
plt.rcParams.update({'mathtext.fontset':'stix'})
colors = ['#d94545','#519872','#73a8d4','#f7b76d','#545066']

counter=-1

# data from analysis
for i in iCL:
    counter+=1
    
    d = open('CL_'+str(i)+'_DFree.out')
    F = open('CL_'+str(i)+'_RBase.out')
    stress_strain_file = open('CL_'+str(i)+'_StressStrain4.out')
    
    lines_stress_strain = stress_strain_file.readlines()
    linesd = d.readlines()
    linesf = F.readlines()
    
    x = [line.split()[1] for line in linesd]
    y = [line.split()[1] for line in linesf]
    
    
    
    
    
    X=[-1*float(i) for i in x]
    Y=[float(i) for i in y]
    
    maxDisp = max(X)
    minDisp = min(X)
    if maxDisp > abs(minDisp):
        AbsMaxDisp = maxDisp
    elif maxDisp < abs(minDisp):
        AbsMaxDisp = minDisp
    maxDispPoss = X.index(AbsMaxDisp)
    maxForce_at_maxDisp = Y[maxDispPoss]
    Keff = abs(maxForce_at_maxDisp) / abs(AbsMaxDisp)
    
    keff_x = [0, AbsMaxDisp, AbsMaxDisp]
    keff_y = [0, maxForce_at_maxDisp,0]
    
    plt.figure(1,figsize=(5,4))
    plt.plot(X,Y,color=colors[counter],linestyle=linestyle_str[counter])
    plt.plot(keff_x,keff_y,color='black', linestyle='dashed' )
    plt.axhline(y = 0.0, color = 'black', linestyle = '-',linewidth=0.5)
    plt.axvline(x = 0.0, color = 'black', linestyle = '-',linewidth=0.5)
    plt.xlim(-30,30)
    plt.ylim(-250,250)
    plt.xlabel('Diplacement (in)', fontsize=11)
    plt.ylabel('BaseShear (kip)', fontsize=11)
    plt.tick_params(direction='out',axis='both',labelsize=10)
    plt.annotate(r'$(d_{max},F(d_{max}))$', xy=(AbsMaxDisp, maxForce_at_maxDisp),
                 xytext=(AbsMaxDisp+0.5, maxForce_at_maxDisp+10))
    plt.text(AbsMaxDisp+1, 120, r'$k_{eff}=\frac{F(d_{max})}{d_{max}}$', fontsize=11,
             bbox={'facecolor': 'white', 'alpha': 1.0, 'pad': 3})

    #plt.xticks(np.linspace(-20,20,9))
    plt.grid()
    plt.show()
    
    plt.savefig("Force_Diplacement_Keff_Calc.pdf",dpi=600,bbox_inches='tight', pad_inches=0)
    