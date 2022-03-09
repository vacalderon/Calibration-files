# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:25:35 2022

@author: VACALDER
"""

import matplotlib.pyplot as plt
import numpy as np

iCL = [0, 5, 10, 15, 20]

#plotting commands
linestyle_str = ['solid', 'dotted','dashed','dashdot','dashed',(0, (3, 5, 1, 5, 1, 5))] 
plt.rcParams.update({'font.family':'serif'})
plt.rcParams.update({'font.serif':'Times New Roman'})
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
    
    
    plt.figure(1,figsize=(5,4))
    plt.plot(X,Y,color=colors[counter],linestyle=linestyle_str[counter])
    
    plt.xlim(-30,30)
    plt.ylim(-225,225)
    plt.xlabel('Diplacement (in)', fontsize=11)
    plt.ylabel('BaseShear (kip)', fontsize=11)
    plt.tick_params(direction='out',axis='both',labelsize=10)
    plt.xticks(np.linspace(-35,35,15))
    plt.legend(['CL= 0%','CL= 5%','CL= 10%','CL= 15%','CL= 20%'],edgecolor='black',loc='upper left',fontsize=9)
    plt.grid()
    plt.show()
    
    plt.savefig("Force_Diplacement_RSN1505.pdf",dpi=600,bbox_inches='tight', pad_inches=0)
    strain =  [line.split()[2] for line in lines_stress_strain]
    sigma =  [line.split()[1] for line in lines_stress_strain]
    
    strain_plt = [float(i) for i in strain]
    stress_plt = [float(i) for i in sigma]
    
    
    plt.figure(2,figsize=(5,4))
    plt.plot(strain_plt,stress_plt,color=colors[counter],linestyle=linestyle_str[counter])
    plt.xlabel('Strain (in/in)', fontsize=11)
    plt.ylabel('Stress (ksi)', fontsize=11)
    plt.tick_params(direction='out',axis='both',labelsize=10)
    plt.legend(['CL= 0%','CL= 5%','CL= 10%','CL= 15%','CL= 20%'],edgecolor='black',loc='lower right',fontsize=9)
    plt.grid()
    plt.show()
    
    plt.savefig("Stress_Strain_RSN1505.pdf",dpi=600,bbox_inches='tight', pad_inches=0)
    
    plt.figure(3,figsize=(5,4))
    displacements = [float(i) for i in x]
    
    original_list = strain_plt
    original_list_copy = original_list.copy()
    original_list_copy.reverse()
    reversed_list = original_list_copy.copy()
    
    
    plt.plot(X,strain_plt,color=colors[counter],linestyle=linestyle_str[counter])
    
    plt.xlim(-30,30)
    plt.ylim(-0.04,0.10)
    plt.xlabel('Dsplacement (in)', fontsize=11)
    plt.ylabel('Strain(in/in)', fontsize=11)
    plt.tick_params(direction='out',axis='both',labelsize=10)
    plt.xticks(np.linspace(-35,35,15))
    plt.legend(['CL= 0%','CL= 5%','CL= 10%','CL= 15%','CL= 20%'],edgecolor='black',loc='upper left',fontsize=9)
    plt.grid()
    plt.show()
    plt.savefig("Diplacement_Strain_RSN1505.pdf",dpi=600,bbox_inches='tight', pad_inches=0)