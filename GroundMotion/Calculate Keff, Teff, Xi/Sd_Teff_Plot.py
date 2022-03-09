# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:25:35 2022

@author: VACALDER
"""

import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd

iCL = [0]

#plotting commands
linestyle_str = ['solid', 'dotted','dashed','dashdot','dashed',(0, (3, 5, 1, 5, 1, 5))] 
plt.rcParams.update({'font.family':'serif'})
plt.rcParams.update({'font.serif':'Times New Roman'})
plt.rcParams.update({'mathtext.fontset':'stix'})
colors = ['#d94545','#519872','#73a8d4','#f7b76d','#545066']

counter=-1

# Column Data

AxialLoad=192.0 #kip
Diameter = 48 #in
Height = Diameter*4
g=386.2197 # in/s^2
Keff = 12.08 # kip/in
AbsMaxDisp = 16.55 #in
Es = 29000 # ksi
YieldStress_Long = 60 #ksi
meff = AxialLoad / g

# Calculating Teff and ductility


Teff = (2 * math.pi) * (math.sqrt((meff / Keff)))

Lsp = 0.15 * YieldStress_Long * 0.875

e_steel_yield = YieldStress_Long/Es
phi_y = 2.25 * e_steel_yield/Diameter
delta_y = phi_y * (Height + Lsp) ** 2 / 3
delta_u = AbsMaxDisp
mu = abs(delta_u) / delta_y

SpectrumFile = open('RSN1492_CHICHI_TCU052-E.AT2.g3.csv')
SpectrumContent = SpectrumFile.readlines()
SDC = SpectrumContent[12:109]
SDC_cols = ['Period', 'SD', 'PSV', 'PSA']
SDC_Data = [line.split(',') for line in SDC[:]]
SDC_DF = pd.DataFrame(columns=SDC_cols, data=SDC_Data)
PeriodStringList = list(SDC_DF['Period'])
SpectralDisplacementStringList = list(SDC_DF['SD'])
PGD = float(SpectralDisplacementStringList[-1])
T = [float(i) for i in PeriodStringList]
SpectralDisplacementList = list(SDC_DF['SD'])
SD_Float = [float(i) for i in SpectralDisplacementList]
SD_at_Teff = np.interp(Teff, T, SD_Float)

if mu > 1:
    xi_eq = 0.05 + 0.565 * (mu - 1) / (mu * np.pi)
    DF = np.sqrt((0.07) / (0.05 + xi_eq))
    SD_Teff_xi_eq = DF * SD_at_Teff
    SD_Teff_xi_eq_all = list(np.array(SD_Float)*DF)
elif mu <= 1:
    SD_Teff_xi_eq = SD_at_Teff
    SD_Teff_xi_eq_all = list(np.array(SD_Float)*DF)
Teff_x = [Teff,Teff,0]
SD_y = [0,SD_at_Teff,SD_at_Teff]
SD_y_2 = [0,SD_Teff_xi_eq,SD_Teff_xi_eq]

plt.figure(1,figsize=(5,4))
plt.plot(T,SD_Float,color=colors[0],linestyle=linestyle_str[0])
plt.plot(T,SD_Teff_xi_eq_all,color=colors[2],linestyle=linestyle_str[0])
plt.plot(Teff_x,SD_y,color='black', linestyle='dashed' )
plt.plot(Teff_x,SD_y_2,color='black', linestyle='dashed' )
plt.annotate(r'$(T_{eff},SD(T_{eff},\xi=5\%))$', xy=(Teff, SD_at_Teff),
              xytext=(Teff-1.27, SD_at_Teff+5))
plt.annotate(r'$(T_{eff},SD(T_{eff},\xi_{eff}))$', xy=(Teff, SD_Teff_xi_eq),
              xytext=(Teff-1.27, SD_Teff_xi_eq+5))
plt.annotate(r'$\xi=5\%$', xy=(3, 115),fontweight='bold',color=colors[0],fontsize=11)
plt.annotate(r'$\xi_{eff}$', xy=(3, 62.5),fontweight='bold',color=colors[2],fontsize=11)

plt.scatter(Teff, SD_at_Teff, s=40, marker='^',color=colors[1])
plt.scatter(Teff, SD_Teff_xi_eq, s=40, marker='^',color=colors[1])

plt.xlim(0,4)
plt.ylim(0,130)
plt.xlabel('Period, T (s)', fontsize=11)
plt.ylabel('Spectral Displacement, SD (in)', fontsize=11)
plt.tick_params(direction='out',axis='both',labelsize=10)
plt.grid()
plt.show()

plt.savefig("SpectralDisplacement_SD(Teff,xi)_Calc.pdf",dpi=600,bbox_inches='tight', pad_inches=0)