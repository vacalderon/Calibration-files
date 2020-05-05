# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:12:06 2019

@author: pchi893
"""

from openseespy.opensees import *
# import the os module
#    import os
import math
import numpy as np
import matplotlib.pyplot as plt
from LibUnitsMUS import *
import ManderCC_Rect
from __main__ import *



# -----------------------------------------------------------------------------
#	OpenSees (Tcl) code by:	Silvia Mazzoni & Frank McKenna, 2006

#
#    ^Y
#    |
#    3       __ 
#    |          | 
#    |          |
#    |          |
#  (2)       LCol
#    |          |
#    |          |
#    |          |
#  =2=      _|_  -------->X
#  =1=      ZeroLength

# ------------------------------------------------------------------------------
# |                      IMPORTS
# ------------------------------------------------------------------------------
from typing import Any

wipe()

# ------------------------------------------------------------------------------
#                           GENERATE GEOMETRY
# ------------------------------------------------------------------------------


model('basic', '-ndm', 2, '-ndf', 3)


# define section geometry
HCol =  300 * mm  # Column Height
BCol =  300 * mm  # Column Base
LCol = 1500 * mm  # column length
ACol = HCol * BCol  # cross-sectional area, make stiff
IzCol = 1/12 * BCol * HCol ** 3  # Column moment of inertia
fpc28 = 20*MPa  # Strength of the concrete
ALR = 0.22        # Axial Load Ratio
PCol = ALR*ACol*fpc28  # Axial Load
Weight = PCol     # superstructure weight

#PCol = Weight  # nodal dead-load weight per column
Mass = PCol / g

node(1, 0.0, 0.0)
node(2, 0.0, 0.0)
node(3, 0.0, LCol)
#  Node, Dx, Dy, Rz
fix(1, 1, 1, 1)
fix(2, 1, 0, 0)

mass(3, Mass, 1e-9, 0.0)



# MATERIAL parameters
IDconcC = 1  # material ID tag -- confined cover concrete
IDconcU = 2  # material ID tag -- unconfined cover concrete
IDreinf = 3  # material ID tag -- reinforcement
IDSP    = 4  # material ID tag -- Strain Penetration
# Define materials for nonlinear columns
# ------------------------------------------
# Longitudinal steel properties

CL = 20.0

Fy = 520.0 * MPa * (1-0.021*CL) # STEEL yield stress
Fu = 620.0 * MPa * (1-0.021*CL) # Steel Ultimate Stress
Es = 290000.0 * MPa  # modulus of steel
Bs = 0.01  # strain-hardening ratio
R0 = 20.0  # control the transition from elastic to plastic branches
cR1 = 0.925  # control the transition from elastic to plastic branches
cR2 = 0.15  # control the transition from elastic to plastic branches
c = 30 * mm  # Column cover to reinforcing steel NA.
numBarsSec = 2  # number of uniformly-distributed longitudinal-reinforcement bars Extreme Fibers
#numBarsMid= 2   # number of uniformly-distributed longitudinal-reinforcement bars Middle Portion
barAreaSec = 200 * mm2 *(1-CL/100)  # area of longitudinal-reinforcement bars
dbl = (barAreaSec*4/math.pi)**0.5
 

# Transverse Steel Properties
fyt = 520 * MPa  # Yield Stress of Transverse Steel
Ast = 50.3 * mm2  # Area of transverse steel
dbt = 8 * mm  # Diameter of transverse steel
st = 300 * mm  # Spacing of spiral
nbx =2
nby =2
bx=BCol-2*c
by=HCol-2*c
Ashx = nbx*Ast
Ashy = nby*Ast
roh_x = 2*Ashx/(st*bx)
roh_y = 2*Ashy/(st*by)

#Dprime = DCol - 2 * c - dbt*0.5  # Inner core diameter
# print('Dprime= ',Dprime)
#Rbl = Dprime * 0.5 - dbt*0.5 - dbl * 0.5  # Location of longitudinal bar
# print('Rbl =', Rbl)

# nominal concrete compressive strength
fpc = 32.4 * MPa  # CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
Ec = 57.0 * ksi * math.sqrt(fpc / psi)  # Concrete Elastic Modulus

# unconfined concrete
fc1U    = -fpc;  # UNCONFINED concrete (todeschini parabolic model), maximum stress
eps1U   = -0.002  # strain at maximum strength of unconfined concrete
fc2U    = 0.0 * MPa  # ultimate stress
eps2U   = -0.0064  # strain at ultimate stress
lambdac = 0.1  # ratio between unloading slope at $eps2 and initial slope $Ec

mand = ManderCC_Rect.ManderCC_Rect(fpc, fyt, roh_x, roh_y)

fc   = mand[0]
eps1 = mand[1]
fc2  = mand[2]
eps2 = mand[3]

# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
uniaxialMaterial('Concrete01', IDconcC, fc, eps1, fc2, eps2)

# Cover concrete (unconfined)
uniaxialMaterial('Concrete01', IDconcU, fc1U, eps1U, fc2U, eps2U)

# STEEL
# Reinforcing steel 
params = [R0, cR1, cR2]
#                        tag  fy E0    b
uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs, R0, cR1, cR2)


# STRAIN PENETRATION MATERIAL
SPalpha=0.4
SPsy=0.1*((dbl*Fy)*(2*SPalpha+1)/(4000*((-fc)**0.5)))**(1/SPalpha)+0.013
SPsu=35*SPsy
SPb=0.45
SPR=1.01


#uniaxialMaterial StrPen01   Tag  fy  sy fu su b  R  
uniaxialMaterial('Bond_SP01',IDSP, Fy,SPsy,Fu,SPsu,SPb,SPR)

# Writing Material data to file
with open("mat.out", 'w') as matfile:
    matfile.write("%s %s %s %s %s %s %s %s %s %s %s %s \n" %(Fy, fyt, Ast, st, PCol, HCol, barAreaSec, fc,SPsy,SPsu,SPb,SPR))
matfile.close

#-------------------------------------------------------------------------
#               DEFINE PLASTICE HIGE PROPERTIES
#-------------------------------------------------------------------------    

k=0.2*(Fu/Fy - 1)
if k > 0.08:
    k=0.08
Leff=LCol
Lpc=k*Leff + 0.4*HCol
gamma=0.33 #Assuming unidirectional action
Lpt=Lpc+gamma*HCol
# FIBER SECTION properties -------------------------------------------------------------
# Define cross-section for nonlinear columns
# --------------------------------------------------------------------------------------
#Create Rectangular Sections
# FIBER SECTION properties -------------------------------------------------------------
# symmetric section
#
#
#
#             ------------------------------  --
#             |   o         o       o    |     |
#             |                          |     |
#             |                          |     |
#    y <---   |   o         +       o    |     B
#             |                          |     |
#             |                          |     |
#             |   o        o        o    |     |
#             -----------------------------  --
#             |--------    H    --------|
#                          |
#                          V
#                          z

y1= HCol/2
z1= BCol/2
ColSecTag=1

section('Fiber', ColSecTag, '-GJ', 1e+10)

# Create the concrete core fibers
patch('rect', IDconcC, 10, 1 , c-y1, c-z1, y1-c, z1-c)

# Create the concrete cover fibers (top, bottom, left, right)
patch('rect', IDconcU, 10, 1,  -y1, z1-c,   y1,  z1)
patch('rect', IDconcU, 10, 1,  -y1,  -z1,   y1, c-z1)
patch('rect', IDconcU, 2,  1,  -y1, c-z1, c-y1, z1-c)
patch('rect', IDconcU, 2,  1, y1-c, c-z1,   y1, z1-c)


# Create the reinforcing fibers (left, middle, right)
layer('straight', IDreinf, numBarsSec, barAreaSec, y1-c, z1-c, y1-c, c-z1)
#layer('straight', IDreinf, 2, As,  0.0, z1-c,  0.0, c-z1)
layer('straight', IDreinf, numBarsSec, barAreaSec, c-y1, z1-c, c-y1, c-z1)

# Section for bond slip
SecTag2 = 2
section('Fiber',SecTag2,'-GJ',1e+10)


# Create the concrete core fibers
patch('rect',IDconcC,10,1 ,c-y1, c-z1, y1-c, z1-c)

# Create the concrete cover fibers (top, bottom, left, right)
patch('rect', IDconcU, 10, 1,  -y1, z1-c,   y1,  z1)
patch('rect', IDconcU, 10, 1,  -y1,  -z1,   y1, c-z1)
patch('rect', IDconcU, 2,  1,  -y1, c-z1, c-y1, z1-c)
patch('rect', IDconcU, 2,  1, y1-c, c-z1,   y1, z1-c)


# Create the reinforcing fibers (left, middle, right)
layer('straight', IDSP, numBarsSec, barAreaSec, y1-c, z1-c, y1-c, c-z1)
#layer('straight', IDSP, 2, As,  0.0, z1-c,  0.0, c-z1)
layer('straight', IDSP, numBarsSec,  barAreaSec, c-y1, z1-c, c-y1, c-z1)

# Creating Elements

ColTransfTag = 1
geomTransf('Linear', ColTransfTag)

ZL_eleTag=1
element('zeroLengthSection', ZL_eleTag, 1, 2, SecTag2, '-orient', 0., 1., 0., 1., 0., 0.)

ColeleTag = 2

# Defining Fiber Elements as ForceBeamColumn
#element('nonlinearBeamColumn', eleTag, 1, 2, numIntgrPts, ColSecTag, ColTransfTag)
ColIntTag=1
# beamIntegration('Lobatto',ColIntTag,ColSecTag,numIntgrPts)
beamIntegration('HingeRadau', ColIntTag, ColSecTag, Lpt, ColSecTag, 1e-10, ColSecTag)
element('forceBeamColumn', ColeleTag, 2, 3, ColTransfTag, ColIntTag, '-mass', 0.0)

# Setting Recorders


recorder('Node', '-file','DFree.out', '-time','-node', 3, '-dof', 1, 2, 3, 'disp')
recorder('Node', '-file','/DBase.out', '-time', '-node', 1, '-dof', 1, 2, 3, 'disp')
recorder('Node', '-file', 'RBase.out', '-time', '-node', 2, '-dof', 1, 2, 3, 'reaction')
recorder('Element', '-file', 'StressStrain.out', '-time','-ele', 2, 'section', '1', 'fiber', str(y1-c)+','+str(z1-c),'mat','3','stressStrain')  #Rbl,0, IDreinf
#recorder('Element', '-file', 'StressStrain2.out','-time','-ele', 2, 'section', '1', 'fiber', str(-Dprime)+', 0.0','mat','1','stressStrain')  #Rbl,0, IDreinf
#recorder('Element', '-file', 'StressStrain3.out','-time','-ele', 2, 'section', '1', 'fiber', str(-DCol)+', 0.0','mat','2','stressStrain')
# recorder('Element', '-file', datadir+'Data-2c/DCol.out','-time', '-ele', 1, 'deformations')

#------------------------------------------------------------------------------ 
#|                      DISPLACEMENT CONTROL RUN
#------------------------------------------------------------------------------


#defining gravity loads
timeSeries('Linear', 1)
pattern('Plain', 1, 1)
load(3, 0.0, -PCol, 0.0)

Tol = 1e-8 # convergence tolerance for test
NstepGravity = 10
DGravity = 1/NstepGravity
integrator('LoadControl', DGravity) # determine the next time step for an analysis
numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
system('BandGeneral','-piv') # how to store and solve the system of equations in the analysis
constraints('Plain') # how it handles boundary conditions
test('NormDispIncr', Tol, 10) # determine if convergence has been achieved at the end of an iteration step
algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration
analysis('Static') # define type of analysis static or transient
analyze(NstepGravity) # apply gravity

loadConst('-time', 0.0) #maintain constant gravity loads and reset time to zero
 
#applying Dynamic Ground motion analysis
timeSeries('Linear', 2)
pattern('Plain', 2, 2)
load(3, 1, 0, 0)
constraints('Plain') # how it handles boundary conditions
numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
system('BandGeneral') # how to store and solve the system of equations in the analysis
test('EnergyIncr', Tol, 10) # determine if convergence has been achieved at the end of an iteration step
algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration


dy = 13.65 * mm
dy1 = 12.4 * mm
ddmax=10
nIncres=100
FracStr=0.3
stDisp=0.
pdisp=[]
disp=[]

for i in range(0,4):
    pdisp.append((i+1)*0.25*dy1)
    pdisp.append(-(i+1)*0.25*dy1)
    
    
for i in range(2,5):
    disp.append(i*dy/2)
    disp.append(-i*dy/2)
    disp.append(i*dy/2)
    disp.append(-i*dy/2)
    disp.append(i*dy/2)
    disp.append(-i*dy/2)
    
    
for i in range (3,ddmax+1):
    disp.append(i*dy)
    disp.append(-i*dy)
    disp.append(i*dy)
    disp.append(-i*dy)
    disp.append(i*dy)
    disp.append(-i*dy)    
    
    
TolStatic=1e-6
maxNumIterStatic=6
testTypeStatic='EnergyIncr'
algorithTypeStatic='Newton'
fmt1='%s Pushover analysis: 3 %3i, dof %1i, Curv=%4f /%s'
count=0

Atest = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
Algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

for xdisp in pdisp:
    count=count+1
    netdisp=xdisp-stDisp
    dD=netdisp/nIncres
    integrator('DisplacementControl',3,1,dD)
    analysis('Static')
    ok=analyze(nIncres)       
    
    for j in Algorithm:
        

        if ok != 0:
            if j < 4:
                algorithm(Algorithm[j], '-initial')
                
            else:
                algorithm(Algorithm[j])
                
            test(testTypeStatic, TolStatic, maxNumIterStatic)
            ok = analyze(1)
            algorithm(algorithTypeStatic)                            
            if ok == 0:
                break
        else:
            continue
    stDisp=xdisp
    u3 = nodeDisp(3, 1)
    print(u3)
    
    
for xdisp in disp:
    count=count+1
    netdisp=xdisp-stDisp
    dD=netdisp/nIncres
    integrator('DisplacementControl',3,1,dD)
    analysis('Static')
    ok=analyze(nIncres)       
    
    for j in Algorithm:
        

        if ok != 0:
            if j < 4:
                algorithm(Algorithm[j], '-initial')
                
            else:
                algorithm(Algorithm[j])
                
            test(testTypeStatic, TolStatic, maxNumIterStatic)
            ok = analyze(5)
            algorithm(algorithTypeStatic)                            
            if ok == 0:
                break
        else:
            continue
    stDisp=xdisp
    u3 = nodeDisp(3, 1)
    print(u3)

#Force Displacement Plot
#    
#    
#d=open("DFree.out")
#F=open("RBase.out")
#
#
#linesd = d.readlines()
#linesf = F.readlines()
#x = [line.split()[1] for line in linesd]
#y = [line.split()[-1] for line in linesf]
#
#X=[float(i) for i in x]
#Y=[float(i) for i in y]
#
#plt.figure(1)    
#plt.plot(X,Y)
#plt.title('Example for ChiChi EQ w/c=0.4', fontsize=32)
#plt.xlabel('Diplacement (in)', fontsize=24)
#plt.ylabel('BaseShear (kip)', fontsize=24)
#plt.tick_params(direction='out',axis='both',labelsize=20)
#plt.grid()
#plt.show()

# Steel Stress Strain Analysis


SteelStressStrain=open("StressStrain.out")
linesSteelStressStrain=SteelStressStrain.readlines()
StlStress=[line.split()[1] for line in linesSteelStressStrain]
StlStrain=[line.split()[2] for line in linesSteelStressStrain]
sigmaStl=[float(i) for i in StlStress]
epsilonStl=[float(i) for i in StlStrain]

plt.figure(2)
plt.plot(epsilonStl,sigmaStl)
plt.title('Example of Steel Response for ChiChi EQ w/c=0.4', fontsize=32)
plt.xlabel('strain (in/in)', fontsize=24)
plt.ylabel('Stress (ksi)', fontsize=24)
plt.tick_params(direction='out',axis='both',labelsize=20)
plt.grid()
plt.show()