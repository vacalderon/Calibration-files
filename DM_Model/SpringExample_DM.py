# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 14:32:32 2021

@author: VACALDER
"""
from openseespy.opensees import *
import matplotlib.pyplot as plt
# 1. Use variables to define units

from LibUnitsMUS import *

# ------------------------------
# Start of model generation
# -----------------------------

# remove existing model
wipe()

# set modelbuilder
model('basic', '-ndm', 1, '-ndf', 1)

# Define the nodes
node(1,0)
node(2,0)

#Restraints
fix(1, 1)

# 2. Define your model using nonlinear elements with elastic sections
#    in this case it will all be a nonlinear spring
# Define spring properties

# NL Material properties
Fy = 60 * ksi
E = 29000 * ksi
b = 0.1
ID_k1 = 101
lsr=6
CL=0.2

# uniaxialMaterial('Steel01', ID_k1, Fy, E, b)
uniaxialMaterial('ReinforcingSteel', ID_k1, Fy*(1-0.014*100*CL), 1.35*Fy, E, 0.1*E, 0.0074, 0.09, 
                 '-CMFatigue', 0.489, 0.54*(1+0.004*100*CL), 0.6) 
 

# translational spring
k1_elemTag = 1
element('zeroLength', k1_elemTag, 1, 2, '-mat', ID_k1, '-dir', 1)

# 3. Apply gravity loads (Only a translational load in this case)

# P= 100*lbf
# timeSeries('Linear', 1)
# pattern('Plain', 1, 1)
# Create the nodal load - command: load nodeID xForce yForce
# load(2, P)

#Set Recorders 

recorder('Node', '-file','DFree.out', '-time','-node', 2, '-dof', 1, 'disp')
recorder('Node', '-file', 'RBase.out', '-time','-node', 1, '-dof', 1, 'reaction')

# ------------------------------
# 7. Run analysis
# ------------------------------

# create SOE

disp = [0.002,-0.002, 0.002, -0.002, 0.002, -0.002, 
        0.01,-0.01, 0.01,-0.01, 0.01,-0.01, 
        0.02,-0.02, 0.02,-0.02, 0.02,-0.02,
        0.035,-0.035, 0.035,-0.035, 0.035,-0.035,
        0.06,-0.06, 0.06,-0.06, 0.06,-0.06,
        0.15]

#applying Dynamic Ground motion analysis
Tol=1e-8
stDisp=0.
nIncres=100
timeSeries('Linear', 1)
pattern('Plain', 1, 1)
load(2, 1)
constraints('Plain') # how it handles boundary conditions
numberer('Plain') # renumber dof's to minimize band-width (optimization), if you want to
system('BandGeneral') # how to store and solve the system of equations in the analysis
test('EnergyIncr', Tol, 10) # determine if convergence has been achieved at the end of an iteration step
algorithm('Newton') # use Newton's solution algorithm: updates tangent stiffness at every iteration

TolStatic=1e-6
maxNumIterStatic=10
testTypeStatic='EnergyIncr'
algorithTypeStatic='Newton'
fmt1='%s Pushover analysis: 3 %3i, dof %1i, Curv=%4f /%s'
count=0

Atest = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
Algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}



for xdisp in disp:
    count=count+1
    netdisp=xdisp-stDisp
    dD=netdisp/nIncres
    integrator('DisplacementControl',2,1,dD)
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
    u3 = nodeDisp(2, 1)
    print(u3)



wipe()
#Force Displacement Plot


d_file=open("DFree.out")
F=open("RBase.out")


linesd = d_file.readlines()
linesf = F.readlines()
x = [line.split()[1] for line in linesd]
y = [line.split()[0] for line in linesf]

X=[float(i) for i in x]
Y=[float(i) for i in y]

number_of_rows=min(len(X),len(Y))

plt.figure(1)
plt.plot(X,Y)
plt.title('Example for ChiChi EQ w/c=0.4', fontsize=32)
plt.xlabel('Diplacement (in)', fontsize=24)
plt.ylabel('BaseShear (kip)', fontsize=24)
plt.tick_params(direction='out',axis='both',labelsize=20)
plt.grid()
plt.show()