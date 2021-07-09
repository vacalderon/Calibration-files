# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:12:06 2019

@author: pchi893
"""
# ------------------------------------------------------------------------------
# |                      IMPORTS
# ------------------------------------------------------------------------------

from openseespy.opensees import *
import math

# ------------------------------------------------------------------------------
# |                      UNITS
# ------------------------------------------------------------------------------
inch=1.                 # define basic units
sec=1.0
kip=1.0
rad=1.0
LunitTXT="Inches"       # define basic-unit text for output
FunitTXT="Kips"                  
TunitTXT="Seconds"                 

ft=inch*12.0            # define dependent units
lbf=kip/1000.0
ksi=kip/(inch**2)
psi=lbf/(inch**2)
in2=inch*inch

m=inch*39.3701          # define metric basic units
N=kip*2.2481e-4

mm=m/1000.0             # define dependent metric units
cm=m/100.0          
MPa=N/(mm**2.0)
Pa=MPa/1000.0
GPa=MPa*1000.0
kN=N*1000.0
mm2=mm**2.0
m2=m*m
PI=2*math.asin(1.0)     # define constants
g=9.81*m/(sec**2.0)
deg=rad/180.0*PI

U=1e10                  # a really large number
u=1/U                   # a really small number mm



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


wipe()

# ------------------------------------------------------------------------------
#                           GENERATE GEOMETRY
# ------------------------------------------------------------------------------

model('basic', '-ndm', 2, '-ndf', 3)
LCol = 8.0 * ft  # column length
PCol=266.4573*kip
Weight = PCol  # superstructure weight

# define section geometry
DCol = 24.0 * inch  # Column Diameterepth


#PCol = Weight  # nodal dead-load weight per column
Mass = PCol / g

#  Node, Dx, Dy, Rz
node(1, 0.0, 0.0)
node(2, 0.0, 0.0)
node(3, 0.0, LCol)
#  Node Fixity
fix(1, 1, 1, 1)
fix(2, 1, 0, 0)

#Setting mass for analysis
mass(3, Mass, 1e-9, 0.0)

# ------------------------------------------------------------------------------
#                           MATERIAL PARAMETERS
# ------------------------------------------------------------------------------

IDconcC = 1  # material ID tag -- confined cover concrete
IDconcU = 2  # material ID tag -- unconfined cover concrete
IDreinf = 3  # material ID tag -- reinforcement
IDSP    = 4  # material ID tag -- Strain Penetration
# Define materials for nonlinear columns
# ------------------------------------------
# Longitudinal steel properties
Fy = 574.0 * MPa  # STEEL yield stress
Fu = 753.3 * MPa  # Steel Ultimate Stress
Es = 160000.0 * MPa  # modulus of steel
Bs = 0.012  # strain-hardening ratio
R0 = 20.0  # control the transition from elastic to plastic branches
cR1 = 0.9  # control the transition from elastic to plastic branches
cR2 = 0.08  # control the transition from elastic to plastic branches
a1=0.039
a2=1.0
a3=0.029
a4=1.0
c = 1.25 * inch  # Column cover to reinforcing steel NA.
numBarsSec = 16  # number of uniformly-distributed longitudinal-reinforcement bars
barAreaSec = 0.60 * in2  # area of longitudinal-reinforcement bars
dbl = (7/8) * inch


# Transverse Steel Properties
fyt = 574.0 * ksi  # Yield Stress of Transverse Steel
Ast = 0.11 * in2  # Area of transverse steel
dbt = 0.375 * inch  # Diameter of transverse steel
st = 2.0 * inch  # Spacing of spiral
Dprime = DCol - 2 * c - dbt*0.5  # Inner core diameter
Rbl = Dprime * 0.5 - dbt*0.5 - dbl * 0.5  # Location of longitudinal bar

# nominal concrete compressive strength
fpc = 39.8 * MPa  # CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
Ec = 57.0 * ksi * math.sqrt(fpc / psi)  # Concrete Elastic Modulus

# unconfined concrete
fc1U = -fpc;  # UNCONFINED concrete (todeschini parabolic model), maximum stress
eps1U = -0.002  # strain at maximum strength of unconfined concrete
fc2U = 0.0 * MPa  # ultimate stress
eps2U = -0.0064  # strain at ultimate stress
lambdac = 0.1  # ratio between unloading slope at $eps2 and initial slope $Ec


# Confined Concrete Properties

fc = -56.1 * MPa
eps1 = -0.0061
fc2 = -42.2 * MPa
eps2 = -0.0247

# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
uniaxialMaterial('Concrete01', IDconcC, fc, eps1, fc2, eps2)

# Cover concrete (unconfined)
uniaxialMaterial('Concrete01', IDconcU, fc1U, eps1U, fc2U, eps2U)

# STEEL
# Reinforcing steel 
params = [R0, cR1, cR2]
#                        tag  fy E0    b
uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs,R0, cR1, cR2, a1, a2, a3, a4)


# STRAIN PENETRATION MATERIAL
SPalpha=0.4
SPsy=0.1*((dbl*Fy)*(2*SPalpha+1)/(4000*((-fc)**0.5)))**(1/SPalpha)+0.013
SPsu=35*SPsy
SPb=0.45
SPR=1.01


#uniaxialMaterial StrPen01   Tag  fy  sy fu su b  R
uniaxialMaterial('Bond_SP01',IDSP, Fy,SPsy,Fu,SPsu,SPb,SPR)

#-------------------------------------------------------------------------
#               DEFINE PLASTICE HIGE PROPERTIES
#-------------------------------------------------------------------------

k=0.2*(Fu/Fy - 1)
if k > 0.08:
    k=0.08
Leff=LCol
Lpc=k*Leff + 0.4*DCol
gamma=0.33 #Assuming unidirectional action
Lpt=Lpc+gamma*DCol
# FIBER SECTION properties -------------------------------------------------------------
# Define cross-section for nonlinear columns
# ------------------------------------------

# set paramaters for Section 1
ColSecTag = 1
ri = 0.0
ro = DCol / 2.0
nfCoreR = 8
nfCoreT = 8
nfCoverR = 2
nfCoverT = 8
rc = ro - c
theta = 360.0 / numBarsSec

section('Fiber', ColSecTag,'-GJ',1e+10)

# Create the concrete fibers
patch('circ', 1, nfCoreT, nfCoreR, 0.0, 0.0, ri, rc, 0.0, 360.0)  # Define the core patch
patch('circ', 2, nfCoverT, nfCoverR, 0.0, 0.0, rc, ro, 0.0, 360.0)  # Define Cover Patch

# Create the reinforcing fibers
layer('circ', 3, numBarsSec, barAreaSec, 0.0, 0.0, Rbl, theta, 360.0)

# Set parameters for ZeroLength Element

SecTag2 = 2
section('Fiber',SecTag2,'-GJ',1e+10)

# Create the concrete fibers
patch('circ', 1, nfCoreT, nfCoreR, 0.0, 0.0, ri, rc, 0.0, 360.0)  # Define the core patch
patch('circ', 2, nfCoverT, nfCoverR, 0.0, 0.0, rc, ro, 0.0, 360.0)  # Define Cover Patch

# Create the reinforcing fibers for Strain Penetration
layer('circ', IDSP, numBarsSec, barAreaSec, 0.0, 0.0, Rbl, theta, 360.0)

# Creating Elements

ColTransfTag = 1
geomTransf('Linear', ColTransfTag)

ZL_eleTag=1
element('zeroLengthSection', ZL_eleTag, 1, 2, SecTag2, '-orient', 0., 1., 0., 1., 0., 0.)

ColeleTag = 2

# Defining Fiber Elements as ForceBeamColumn

ColIntTag=1
beamIntegration('HingeRadau', ColIntTag, ColSecTag, Lpt, ColSecTag, 1e-10, ColSecTag)
element('forceBeamColumn', ColeleTag, 2, 3, ColTransfTag, ColIntTag, '-mass', 0.0)

# ------------------------------------------------------------------------------
#                           SETTING RECORDERS
# ------------------------------------------------------------------------------

recorder('Node', '-file','DFree.out', '-time','-node', 3, '-dof', 1, 2, 3, 'disp')
recorder('Node', '-file','/DBase.out', '-time', '-node', 1, '-dof', 1, 2, 3, 'disp')
recorder('Node', '-file', 'RBase.out', '-time', '-node', 2, '-dof', 1, 2, 3, 'reaction')
recorder('Element', '-file', 'StressStrain.out', '-time','-ele', 2, 'section', '1', 'fiber', str(Rbl), '0.0','3','stressStrain')

#------------------------------------------------------------------------------
#|                      NLTH Analysis Run
#------------------------------------------------------------------------------

#
# 1. Defining gravity loads
#

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

#
# 2. Applying Dynamic Ground motion analysis
#

GMdirection = 1
GMfile = '95.0_96.0_01.g3'
GMfact = 981.0



Lambda = eigen('-fullGenLapack', 1) # eigenvalue mode 1
Omega = math.pow(Lambda[0], 0.5)
betaKcomm = 2 * (0.02/Omega)
xDamp = 0.04		# 4 % damping ratio
alphaM = 0.0		# M-prop. damping; D = alphaM*M	
betaKcurr = 0.0		# K-proportional damping;      +beatKcurr*KCurrent
betaKinit = 0.0     # initial-stiffness proportional damping      +beatKinit*Kini

rayleigh(alphaM,betaKcurr, betaKinit, betaKcomm) # RAYLEIGH damping

# Uniform EXCITATION: acceleration input
IDloadTag = 400		# load tag
dt = 0.005			# time step for input ground motion
GMfatt = 1.0		# data in input file is in g Unifts -- ACCELERATION TH
maxNumIter = 10
timeSeries('Path', 2, '-dt', dt, '-filePath', GMfile, '-factor', GMfact)
pattern('UniformExcitation', IDloadTag, GMdirection, '-accel', 2) 

wipeAnalysis()

constraints('Transformation')
numberer('Plain')
system('BandGeneral')
test('EnergyIncr', Tol, maxNumIter)
algorithm('ModifiedNewton')
NewmarkGamma = 0.5
NewmarkBeta = 0.25
integrator('Newmark', NewmarkGamma, NewmarkBeta)
analysis('Transient')
DtAnalysis = 0.0001
TmaxAnalysis = 97.58
Nsteps =  int(TmaxAnalysis/ DtAnalysis)
analyze(Nsteps, DtAnalysis)
tCurrent = getTime()
u3 = nodeDisp(3, 1)
print("u = ", u3)

wipe()