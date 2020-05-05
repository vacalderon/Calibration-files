# Fiber-based model for Reinforced Concrete Column subjected to three cycle-set loading

# SET UP MODEL ----------------------------------------------------------------------------
wipe;                               # clear memory of all past model definitions
model basic -ndm 2 -ndf 3;  # Define the model builder, ndm=#dimension, ndf=#dofs
source Units/LibUnitsUS.tcl;                # define units
source Data/OpenSeesInputs.tcl;         # inputs from Cumbia calculation

# MATERIAL parameters -------------------------------------------------------------------
set IDconcCore  401;        # material ID tag -- confined core concrete
set IDconcCover 402;        # material ID tag -- unconfined cover concrete
set IDreinf     403;        # material ID tag -- steel with maximum and minimum strains
set IDstrap     404;        # material ID tag -- strain penetration element

# set material                  tag     f'c     ec0    f'cu  ecu
uniaxialMaterial Concrete01 $IDconcCore $fc1C $eps1C $fc2C $eps2C;          # build core concrete (confined)
uniaxialMaterial Concrete01 $IDconcCover $fc1U $eps1U $fc2U $eps2U;         # build cover concrete (unconfined)
uniaxialMaterial ReinforcingSteel $IDreinf  $Fy  $Fu    $Es $Esh $esh $esult;# -CMFatigue  $Cf  $alpha  $Cd ;#-GABuck   $Isr    $beta   $r      $gama ;
                        #   tag     yield   rebarSlip   ultimate    rebarSlipU  hardeningRatio  pinchingFactor
uniaxialMaterial Bond_SP01 $IDstrap $Fy     0.021320995       $Fu         0.64         0.1             0.6;

# Nodal coordinates
#    tag  X Y
node 101  0 0;          # node 1, x-coor y-coor
node 102  0 0;
node 103  0 $LCol;

# point constraints
#   node DX DY DZ
fix 101  1  1  1;   
fix 102  1  0  0;

equalDOF 101 102 2

# section GEOMETRY -------------------------------------------------------------
set ColSecTag 200;      # tag for column section
set spSecTag  201;      # tag for strain-penetration section 
set ri 0.0;             # inner radius of the section, only for hollow sections
set ro [expr $DSec/2];  # overall (outer) radius of the section
set nfCoreR 8;          # number of radial divisions in the core (number of "rings")
set nfCoreT 8;          # number of theta divisions in the core (number of "wedges")
set nfCoverR 4;         # number of radial divisions in the cover
set nfCoverT 8;         # number of theta divisions in the cover

# Define the fiber section for column
section fiberSec $ColSecTag  {
    set rc [expr $ro-$coverSec];                    # Core radius
    patch circ $IDconcCore $nfCoreT $nfCoreR 0 0 $ri $rc 0 360;     # Define the core patch
    patch circ $IDconcCover $nfCoverT $nfCoverR 0 0 $rc $ro 0 360;  # Define the cover patch
    set theta [expr 360.0/$numBarsSec];     # Determine angle increment between bars
    layer circ $IDreinf $numBarsSec $barAreaSec 0 0 $rc $theta 360; # Define the reinforcing layer
}
# Define the fiber section for the strain penetration section
section fiberSec $spSecTag  {
    set rc [expr $ro-$coverSec];                    # Core radius
    patch circ $IDconcCore $nfCoreT $nfCoreR 0 0 $ri $rc 0 360;     # Define the core patch
    patch circ $IDconcCover $nfCoverT $nfCoverR 0 0 $rc $ro 0 360;  # Define the cover patch
    set theta [expr 360.0/$numBarsSec];     # Determine angle increment between bars
    layer circ $IDstrap $numBarsSec $barAreaSec 0 0 $rc $theta 360; # Define the reinforcing layer
}

# Calculated parameters for column
set Mass [expr $PCol/$g];       # nodal mass
set ACol [expr 0.25*$PI*$DSec*$DSec];   # cross-sectional area of column
set IzCol [expr 1./64.*$PI*$DSec*$DSec*$DSec*$DSec];        # section moment of inertia

# parameters particular to this model
set IDctrlNode 103;         # node where displacement is read for displacement control
set IDctrlDOF 1;            # degree of freedom of displacement read for displacement control
set iSupportNode "101";     # define support node, if needed.

# nodal masses:
mass 103 $Mass  1e-9 0.;        # node#, Mx My Mz, Mass=Weight/g, neglect rotational inertia at nodes

# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
set colTag 300;     # column tag number
set ColTransfTag 1;             # associate a tag to column transformation
set ColTransfType Linear ;          # options, Linear PDelta Corotational 
geomTransf $ColTransfType $ColTransfTag ;   

# define the zero-length strain penetration element
set spTag 301;

# element connectivity:
set numIntgrPts 5;                              # number of integration points for force-based element
#                           tag     iNode   jNode   numberIntPoints sectionTag  transformTag
element forceBeamColumn     $colTag 102     103     $numIntgrPts    $ColSecTag  $ColTransfTag;  

#                           tag     nodeI  nodeJ    materialTag 
element zeroLengthSection   $spTag  101     102     $spSecTag   ;   

# Recorders
recorder Node -file Data/DFree1.out -time -node 103 -dof 1 disp; 
recorder Element -file Data/NorthBarSSIP1.out -time -ele $colTag section 1 fiber $rc 0 $IDreinf stressStrain;

# define GRAVITY -------------------------------------------------------------
pattern Plain 1 Linear {
   load 103 0 -$PCol 0
}

# Gravity-analysis parameters -- load-controlled static analysis
set Tol 1.0e-8;                 # convergence tolerance for test
constraints Plain;              # how it handles boundary conditions
numberer Plain;                 # renumber dof's to minimize band-width (optimization), if you want to
system BandGeneral;             # how to store and solve the system of equations in the analysis
test NormDispIncr $Tol 6 ;      # determine if convergence has been achieved at the end of an iteration step
algorithm Newton;               # use Newton's solution algorithm: updates tangent stiffness at every iteration
set NstepGravity 1;             # apply gravity in 1 step
set DGravity [expr 1./$NstepGravity];   # first load increment;
integrator LoadControl $DGravity;   # determine the next time step for an analysis
analysis Static;                # define type of analysis static or transient
analyze $NstepGravity;          # apply gravity
# ------------------------------------------------- maintain constant gravity loads and reset time to zero
loadConst -time 0.0

puts "Model Built"
puts " "

# testing the model
# Apply three-set cyclic loading to the element


# create load pattern for lateral pushover load
set Hload 1;            # define the lateral load as a 1 so that the pseudo time equals the lateral load 
pattern Plain 200 Linear {;         # define load pattern -- generalized
    load 103 $Hload 0.0 0.0
}

# STATIC-ANALYSIS parameters
constraints Plain
numberer Plain
system BandGeneral
test EnergyIncr 1.e-8 6 0
algorithm Newton      

# Three-Set-Cyclic Loading Parameters
# set NumCycles 3;          # Number of cycles at each displacement 
# set stDisp        0;          # initial displacement 

# Build Parameters for Loading Protocals
# set prey " ";
# for {set i 1} {$i <= 4} {incr i 1} {
    # set prey "$prey [expr $i*0.25*$dy1]"
# }
# set iDisp " ";
# set iDisp "$iDisp [expr 1*$dy]"
# set iDisp "$iDisp [expr 1.5*$dy]"
# for {set i 2} {$i <= $ddmax} {incr i 1} {
    # set iDisp "$iDisp [expr $i*$dy]"
# }


# Three-Set-Cyclic Loading Parameters 
set stDisp      0;          # initial displacement 

# Build displacement history for the model
set pdisp " ";
set disp " ";
for {set i 1} {$i <= 4} {incr i 1} {
    set pdisp "$pdisp [expr $i*0.25*$dy1]"
    set pdisp "$pdisp [expr -$i*0.25*$dy1]"
}
for {set i 2} {$i <= 4} {incr i 1} {
    set disp "$disp [expr $i*$dy/2]"
    set disp "$disp [expr -$i*$dy/2]"
    set disp "$disp [expr $i*$dy/2]"
    set disp "$disp [expr -$i*$dy/2]"
    set disp "$disp [expr $i*$dy/2]"
    set disp "$disp [expr -$i*$dy/2]"
}
for {set i 3} {$i <= $ddmax} {incr i 1} {
    set disp "$disp [expr $i*$dy]"
    set disp "$disp [expr -$i*$dy]"
    set disp "$disp [expr $i*$dy]"
    set disp "$disp [expr -$i*$dy]"
    set disp "$disp [expr $i*$dy]"
    set disp "$disp [expr -$i*$dy]"
}
set disp "$disp 0"

set TolStatic 1.e-6;
set maxNumIterStatic 6;
set testTypeStatic EnergyIncr;
set algorithmTypeStatic Newton;
set fmt1 "%s Pushover analysis: CtrlNode %.3i, dof %.1i, Curv=%.4f /%s";

# Apply pre-yield displacement history to model
set count 0;
foreach xdisp  $pdisp {
    set count [expr $count+1];
    puts [format "Completed cycle $count for %.2f  F'y" [expr $xdisp/$dy1]];
    set netdisp [expr $xdisp-$stDisp];
    set dD [expr $netdisp/$nIncres];
                                  # $node       $dof    $incr   <$numIter   $Î”Umin         $Î”Umax>    
    integrator DisplacementControl  103         1       $dD     ;#1         [expr $dD/4]    $dD;
    analysis Static;
    set ok [analyze $nIncres];
    
    if {$ok != 0} {  
            puts "Trying Newton with Initial Tangent .."
            test NormDispIncr   $TolStatic      10 0
            algorithm Newton -initial
            set ok [analyze 1]
            test $testTypeStatic $TolStatic      $maxNumIterStatic    0
            algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
            puts "Trying Broyden .."
            algorithm Broyden 8
            set ok [analyze 1 ]
            algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
            puts "Trying NewtonWithLineSearch .."
            algorithm NewtonLineSearch 0.8 
            set ok [analyze 1]
            algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {;                # stop if still fails to converge
            puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
            return -1
        }; # end if
    set stDisp $xdisp;
}

# Apply plastic displacement history to model
foreach xdisp  $disp {
    set count [expr $count+1];
    puts [format "Completed cycle $count for %.1f dy, %.2f inches" [expr $xdisp/$dy] $xdisp];
    set netdisp [expr $xdisp-$stDisp];
    set dD [expr $netdisp/$nIncres];
                                #   $node       $dof    $incr       <$numIter   $Î”Umin         $Î”Umax>    
    integrator DisplacementControl  103         1       $dD         ;#2         [expr $dD/4]    $dD;
    analysis Static;
    set ok [analyze $nIncres];
    
    if {$ok != 0} {  
            puts "Trying Newton with Initial Tangent .."
            test NormDispIncr   $TolStatic      10 0
            algorithm Newton -initial
            set ok [analyze 4]
            test $testTypeStatic $TolStatic      $maxNumIterStatic    0
            algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
            puts "Trying Broyden .."
            algorithm Broyden 8
            set ok [analyze 5 ]
            algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {
            puts "Trying NewtonWithLineSearch .."
            algorithm NewtonLineSearch 0.8 
            set ok [analyze 4]
            algorithm $algorithmTypeStatic
        }
        if {$ok != 0} {;                # stop if still fails to converge
            puts [format $fmt1 "PROBLEM" $IDctrlNode $IDctrlDOF [nodeDisp $IDctrlNode $IDctrlDOF] $LunitTXT]
            return -1
        }; # end if
    set stDisp $xdisp;
}

wipe;

