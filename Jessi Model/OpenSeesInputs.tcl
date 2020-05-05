# OpenSees inputs for Material Models

# Confined Concrete 
set fc1C 	[expr -56.1*$MPa];
set eps1C	-0.0061; 
set fc2C 	[expr -42.2*$MPa];
set eps2C 	-0.0247; 

 # Unconfined Concrete 
set fc1U 	[expr -39.8*$MPa];
set eps1U	-0.0020; 
set fc2U 	[expr -0.0*$MPa];
set eps2U	-0.0064; 

 # ReinforcingSteel Properties 
set Fy 		[expr 574.0*$MPa];
set Fu		[expr 753.3*$MPa];
set Es		[expr 160000.0*$MPa];
set Esh      [expr 8609.9*$MPa];
set esh		0.01; 
set esult	0.09500; 
set Isr		2.67; 
set beta		1.00; 
set r		0.55; 
set gama		0.50; 
set a1		4.30; 
set limit	0.90; 
set Cf    	0.500; 
set alpha	0.450; 
set Cd   	0.600; 

 # Section Geometry 
set DSec     [expr 23.0*$in]; 
set coverSec [expr 1.2500*$in]; 
set numBarsSec  16; 
set dbar [expr 0.750*$in]; 
set barAreaSec [expr 0.442*$in2]; 

 # Column Input Parameters 
set LCol [expr  8.0*$ft]; 
set PCol [expr 580.6*$kN]; 

 # Other Parameters 
set dy 			[expr 1.050000*$in]; 
set dy1			[expr 0.780000*$in]; 
set ddmax		 6; 
set nIncres		100; 
set FracStr		0.300000; 

