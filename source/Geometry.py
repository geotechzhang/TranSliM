# Coded by Wangcheng Zhang on 29/04/21 for Eulerian depth-averaged MPM progressive landslide failure - staggerd meshes
# ========import lines==========
import numpy as np
import math

#--------------------geometry-----------------------
Length=16000.0   
Width=20000.0
HIni=8.0
Slope_YL=2000.0
suS_XL=1000.0
AngC=12.0*math.pi/180.0
Sof_X=400.0
Sof_Y=400.0
ns=2
H_slope=1000.0

dx=20.0
dy=20.0
Nx=int(Length/dx)+1   #Node No limit (from 0)
Ny=int(Width/dy)+1   #Node No limit (from 0)
dt=0.1   #time step
dtdx=dt/dx
dtdy=dt/dy
T=300.0 #total time
Nt=int(T/dt)+1
ParNoCell=1  #Material points
dx_m=dx/ParNoCell
dy_m=dy/ParNoCell
Nx_m=int(Length/dx_m)+1
Ny_m=int(Width/dy_m)+1
output_No=int(1.0/dt)

		
		









