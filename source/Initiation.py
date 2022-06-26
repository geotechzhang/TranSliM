# Coded by Wangcheng Zhang on 29/04/21 for Eulerian depth-averaged MPM progressive landslide failure - staggerd meshes
# ========import lines==========
import numpy as np
import math
import sys
from Geometry import *
#------------------nodal point initiation------------
x_N=np.linspace(-Length/2, Length/2, Nx,dtype='float32')
y_N=np.linspace(-Width*2.0/3.0, Width/3.0, Ny,dtype='float32')
X_N,Y_N=np.meshgrid(x_N, y_N) 
X_N=X_N.astype(np.float32)
Y_N=Y_N.astype(np.float32)

VelUN_HX=np.zeros((Ny,Nx+1),dtype='float32') #vel of each node, HALF+1 to have gohst nodes boundary conditions
VelUC_HX=np.zeros((Ny,Nx+1),dtype='float32')
VelVN_HY=np.zeros((Ny+1,Nx),dtype='float32') 
VelVC_HY=np.zeros((Ny+1,Nx),dtype='float32')
VelUC_N=np.zeros((Ny,Nx),dtype='float32') #vel of each node
VelVC_N=np.zeros((Ny,Nx),dtype='float32')
QX_HX=np.zeros((Ny,Nx+1),dtype='float32') #FLUX
QX_N=np.zeros((Ny,Nx),dtype='float32')
QXU_N=np.zeros((Ny,Nx),dtype='float32')
QX_HXY=np.zeros((Ny+1,Nx+1),dtype='float32')
QXV_HXY=np.zeros((Ny+1,Nx+1),dtype='float32')
QY_HY=np.zeros((Ny+1,Nx),dtype='float32')
QY_N=np.zeros((Ny,Nx),dtype='float32')
QYV_N=np.zeros((Ny,Nx),dtype='float32')
QY_HXY=np.zeros((Ny+1,Nx+1),dtype='float32')
QYU_HXY=np.zeros((Ny+1,Nx+1),dtype='float32')
HN_N=np.zeros((Ny,Nx),dtype='float32')  #height of each node
HC_N=np.zeros((Ny,Nx),dtype='float32')
HC_HX=np.zeros((Ny,Nx+1),dtype='float32')
HN_HX=np.zeros((Ny,Nx+1),dtype='float32')
HC_HY=np.zeros((Ny+1,Nx),dtype='float32')
HN_HY=np.zeros((Ny+1,Nx),dtype='float32')
Ang_N=np.zeros((Ny,Nx),dtype='float32')
DUP_N=np.zeros((Ny,Nx,2),dtype='float32')
tauW_N=np.zeros((Ny,Nx,2),dtype='float32')
tauW_N_Tri=np.zeros((Ny,Nx,2),dtype='float32')
tauWX_HX=np.zeros((Ny,Nx+1),dtype='float32')
tauWY_HY=np.zeros((Ny+1,Nx),dtype='float32')
suW_N=np.zeros((Ny,Nx),dtype='float32')
suS_N=np.zeros((Ny,Nx),dtype='float32')
DEpsi_N=np.zeros((Ny,Nx,4),dtype='float32')
SIJ_N=np.zeros((Ny,Nx,4),dtype='float32')
Deviaq_N_Tri=np.zeros((Ny,Nx),dtype='float32')
PlaCorS_N=np.zeros((Ny,Nx,1),dtype='float32')
DPEEQ_N=np.zeros((Ny,Nx),dtype='float32')
PEEQ_N=np.zeros((Ny,Nx),dtype='float32')
p_N=np.zeros((Ny,Nx,1),dtype='float32') 
COORD_N=np.zeros((Ny*Nx,2),dtype='float32') 
Sigma_N=np.zeros((Ny,Nx,4),dtype='float32')
PX_N=np.zeros((Ny,Nx),dtype='float32')
PY_N=np.zeros((Ny,Nx),dtype='float32')
TX_HXY=np.zeros((Ny+1,Nx+1),dtype='float32')
TY_HXY=np.zeros((Ny+1,Nx+1),dtype='float32')
taugX_HX=np.zeros((Ny,Nx+1),dtype='float32')
taugY_HY=np.zeros((Ny+1,Nx),dtype='float32')
tauWX_HX=np.zeros((Ny,Nx+1),dtype='float32')
tauWY_HY=np.zeros((Ny+1,Nx),dtype='float32')
PlaCorW_N=np.zeros((Ny,Nx,1),dtype='float32')
DampPX_N=np.zeros((Ny,Nx),dtype='float32')
DampPY_N=np.zeros((Ny,Nx),dtype='float32')
temp=np.zeros((Ny,Nx),dtype='float32')
