# Coded by Wangcheng Zhang on 29/04/21 for Eulerian depth-averaged MPM progressive landslide failure - staggerd meshes
# ========import lines==========
import numpy as np
import math
import sys
from Geometry import *
from Initiation import *

Ang_N[:,:]=np.arctan(math.tan(AngC)*np.exp(-(abs(Y_N[:,:])/H_slope)*math.tan(AngC)))

COORD_N=np.reshape((np.array((X_N,Y_N),dtype='float32').T),(Ny*Nx,2))
suS_N[:,:]=TpS
suW_N[:,:]=np.minimum(((abs(X_N[:,:])/Sof_X)**ns+(abs(Y_N[:,:])/Sof_Y)**ns)*T0/2.0,TpW*2.0)

HC_N[:,:]=HIni*np.cos(Ang_N[:,:])
HN_N[:,:]=HIni*np.cos(Ang_N[:,:])
HC_HX[:,:]=HIni
HC_HX[:,1:-1]=(HC_N[:,:-1]+HC_N[:,1:])/2
HC_HY[1:-1]=(HC_N[:-1]+HC_N[1:])/2
HN_HX[:]=HC_HX[:]
HN_HY[:]=HC_HY[:]
Sigma_N[:,:,2]=DensNet*g*np.cos(Ang_N[:,:])*HC_N[:,:]/2
Sigma_N[:,:,0]=Sigma_N[:,:,2]*K0
Sigma_N[:,:,1]=Sigma_N[:,:,2]*K0

p_N[:,:,0]=np.sum(Sigma_N[:,:,:3]/3,axis=-1)
SIJ_N[:,:,:3]=Sigma_N[:,:,:3]-p_N[:,:]

PX_N[:,:]=Sigma_N[:,:,0]*HC_N[:,:]*dy
PY_N[:,:]=Sigma_N[:,:,1]*HC_N[:,:]*dx
TX_HXY[1:-1,1:-1]=(Sigma_N[:-1,:-1,3]*HC_N[:-1,:-1]+Sigma_N[:-1,1:,3]*HC_N[:-1,1:]+Sigma_N[1:,:-1,3]*HC_N[1:,:-1]+Sigma_N[1:,1:,3]*HC_N[1:,1:])*dx/4
TY_HXY[1:-1,1:-1]=(Sigma_N[:-1,:-1,3]*HC_N[:-1,:-1]+Sigma_N[:-1,1:,3]*HC_N[:-1,1:]+Sigma_N[1:,:-1,3]*HC_N[1:,:-1]+Sigma_N[1:,1:,3]*HC_N[1:,1:])*dy/4
taugX_HX[:,1:-1]=DensNet*gx*HC_HX[:,1:-1]
taugY_HY[1:-1,:]=DensNet*(-g)*np.sin(Ang_N[1:,:])*HC_HY[1:-1,:]
tauWX_HX[:]=-taugX_HX[:]
tauWY_HY[:]=-taugY_HY[:]
tauWX_HX[:,1:-1]=-((PX_N[:,:-1]-PX_N[:,1:])+(TX_HXY[1:,1:-1]-TX_HXY[:-1,1:-1]))/dx/dy-taugX_HX[:,1:-1]
tauWY_HY[1:-1,:]=-((PY_N[:-1,:]-PY_N[1:,:])+(TY_HXY[1:-1,1:]-TY_HXY[1:-1,:-1]))/dx/dy-taugY_HY[1:-1,:]
tauW_N[:,:,0]=(tauWX_HX[:,:-1]+tauWX_HX[:,1:])/2
tauW_N[:,:,1]=(tauWY_HY[:-1,:]+tauWY_HY[1:,:])/2