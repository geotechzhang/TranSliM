# Coded by Wangcheng Zhang on 29/04/21 for Eulerian depth-averaged MPM progressive landslide failure - staggerd meshes
# ========import lines==========
import numpy as np
import math
import sys
import datetime
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy import interpolate

from Geometry import *
from Initiation import *
from FirstStep import *
# ========parameter input==========
#----------------------material properties----------------
DensTol=1.8
DensNet=DensTol-1.0
GW=1656.0    #stiffness
TpW=18.0
TrW=6.0
T0=2.0
StW=TpW/TrW
UpLimW=0.2
SofSpeedW=-(TpW-TrW)/UpLimW  #hardening
SofTime=-0.0
GS=500.0
TpS=18.0
TrS=3.0
StS=TpS/TrS
PEEQLimS=0.2/min(dx,dy)
SofSpeedS=-(TpS-TrS)/PEEQLimS
#-----------------load properties---------
g=9.81
gx=0.0
K0=0.5

zeta=0.2 #damping
DampR=zeta*100.0
# ========Main program==========initiation: P, tauW,taug
for j in range (Nt):
	#--------------MAIN LOOP------------------
	QX_HX[:,1:-1]=np.where(VelUC_HX[:,1:-1]>0.0,VelUC_HX[:,1:-1]*HC_N[:,:-1],VelUC_HX[:,1:-1]*HC_N[:,1:]) 
	QY_HY[1:-1]=np.where(VelVC_HY[1:-1]>0.0,VelVC_HY[1:-1]*HC_N[:-1],VelVC_HY[1:-1]*HC_N[1:]) 
	HN_N[:,:]=(QX_HX[:,:-1]-QX_HX[:,1:])*dtdx+(QY_HY[:-1,:]-QY_HY[1:,:])*dtdy+HC_N[:,:]
	HN_HX[:,1:-1]=(HN_N[:,:-1]+HN_N[:,1:])/2
	HN_HY[1:-1]=(HN_N[:-1]+HN_N[1:])/2
	QX_N[:,:]=(QX_HX[:,:-1]+QX_HX[:,1:])/2
	QXU_N[:,:]=np.where(QX_N[:,:]>0.0,QX_N[:,:]*VelUC_HX[:,:-1],QX_N[:,:]*VelUC_HX[:,1:])
	QY_HXY[:,1:-1]=(QY_HY[:,:-1]+QY_HY[:,1:])/2
	QYU_HXY[1:-1]=np.where(QY_HXY[1:-1]>0.0,QY_HXY[1:-1]*VelUC_HX[:-1],QY_HXY[1:-1]*VelUC_HX[1:])
	VelUN_HX[:,1:-1]=(((QYU_HXY[:-1,1:-1]-QYU_HXY[1:,1:-1])+(TX_HXY[1:,1:-1]-TX_HXY[:-1,1:-1])/DensTol/dx)*dtdy+((QXU_N[:,:-1]-QXU_N[:,1:])+(PX_N[:,:-1]+DampPX_N[:,:-1]-PX_N[:,1:]-DampPX_N[:,1:])/DensTol/dy)*dtdx+(taugX_HX[:,1:-1]+tauWX_HX[:,1:-1])/DensTol*dt+HC_HX[:,1:-1]*VelUC_HX[:,1:-1])/HN_HX[:,1:-1]
	QY_N[:,:]=(QY_HY[:-1]+QY_HY[1:])/2
	QYV_N[:]=np.where(QY_N[:]>0.0,QY_N[:]*VelVC_HY[:-1],QY_N[:]*VelVC_HY[1:])
	QX_HXY[1:-1]=(QX_HX[:-1]+QX_HX[1:])/2
	QXV_HXY[:,1:-1]=np.where(QX_HXY[:,1:-1]>0.0,QX_HXY[:,1:-1]*VelVC_HY[:,:-1],QX_HXY[:,1:-1]*VelVC_HY[:,1:])
	VelVN_HY[1:-1,:]=(((QXV_HXY[1:-1,:-1]-QXV_HXY[1:-1,1:])+(TY_HXY[1:-1,1:]-TY_HXY[1:-1,:-1])/DensTol/dy)*dtdx+((QYV_N[:-1,:]-QYV_N[1:,:])+(PY_N[:-1,:]+DampPY_N[:-1,:]-PY_N[1:,:]-DampPY_N[1:,:])/DensTol/dx)*dtdy+(taugY_HY[1:-1,:]+tauWY_HY[1:-1,:])/DensTol*dt+HC_HY[1:-1,:]*VelVC_HY[1:-1,:])/HN_HY[1:-1,:]
	HC_N[:]=HN_N[:]
	HC_HX[:]=HN_HX[:]
	HC_HY[:]=HN_HY[:]
	VelUC_HX[:]=VelUN_HX[:]
	VelVC_HY[:]=VelVN_HY[:]
	#---------------UPDATE STRESS AND FORCE-------------
	VelUC_N[:,:]=(VelUC_HX[:,:-1]+VelUC_HX[:,1:])/2
	VelVC_N[:]=(VelVC_HY[:-1]+VelVC_HY[1:])/2
	tauW_N_Tri[:]=tauW_N[:]-GW*np.stack((VelUC_N*dt,VelVC_N*dt),axis=-1)
	temp=np.maximum(np.sqrt(np.einsum('ijk,ijk->ij', tauW_N_Tri,tauW_N_Tri)),np.ones((Ny,Nx))*0.00001,dtype='float32')
	PlaCorW_N=np.reshape(np.minimum(suW_N/temp,np.ones((Ny,Nx)),dtype='float32'),(Ny,Nx,1))
	tauW_N[:]=tauW_N_Tri[:]*PlaCorW_N[:]
	tauWX_HX[:,1:-1]=(tauW_N[:,:-1,0]+tauW_N[:,1:,0])/2
	tauWY_HY[1:-1,:]=(tauW_N[:-1,:,1]+tauW_N[1:,:,1])/2
	DUP_N=(1-np.reshape(PlaCorW_N,((Ny,Nx))))*np.sqrt(np.einsum('ijk,ijk->ij', tauW_N_Tri,tauW_N_Tri))/GW
	suW_N=np.where(((X_N[:,:]/Sof_X)**ns+(Y_N[:,:]/Sof_Y)**ns)<1,np.maximum(suW_N+SofSpeedW*DUP_N+SofTime*dt,TrW),np.maximum(suW_N+SofSpeedW*DUP_N,TrW))

	DEpsi_N[:,:,0]=(VelUC_HX[:,:-1]-VelUC_HX[:,1:])*dtdx  #COMPRESSION IS POSITIVE
	DEpsi_N[:,:,1]=(VelVC_HY[:-1,:]-VelVC_HY[1:,:])*dtdy
	DEpsi_N[:,:,2]=-DEpsi_N[:,:,0]-DEpsi_N[:,:,1]
	DEpsi_N[1:-1,1:-1,3]=((VelVC_N[1:-1,2:]-VelVC_N[1:-1,:-2])*dtdx/2+(VelUC_N[2:,1:-1]-VelUC_N[:-2,1:-1])*dtdy/2)/2
	SIJ_N[:]=SIJ_N[:]+2*GS*DEpsi_N[:]
	Deviaq_N_Tri[:,:]=np.sqrt((((SIJ_N[:,:,0]-SIJ_N[:,:,1])**2+(SIJ_N[:,:,1]-SIJ_N[:,:,2])**2+(SIJ_N[:,:,2]-SIJ_N[:,:,0])**2)/6+SIJ_N[:,:,3]**2)*3)
	temp=np.maximum(Deviaq_N_Tri,np.ones((Ny,Nx))*0.00001,dtype='float32')
	PlaCorS_N=np.reshape(np.minimum(2*suS_N/temp,np.ones((Ny,Nx)),dtype='float32'),(Ny,Nx,1))
	SIJ_N[:]=SIJ_N[:]*PlaCorS_N[:]
	DPEEQ_N=(1-np.reshape(PlaCorS_N,((Ny,Nx))))*Deviaq_N_Tri/(3*GS)/2
	PEEQ_N[:]=PEEQ_N[:]+DPEEQ_N[:]
	suS_N=np.maximum(suS_N+SofSpeedS*DPEEQ_N,TrS,dtype='float32')
	p_N[:,:,0]=DensNet*g*np.cos(Ang_N[:,:])*HC_N[:,:]/2-SIJ_N[:,:,2]
	Sigma_N[:]=SIJ_N[:]
	Sigma_N[:,:,:3]=SIJ_N[:,:,:3]+p_N[:,:]
	PX_N[:,:]=Sigma_N[:,:,0]*HC_N[:,:]*dy
	PY_N[:,:]=Sigma_N[:,:,1]*HC_N[:,:]*dx
	TX_HXY[1:-1,1:-1]=(Sigma_N[:-1,:-1,3]*HC_N[:-1,:-1]+Sigma_N[:-1,1:,3]*HC_N[:-1,1:]+Sigma_N[1:,:-1,3]*HC_N[1:,:-1]+Sigma_N[1:,1:,3]*HC_N[1:,1:])*dy/4
	TY_HXY[1:-1,1:-1]=(Sigma_N[:-1,:-1,3]*HC_N[:-1,:-1]+Sigma_N[:-1,1:,3]*HC_N[:-1,1:]+Sigma_N[1:,:-1,3]*HC_N[1:,:-1]+Sigma_N[1:,1:,3]*HC_N[1:,1:])*dx/4
	taugX_HX[:,1:-1]=DensNet*gx*HC_HX[:,1:-1]
	taugY_HY[1:-1,:]=DensNet*(-g)*np.sin(Ang_N[1:,:])*HC_HY[1:-1,:]
	DampPX_N[:,:]=-(VelUC_HX[:,1:]-VelUC_HX[:,:-1])*DampR*HC_N[:,:]*dy
	DampPY_N[:,:]=-(VelVC_HY[1:,:]-VelVC_HY[:-1,:])*DampR*HC_N[:,:]*dx
	#============output==============
	if j%output_No==0:
		SSA=0.0
		SSX=0
		SSY=0
		for x in np.nditer(suW_N):
			if x<TrW*1.1:
				SSA=SSA+1
		SSA=SSA*dx*dy
		for x in np.nditer(suW_N[:,int(Nx/2)]):
			if x<TrW*1.1:
				SSX=SSX+1
		SSX=(SSX+1)*dx
		for x in np.nditer(suW_N[int(Ny/2),:]):
			if x<TrW*1.1:
				SSY=SSY+1
		SSY=(SSY+1)*dy
		print (j*dt,datetime.datetime.now(),SSA,SSX,SSY)
		
		plt.figure(figsize=(8, 8))
		levels = np.linspace(4, 20, 9)
		CS = plt.contourf(X_N, Y_N, np.reshape(PlaCorS_N,((Ny,Nx)))*Deviaq_N_Tri,levels=levels, cmap='RdGy', extend='both')
		#plt.contourf(X_N, Y_N, np.reshape(PlaCorS_N,((Ny,Nx)))*Deviaq_N_Tri, cmap='RdGy')
		plt.colorbar(CS)
		fname='plotDQ'+str(j*dt*1000)+'.png'
		plt.savefig(fname,format='png', dpi=300)
		plt.close()
		plt.clf()
		
		plt.figure(figsize=(8, 8))
		levels = np.linspace(0.0, 0.2, 11)
		CS = plt.contourf(X_N, Y_N, PEEQ_N,levels=levels, cmap='RdBu', extend='both')
		plt.colorbar(CS)
		fname='plotPEQ'+str(j*dt*1000)+'.png'
		plt.savefig(fname,format='png', dpi=300)
		plt.close()
		plt.clf()

		
		plt.figure(figsize=(8, 8))
		levels = np.linspace(0.5, 1.5, 6)
		CS = plt.contourf(X_N, Y_N, HC_N/(HIni*np.cos(Ang_N[:,:])), levels=levels, cmap='RdBu', extend='both')
		plt.colorbar(CS)
		fname='plotH'+str(j*dt*1000)+'.png'
		plt.savefig(fname,format='png', dpi=300)
		plt.close()
		plt.clf()
		
		plt.figure(figsize=(8, 8))
		levels = np.linspace(3.0, 15.0, 7)
		CS = plt.contourf(X_N, Y_N, suW_N,levels=levels, cmap='RdBu', extend='both')
		plt.colorbar(CS)
		fname='plotSU'+str(j*dt*1000)+'.png'
		plt.savefig(fname,format='png', dpi=300)
		plt.close()
		plt.clf()

		
		









