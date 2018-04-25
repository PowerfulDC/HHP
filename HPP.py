#!/usr/bin/python
import numpy as np
from subprocess import call
import os
import math
import shutil
import sys
import argparse 

from eos.eos import Eos
#################################################
#the grids of the  hydro profile:
#   NX0=NY0=70
#	NZ = 41
#	DX=DY=DZ=0.3
#	TAU0 = 0.6
#	DT = 0.3
#################################################
class HPP(object):
	def __init__(self,path):
		data_path = path
		print('Loading data file,please wait for a minuts!')
		datafile = os.path.join(data_path,'bulk3D.dat')
		self.data_t = np.loadtxt(datafile)
		print('Data file loading complete!')
		self.NX0 = 70
		self.NY0 = 70
		self.NZ0 = 41
		self.TAU0 = 0.6
		self.DX0 = 0.3
		self.DY0 = 0.3
		self.DZ0 = 0.3
		self.DT = 0.3
		self.NT = self.data_t.shape[0]//(self.NX0*self.NY0*self.NZ0)
		print("steps of Time is %i"%self.NT)
		self.NX = self.NX0 
		self.NY = self.NY0
		self.NZ = self.NZ0
		self.DX = self.DX0	
		self.DY = self.DY0
		self.DZ = self.DZ0
			
		#Switchs		
		self.Dim_Switch = 3
		self.Grids_Switch = False
		self.sEd = True
		self.sT = False
		self.sVt = True
		self.sVz = True
		self.sFrac = False
			
		self.IEOS = 1
		self.eos = Eos(self.IEOS)

		self.OutPutPath = None

		# self.Ed_txyz = self.data_t[:,0].reshape(self.NT,self.NX0,self.NY0,self.NZ0)
		# self.Vx_txyz = self.data_t[:,1].reshape(self.NT,self.NX0,self.NY0,self.NZ0)
		# self.Vy_txyz = self.data_t[:,2].reshape(self.NT,self.NX0,self.NY0,self.NZ0)
		# self.Vz_txyz = self.data_t[:,3].reshape(self.NT,self.NX0,self.NY0,self.NZ0)
		self.Block_txyz = np.zeros(self.NT*self.NX0*self.NY0*self.NZ0*4).reshape(self.NT,self.NX0,self.NY0,self.NZ0,4)
		self.Block_txyz[:,:,:,:,0] = self.data_t[:,0].reshape(self.NT,self.NX0,self.NY0,self.NZ0) #Ed
		self.Block_txyz[:,:,:,:,1] = self.data_t[:,1].reshape(self.NT,self.NX0,self.NY0,self.NZ0) #Vx
		self.Block_txyz[:,:,:,:,2] = self.data_t[:,2].reshape(self.NT,self.NX0,self.NY0,self.NZ0) #Vy
		self.Block_txyz[:,:,:,:,3] = self.data_t[:,3].reshape(self.NT,self.NX0,self.NY0,self.NZ0) #Vz

#		self.OutPut_col_shape = []
		
		self.Ed_newGrids = []
		self.T_newGrids = []
		self.Frac_newGrids = []
		self.Vt_newGrids = []
		self.Vz_newGrids = []
		self.Grids = []
		self.Hotel = []
		self.todo = []

# load customized grids
	def FormatCFG(self,NX = 200, NY = 200 , NZ= 200, DeltX = 0.3, DeltY = 0.3, DeltZ = 0.3, \
		Dim_Switch = 3, Grids_Switch = False,sEd = True,sT = False, sFrac = False, sV = True, Outputpath = None):
		self.NX = NX
		self.NY = NY
		self.NZ = NZ
		self.DX = DeltX
		self.DY = DeltY
		self.DZ = DeltZ
		self.Dim_Switch = Dim_Switch
		self.Grids_Switch =Grids_Switch
		self.sEd = sEd
		self.sT = sT
		self.sFrac = sFrac
		self.sVt = sV
		self.sVz = sV
		self.OutPutPath = Outputpath
		if Dim_Switch == 3:
			# self.OutPut_col_shape = np.zeros(self.NT*NX*NY*NZ)
			#self.Vt_newGrids = np.zeros((self.NT*NX*NY*NZ,2))
			self.Hotel = np.zeros((self.NT*self.NX*self.NY*self.NZ,10))
#			if Grids_Switch:
#			    self.Grids = np.zeros((self.NT*NX*NY*NZ,4))
		elif Dim_Switch == 2:
			self.sVz = False
                       # self.OutPut_col_shape = np.zeros(self.NT*NX*NY)
                     #   self.Vt_newGrids = np.zeros((self.NT*NX*NY,2))
			self.Hotel = np.zeros((self.NT*self.NX*self.NY,10))
                        #if Grids_Switch:
                       #     self.Grids = np.zeros((self.NT*NX*NY,3))

		
#give the final time
	def Finaltimestep(self):
		return self.NT,self.TAU0,self.DT

#give the fraction of QGP and Hadron
	def Frac(self,Temp):
		if Temp >0.22:
			frac = 1.0
		elif Temp <0.18:
			frac = 0.0
		else:
			frac = (Temp-0.18)/(0.22-0.18)
		return frac
	
#change 2D to 3D directly
	def change3Dto2D(self):
		OutPut2D = self.Block_txyz[:,:,:,self.NZ0//2,0:3].reshape(self.NT*self.NX0*self.NY0,3)
		#np.savetxt(filepath+'new_bulk2D.dat',OutPut2D,\
                #header = 'ED,T,frac,VX,VY'+' NT=%i'%self.NT+' NX=%i'%self.NX0+' NY=%i'%self.NY0)
		return OutPut2D


#return the hydro imformation of the input location
	def loc(self,t=0.6,x=0,y=0,z=0):
		# (x,y,z)should loacted in the range of (xmin,xmax)&(ymin,ymax)and so on
		# This peculiar part is because the trento's grid cannot get the 0 points at the transverse plane
		X_min = -self.DX*(self.NX0-1)/2.0
		X_max = -X_min
		Y_min = -self.DY*(self.NY0-1)/2.0
		Y_max = -Y_min

		Z_min = -(self.NZ0//2)*self.DZ
		Z_max = -Z_min
		if x>X_max or x<X_min or y>Y_max or y<Y_min or z>Z_max or z<Z_min:
			return np.zeros(4)
		else:	
			L_NX_L = int((x-X_min)/self.DX)
			L_NY_L = int((y-Y_min)/self.DY)
			L_NZ_L = int((z-Z_min)/self.DZ)
			L_NT_L = int((t-self.TAU0)/self.DT)
			rt = abs((t-self.TAU0)/self.DT - L_NT_L)
			if rt<1e-6:
				L_NT = L_NT_L
			elif (rt-1)<1e-6:
				L_NT = L_NT_L+1

			#bi mian ge dian shang de quzhi dao zhi fushu
			rx = abs((x-X_min)/self.DX-L_NX_L)
			ry = abs((y-Y_min)/self.DY-L_NY_L)
			rz = abs((z-Z_min)/self.DZ-L_NZ_L)

			return self.Int3D(rx,ry,rz,L_NX_L,L_NY_L,L_NZ_L,L_NT)
#2D

        def loc2D(self,t=0.6,x=0,y=0):
                # (x,y,z)should loacted in the range of (xmin,xmax)&(ymin,ymax)and so on
                # This peculiar part is because the trento's grid cannot get the 0 points at the transverse plane
                X_min = -self.DX*(self.NX0-1)/2.0
                X_max = -X_min
                Y_min = -self.DY*(self.NY0-1)/2.0
                Y_max = -Y_min

                if x>X_max or x<X_min or y>Y_max or y<Y_min :
                        return np.zeros(3)
                else:
                        L_NX_L = int((x-X_min)/self.DX)
                        L_NY_L = int((y-Y_min)/self.DY)
                        L_NT_L = int((t-self.TAU0)/self.DT)
                        rt = abs((t-self.TAU0)/self.DT - L_NT_L)
                        if rt<1e-6:
                                L_NT = L_NT_L
                        elif (rt-1)<1e-6:
                                L_NT = L_NT_L+1

                        #bi mian ge dian shang de quzhi dao zhi fushu
                        rx = abs((x-X_min)/self.DX-L_NX_L)
                        ry = abs((y-Y_min)/self.DY-L_NY_L)

                        return self.Int2D(rx,ry,L_NX_L,L_NY_L,L_NT)

#make 3D chazhi for different observer
	def Int3D(self,rx,ry,rz,L_NX_L,L_NY_L,L_NZ_L,L_NT):
		int100 = self.Block_txyz[L_NT,L_NX_L,L_NY_L,L_NZ_L,:]*(1-rx)+self.Block_txyz[L_NT,L_NX_L+1,L_NY_L,L_NZ_L,:]*rx
		int101 = self.Block_txyz[L_NT,L_NX_L,L_NY_L,L_NZ_L+1,:]*(1-rx)+self.Block_txyz[L_NT,L_NX_L+1,L_NY_L,L_NZ_L+1,:]*rx
		int110 = self.Block_txyz[L_NT,L_NX_L,L_NY_L+1,L_NZ_L,:]*(1-rx)+self.Block_txyz[L_NT,L_NX_L+1,L_NY_L+1,L_NZ_L,:]*rx
		int111 = self.Block_txyz[L_NT,L_NX_L,L_NY_L+1,L_NZ_L+1,:]*(1-rx)+self.Block_txyz[L_NT,L_NX_L+1,L_NY_L+1,L_NZ_L+1,:]*rx
		intA = int101*rz + int100*(1-rz)
		intB = int111*rz + int110*(1-rz)
		intF = intB*ry + intA*(1-ry)

		return intF # intF[0]=Ed , intF[1]=Vx .........
#2D
        def Int2D(self,rx,ry,rz,L_NX_L,L_NY_L,L_NT):
		int10 =  self.Block_txyz[L_NT,L_NX_L,L_NY_L,self.NZ0//2,:]*(1-rx)+self.Block_txyz[L_NT,L_NX_L+1,L_NY_L,self.NZ0//2,:]
		int11 =  self.Block_txyz[L_NT,L_NX_L,L_NY_L+1,self.NZ0//2,:]*(1-rx)+self.Block_txyz[L_NT,L_NX_L+1,L_NY_L+1,self.NZ0//2,:]
		intF2D = int11*ry +int10*(1-ry)
		return intF2D[0:3]

# save data file with required format
	def save(self):#, sGrids = False, sEd = True, sT = False, sFrac = False, sVt = True, sVz = True):
		m=0
		qmark = np.zeros(10,bool)
		sGrids_t = False
		sGrids_z = False
		if self.Grids_Switch:
			sGrids_t = True
			sGrids_z = True
		if self.Dim_Switch == 2:
			sGrids_z = False
			self.sVz = False
		for quant in (sGrids_t,sGrids_t,sGrids_t,sGrids_z, self.sEd , self.sT, self.sFrac,self.sVt,self.sVt,self.sVz): #(t,x,y,z,ed,T,frac,vx,vy,vz)
			qmark[m] = quant 
			m+=1
		if self.Grids_Switch:		
		        if self.NX%2 == 1:
        		    xline = np.linspace(-np.floor(self.NX/2)*self.DX, np.floor(self.NX/2)*self.DX, self.NX, endpoint=True)
            		    yline = np.linspace(-np.floor(self.NY/2)*self.DY, np.floor(self.NY/2)*self.DY, self.NY, endpoint=True)
                   	    
	      		elif self.NX%2 == 0:
	            	    xline = np.linspace(-((self.NX-1)/2.0)*self.DX, ((self.NX-1)/2.0)*self.DX, self.NX, endpoint=True)
       		    	    yline = np.linspace(-((self.NY-1)/2.0)*self.DY, ((self.NY-1)/2.0)*self.DY, self.NY, endpoint=True)
			
			tau = np.linspace(self.TAU0, self.TAU0+(self.NT-1)*self.DT,self.NT,endpoint=True)
			print(tau.shape)
			x_t = np.tile(xline,self.NT)
			y_tx = np.tile(yline,self.NT*self.NX)
			if self.Dim_Switch == 2:
				self.Hotel[:,0]= np.repeat(tau,self.NX*self.NY)
				self.Hotel[:,1] = np.repeat(x_t,self.NY)
				self.Hotel[:,2] = y_tx
			if self.Dim_Switch == 3:
				zline = np.linspace(-np.floor(self.NZ/2)*self.DZ, np.floor(self.NZ/2)*self.DZ, self.NZ, endpoint=True)
				self.Hotel[:,0] = np.repeat(tau,blocksize/self.NT)
				self.Hotel[:,1] = np.repeat(x_t,self.NY*self.NZ)
				self.Hotel[:,2] = np.repeat(y_tx,self.NZ)
				self.Hotel[:,3] = np.tile(zline,blocksize/self.NZ)
		if self.sEd:
			self.Hotel[:,4] = self.todo[:,0]
		if self.sT:
			self.Hotel[:,5] =self.eos.f_T(self.Hotel[:,4])   #T
		if self.sFrac:
			self.Hotel[:,6] = np.array(map(self.Frac,self.Hotel[:,5]))  #Frac
		if self.sVt:
			self.Hotel[:,7] = self.todo[:,1]
			self.Hotel[:,8] = self.todo[:,2]
		if self.sVz:
			self.Hotel[:,9] = self.todo[:,3]
			
		OutPutData = self.Hotel[:,qmark]
		os.chdir(self.OutPutPath)
		if Dim_Switch == 3:
                	np.savetxt('new_bulk3D.dat',OutPutData,\
                	header = 'ED,VX,VY,VEta'+' NT=%i '%self.NT+' NX=%i'%self.NX+' NY=%i'%self.NY+' NZ=%i'%self.NZ)
			print("new_bulk3D.dat Finished")
                elif Dim_Switch == 2:
                	np.savetxt('new_bulk2D.dat',OutPutData,\
               		header = 'ED,VX,VY'+' NT=%i'%self.NT+' NX=%i'%self.NX+' NY=%i'%self.NY)
			print("new_bulk2D.dat Finished")

		
		

	# def load_hydro_profile(data_path)
	# 	data_t = np.loadtxt('bulk3D.dat')
	# 	# Ed,vx,vy,vz
	# 	return data_t[:,0],data_t[:,1],data_t[:,2],data_t[:,3]


def main(data_root_path, Nx,Ny, Nz, Deltx, Delty,Deltz,Event_number, Dim_switch,Grids_switch,sed,st,sfrac,sv,outputpath):
	OutPutPath = outputpath
	if OutPutPath is None:
		OutPutPath  = data_root_path
	os.chdir(OutPutPath)
	OutPutPath = os.getcwd()
		
	for i in xrange(0,Event_number):
		if Event_number > 1:
       	        	 fpath='event%d'%i
       			 filepath = os.path.join(data_root_path,fpath)
       			 OutPutPath = os.path.join(OutPutPath,fpath)       			 
    		elif Event_number == 1:
    			filepath = data_root_path
   		else :
    			print('Event_number shoule be a positive int number')
    			exit(0)
    		if not os.path.exists(OutPutPath):
        	    os.makedirs(OutPutPath)

	    	print('Initialized!!!!')
	  	doit = HPP(path=filepath)
		
		doit.FormatCFG(NX=Nx,NY=Ny,NZ=Nz,DeltX=Deltx,DeltY=Delty, \
		DeltZ=Deltz,Dim_Switch=Dim_switch,Grids_Switch=Grids_switch,\
		sEd=sed,sT=st,sFrac=sfrac,sV=sv,Outputpath=OutPutPath)


	    	NT,TAU0,DeltT = doit.Finaltimestep()
		if Dim_Switch == 2 and NX == doit.NX0 and NY ==doit.NY0 \
		   and DeltX == doit.DX0 and DeltY == doit.DY0:
			print('Start to change file directly!')
		   	doit.todo = doit.change3Dto2D()
		   	doit.save()
		   	continue
		if Dim_Switch ==3 and NX == doit.NX0 and NY ==doit.NY0 and NZ == self.NZ0 \
		   and DeltX == doit.DX0 and DeltY == doit.DY0 and DeltZ == self.DeltZ:
			print('Start to add colums to file directly!')
			doit.todo = doit.Block_txyz.reshape(self.NT*self.NX0*self.NY0*self.NZ0,4)
			doit.save()
			continue

    		Fdata3D = np.zeros(NT*NX*NY*NZ*4).reshape(NT,NX,NY,NZ,4)
		Fdata2D = np.zeros(NT*NX*NY*4).reshape(NT,NX,NY,4)

		print('start to calculate New Hydro profile!')
	    	for ti in xrange(0,NT):
	    		tau = TAU0+ti*DeltT
			print('NTstep:%i'%ti)
    			for xi in xrange(0,NX):
				if NX%2 == 1:
    					x = (xi-NX//2)*DeltX
				elif NX%2 == 0:
					x = (xi-(NX-1)/2.0)*DeltX
    		 		for yi in xrange(0,NY):
					if NY%2 == 1:
        	                	  	y = (yi-NY//2)*DeltY
	       	     	        	elif NY%2 == 0:
                       			 	y = (yi-(NY-1)/2.0)*DeltY
					if Dim_Switch == 3:
    						for zi in xrange(0,NZ):
    							z = (zi-NZ//2)*DeltZ
    							Fdata3D[ti,xi,yi,zi,:]=doit.loc(tau,x,y,z)
					elif Dim_Switch == 2:
						Fdata2D[ti,xi,yi,:]=doit.loc2D(tau,x,y,)
		
		
		if Dim_Switch == 3:
        	    doit.todo=Fdata3D.reshape(NT*NX*NY*NZ,4)
 	      	elif Dim_Switch == 2:
       		    doit.todo=Fdata2D.reshape(NT*NX*NY,3) 
		doit.save()
	print('Done!')
		
#		if Dim_Switch == 3:
#	    		np.savetxt(filepath+'new_bulk3D.dat',Fdata3D.reshape(NT*NX*NY*NZ,4),\
#			header = 'ED,VX,VY,VEta'+' NT=%i '%NT+' NX=%i'%NX+' NY=%i'%NY+' NZ=%i'%NZ)
#		elif Dim_Switch == 2:
#			np.savetxt(filepath+'new_bulk2D.dat',Fdata2D.reshape(NT*NX*NY,0:3),\
#			header = 'ED,VX,VY'+' NT=%i'%NT+' NX=%i'%NX+' NY=%i'%NY)





if __name__ == '__main__':
    print('This program can only be used for the hydro profiles with the TRENTO initial condition.')
    # if len(sys.argv) != 10:
    # 	print('usage:python HPP.py PATH NX NY NZ DeltX DeltY DeltZ Event_number Dim_Switch')
   	# 	print('PATH is the root path of hydro profiles')
   	# 	print("Dim_Switch if a switch for 2D or 3D ,please use '2' or '3' ")
    # 	exit(0)
    parser = argparse.ArgumentParser()

    parser.add_argument("--path",required =True, help="Give the root path for event by event;Give the file pathe for 1 event")
    parser.add_argument("--NX", nargs='?', type = int, \
    default = 70, const= 70,help = "NX, default = 70")
    parser.add_argument("--NY", nargs='?', type = int, \
    default = 70, const= 70,help = "NY, default = 70")
    parser.add_argument("--NZ", nargs='?', type = int, \
    default = 41, const= 41,help = "NZ, default = 41")
    parser.add_argument("--DX", nargs='?', type = float, \
    default = 0.3, const= 0.3,help = "DX, default = 0.3")
    parser.add_argument("--DY", nargs='?', type = float, \
    default = 0.3, const= 0.3,help = "DY, default = 0.3")
    parser.add_argument("--DZ", nargs='?', type = float, \
    default = 0.3, const= 0.3,help = "DZ, default = 0.3")
    parser.add_argument('-NE',"--EBE_N", nargs='?', type = int, \
    default = 1, const= 1,help = "Event number, default = 1")
    parser.add_argument("-dim","--Dim_Switch", nargs='?', type = int, \
    default = 3, const= 3,help = "Dim_Switch, default = 3,should be '2' or '3'")
    parser.add_argument('-g',"--Grids_Switch", nargs='?', type = bool, \
    default = False, const= False,help = "Grids, default = False") 
    parser.add_argument('-ed',"--sEd", nargs='?', type = bool, \
    default = True, const= True,help = "Energy, default = True") 
    parser.add_argument('-T',"--sT", nargs='?', type = bool, \
    default = False, const= False,help = "Temperature, default = False")
    parser.add_argument('-fqgp',"--sFrac", nargs='?', type = bool, \
    default = False, const= False,help = "Fraction of QGP, default = False")
    parser.add_argument('-v',"--sV", nargs='?', type = bool, \
    default = True, const= True,help = "Velocity of flow, default = True")
    parser.add_argument('-o',"--OutPutPath", nargs='?', \
    default = None, const= None,help = "Out Put Path, default = The same with the origin Data ")

    cfg=parser.parse_args()


    data_root_path = cfg.path
    OutPutPath = cfg.OutPutPath
    NX = cfg.NX
    NY = cfg.NY
    NZ = cfg.NZ
    DeltX = cfg.DX
    DeltY = cfg.DY
    DeltZ = cfg.DZ
    Event_Number = cfg.EBE_N
    Dim_Switch = cfg.Dim_Switch
    Grids_Switch =cfg.Grids_Switch
    Switch_Ed = cfg.sEd
    Switch_T = cfg.sT
    Switch_Fqgp = cfg.sFrac
    Switch_V = cfg.sV
    if Dim_Switch == 3 or Dim_Switch == 2:
	    main(data_root_path, Nx=NX, Ny=NY, Nz=NZ, Deltx=DeltX, Delty=DeltY, Deltz=DeltZ,\
	    Event_number=Event_Number, Dim_switch=Dim_Switch, Grids_switch=Grids_Switch,sed=Switch_Ed,\
	    st=Switch_T,sfrac=Switch_Fqgp,sv=Switch_V,outputpath=OutPutPath)
    else:
	print("Dim_Switch if a switch for 2D or 3D ,please use '2' or '3' ")
	exit(0)


