import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy



n = int(input("please enter a Lattice dimensions:")) #Defining lattice dimensions 
errorlimit = float(input("Maximum error wanted:"))
distribution = str(raw_input("please choose type: CENTRAL or WIRE:"))
method = str(raw_input("please choose algorithm: Jacobi, Gauss or GaussRelaxation?:"))
distribution = distribution.upper()
method = method.upper()
psilattice  = np.zeros((n,n,n))
ELattice  = np.zeros((n,n,n))
BXLattice  = np.zeros((n,n,n))
BYLattice  = np.zeros((n,n,n))

if (distribution == "CENTRAL"):
    rholattice = np.zeros((n,n,n)) #Make charge distribution
    mid = int(n/2)
    rholattice[mid,mid,mid] = 1

if (distribution == "WIRE"):
    rholattice= np.zeros((n,n,n))
    for z in range(n):
       mid = int(n/2)
       rholattice[mid,mid,z] =1

                

def psi_update(psi_temp, rholattice,i,j,k):
    psi_new = (1/6)*(psi_temp[i+1,j,k] + psi_temp[i-1,j,k]+psi_temp[i,j+1,k] + psi_temp[i,j-1,k]+psi_temp[i,j,k+1] + psi_temp[i,j,k-1]+ rholattice[i,j,k])
    return psi_new

def psi_Gauss(psi_temp, psi_current, rholattice, i, j, k, omega):
    psi_new = (omega/6)*(psi_current[i-1, j, k] + psi_temp[i+1, j, k]+ psi_current[i, j-1, k]+psi_temp[i, j+1, k]+psi_current[i, j, k-1]+psi_temp[i, j, k+1]+rholattice[i,j,k])+(1-omega)*psi_current[i,j,k]
    return psi_new


def ELattice_update(psilattice, rholattice,dx,i,j,k):
    iterm =  ((psilattice[i+1,j,k]-psilattice[i-1,j,k])/(2*dx))**2
    jterm =  ((psilattice[i,j+1,k]-psilattice[i,j-1,k])/(2*dx))**2
    kterm =  ((psilattice[i,j,k+1]-psilattice[i,j,k-1])/(2*dx))**2
    delterm = (-1)*(np.sqrt(iterm+jterm+kterm))
    return delterm


#define B-field 
def BField_Update(ALattice,i,j,k,dx):
    ALattice = -1*ALattice #Differs from potential by a minus
    Bi = (ALattice[i,j+1,k]-ALattice[i,j-1,k])/(2.*dx)
    Bj = (ALattice[i-1,j,k]-ALattice[i+1,j,k])/(2.*dx)
    #norm = np.sqrt(Bi**2+Bj**2)
    return Bi,Bj
    

dt = 0.5
dx = 1
dy = 1
dx = 1
eps_0 = 1
numsteps =3000
T = numsteps*dt

times = np.linspace(0,T,numsteps)   

t=0
psierror_total = 3
relax_error = 1
omega = 2/(1+np.sin((np.pi/n+1)))

n_iterations = []
errorlist=[]



#Jacobi algorithm

if (method == "JACOBI"): 
    while ((t < numsteps) & (psierror_total > errorlimit)): 
        print(t)
        n_iterations.append(t)
        psi_temp = np.copy(psilattice)
        psierror = 0
        for i in range(1, (n-1)):
            for j in range(1,(n-1)):
                for k in range(1,(n-1)):
                    psilattice[i,j,k] = psi_update(psi_temp, rholattice,i,j,k) # A lattice is actually same as psi but *-1 for this purpose and I = Q =1
                    psierror += abs(psilattice[i,j,k]-psi_temp[i,j,k])

        psierror_total = psierror
        print(psierror_total)
        errorlist.append(psierror_total)
        #print(psierror_total)
        if (t == numsteps-1):
            print("Error did not reach anticipated value")
        t += 1


plt.scatter(n_iterations, errorlist)
plt.show() 
#==============================================================================       
#gauss algorithm
omega = 0.5

if (method == "GAUSS"): 
    while ((t < numsteps) & (psierror_total > errorlimit)): 
        print(t)
        toter = 0
        psi_temp = np.copy(psilattice)
        for i in range(1, (n-1)):
            for j in range(1,(n-1)):
                for k in range(1,(n-1)):
                    #psiTemp = 
                    psilattice[i,j,k] = psi_Gauss(psi_temp,psilattice, rholattice,i,j,k, omega) 
                    
                    psierror = abs(psilattice[i,j,k]-psi_temp[i,j,k])
                    toter += psierror                   
        psierror_total = toter 
        errorlist.append(psierror_total)
        print(psierror_total)
        if (t == numsteps-1):
            print("Error did not reach anticipated value")
        t += 1

print("Final error was {}".format(psierror_total))

#==============================================================================
#relaxation method

iterationVals = []


if (method == 'GAUSSRELAXATION'):
    omegaRange = np.arange(0, 2., 0.1)
    for f in range(len(omegaRange)):
        t=0
        w  = np.zeros((n,n,n))
        
        psierror_total = 2*errorlimit
        while ((psierror_total > errorlimit)): 
            print(t)
            toter = 0
            psi_temp = np.copy(w)
            for i in range(1, (n-1)):
                for j in range(1,(n-1)):
                    for k in range(1,(n-1)):
                        w[i,j,k] = psi_Gauss(psi_temp,w, rholattice,i,j,k, omegaRange[f]) 
                        normFactor = float(n**3)
                        psierror = abs(psilattice[i,j,k]-psi_temp[i,j,k])
                        psierror = psierror/normFactor
                        toter += psierror                   
            psierror_total = toter 
            errorlist.append(psierror_total)
            print(psierror_total)
            #if (t == numsteps-1):
                #print("Error did not reach anticipated value")
            t = t+1
        iterationVals.append(t)
        #print(iterationVals)
    plt.plot(omegaRange, iterationVals)
    plt.show()




#Create Measurment lattices
for i in range(1,(n-1)):
    for j in range(1,(n-1)):
        for k in range(1,(n-1)):
            ELattice[i,j,k] = ELattice_update(psilattice, rholattice,dx,i,j,k)
            #print(ELattice[i,j,k])
            BXLattice[i,j,k], BYLattice[i,j,k] = BField_Update(psilattice,i,j,k,dx) 



oneD = []
half = int(n/2)


rhoLattice = np.copy(rholattice)
values_rho = rhoLattice[:][:][half]

#contour plot of the potential strength in space /cut through the mid-plane

plt.figure()
plt.imshow(psilattice[half][:][:])
plt.colorbar()
plt.show()



#E-field vector plot using quiver
fieldCalc = np.gradient(psilattice)
fieldCalc = np.asarray(fieldCalc)
#print(fieldCalc.shape)
u = -1.*fieldCalc[0]
v = -1*fieldCalc[1]
q = -1*fieldCalc[2]
plt.figure()
plt.quiver(u[:][:][half], v[:][:][half], angles = 'xy')
plt.show()


UE_field = np.arange(0, len(fieldCalc[0]), 1)
VE_field = np.arange(0, len(fieldCalc[0]), 1)
QE_field = np.arange(0, len(fieldCalc[0]), 1)
Xe,Ye,Ze = np.meshgrid(UE_field,VE_field,QE_field,indexing = 'ij')

"""
# Magnetic Field Calc
"""

def GetBField(psilattice):
        bi, bj = np.gradient(psilattice[:][:][1])
        B_field=np.array([-bi,bj,np.zeros((n,n))])
        normB=np.sqrt(B_field[0]**2 + B_field[1]**2)
        B_fieldnorm=B_field/normB
        return B_fieldnorm
 
B_field = GetBField(psilattice)   

plt.figure()
plt.imshow(psilattice[:][:][1])

plt.quiver(B_field[0],B_field[1],angles='xy')
plt.xlabel('$B_{x}$ Field')
plt.ylabel('$B_{y}$ Field')
plt.title('Magnetic Field Plot')
plt.show()

"""
#Magnetic Field Plotting
"""

#B-field potential vs distance from charge
fileL = open("B-potential.txt","w")
Bdistance_values=[]
Bpotential_values = []
for i in range(n):
    for j in range(n):
            distance = np.sqrt(((i-half)**2)+((j-half)**2))
            Bdistance_values.append(np.log(distance))
            pot = psilattice[i,j,half] #fixed k to middle of wire
            Bpotential_values.append(np.log(pot))
            fileL.write('{} {}\n'.format(Bdistance_values,Bpotential_values))
    fileL.write('\n')
        
fileL.close()


plt.figure()
plt.scatter(Bdistance_values, Bpotential_values)
plt.xlabel('Distance from charge')
plt.ylabel('B-field Potential Values')
plt.title('B-field potential vs distance')
plt.show()



#B-field vs distance
dd_values = []
ef_values = []
fileZ = open("E-strength.txt","w")
for i in range(len(Xe)):
    for j in range(len(Ye)):
        for k in range(len(Ze)):
            dist = np.sqrt(((i-half)**2)+((j-half)**2)+((k-half)**2))
            dd_values.append(np.log(dist))
            f = np.sqrt(fieldCalc[0,i,j,k]**2 + fieldCalc[1,i,j,k]**2 + fieldCalc[2,i,j,k]**2)
            ef_values.append(np.log(f))
            fileZ.write('{} {}\n'.format(dd_values,ef_values))
        fileZ.write('\n')
fileZ.close()

plt.figure()
plt.scatter(dd_values, ef_values)
plt.xlabel('Distance from charge log(1/r)')
plt.ylabel('E-field Values')
plt.show()

#==============================================================================
#==============================================================================
#=========================  Electric field Plotting============================
#==============================================================================
#E-field potential vs distance
distance_values = []
potential_values = []

fileR = open("E-potential.txt","w")
for i in range(n):
    for j in range(n):
        for k in range(n):
            distance = np.sqrt(((i-half)**2)+((j-half)**2)+((k-half)**2))
            #print((distance))
            distance_values.append(np.log(distance))
            pot = psilattice[i,j,k]
            #print(((pot)))
            potential_values.append(np.log(pot))
            fileR.write('{} {}\n'.format(distance_values,potential_values))
        fileR.write('\n')
fileR.close()
#print(potential_values)


plotting = np.polyfit(distance_values, potential_values, 1)
print(plotting)
plt.figure()
plt.scatter(distance_values, potential_values)
plt.xlabel('Distance from charge')
plt.ylabel('Potential Values')
plt.title('E-field potential vs distance')

plt.show()


###############################################################################
#E-field vs distance
dd_values = []
ef_values = []
fileZ = open("E-strength.txt","w")
for i in range(len(Xe)):
    for j in range(len(Ye)):
        for k in range(len(Ze)):
            dist = np.sqrt(((i-half)**2)+((j-half)**2)+((k-half)**2))
            dd_values.append(np.log(dist))
            f = np.sqrt(fieldCalc[0,i,j,k]**2 + fieldCalc[1,i,j,k]**2 + fieldCalc[2,i,j,k]**2)
            ef_values.append(np.log(f))
            fileZ.write('{} {}\n'.format(dd_values,ef_values))
        fileZ.write('\n')
fileZ.close()

plt.figure()
plt.scatter(dd_values, ef_values)
plt.xlabel('Distance from charge log(1/r)')
plt.ylabel('E-field Values')
plt.show()


#============================================================================



    
    


