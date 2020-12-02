#Importation des bibliothèques

import math as ma
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import plotly
import numpy as np
from numpy import linalg
import plotly.graph_objects as go
import plotly.offline as py
import json
import xlrd
from xlwt import Workbook, Formula
import xlsxwriter, pandas
import time

# ouverture du fichier Excel 

wb = xlrd.open_workbook('SIGMA_MECA.xlsx')
 
# lecture des données dans chaque feuille

RAYON= wb.sheet_by_name(u'SIGMA_MECA')

nb_ligne=RAYON.nrows
nb_colonne=RAYON.ncols

Rayon=np.zeros((nb_ligne,nb_colonne))
for i in range(nb_colonne):
    Rayon[:,i]=RAYON.col_values(i)

#Paramètres géométriques

Rint=Rayon[-1,1] #Rayon intérieur
Rext=50 #Rayon extérieur 

#Paramètres temporels et discrétisation 

time=50 #durée
N_time=500 #STEP de calcul
delta_t=time/N_time #pas de temps
nr=6 #Nombre d'éléments dans l'épaisseur en 2D
na=120 #Nombre d'éléments circonférentiel en 2D

#Paramètres matériaux et physiques

rho=1.225/1000000000 #Masse volumique [kg/mm^3]
lamda=0.025/1000 #Conductivité thermique [W/(mm*K)]
cp=1005 #Capacité thermique [J/Kg*K]
alpha=0/10**6 #Dilatation thermique
q=0 #Flux de chaleur
E=210e3 #Module d'Young
nu=0.3 #Coefficient de Poisson
To=20+273 #Température homogène à l'intant initial [K]
Tf=70+273 #Température finale de chauffe [K]
P=150 #Pression
F_l=P*2*ma.pi*Rint/na #Force linéique appliqué aux noeuds

##Calcul des différentes matrices élémentaires

#Calcul de J Jacobien de la transformation et B le gradient aisni que les fonctions de formes N

def JBN(xi,eta,x1,y1,x2,y2,x3,y3,x4,y4):
    J=0.25*np.dot(np.array([[eta-1,-eta-1,1+eta,1-eta],[xi-1,1-xi,1+xi,-xi-1]]),np.array([[x1,y1],[x2,y2],[x3,y3],[x4,y4]]))
    I_J=np.linalg.inv(J)
    B=0.25*np.dot(I_J,np.array([[eta-1,-eta-1,1+eta,1-eta],[xi-1,1-xi,1+xi,-xi-1]]))
    N=np.array([[0.25*(1-xi)*(1-eta),0.25*(1-xi)*(1+eta),0.25*(1+xi)*(1+eta),0.25*(1+xi)*(1-eta)]])
    return J,B,N

#Pour un élément à 4 noeuds

def MAT_CAPA_THER(x1,y1,x2,y2,x3,y3,x4,y4): #Matrice de capacité ou masse thermique 
    w_k=[1,1,1,1]
    x_k=[-ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3),ma.sqrt(1/3)]
    y_k=[-ma.sqrt(1/3),ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3)]
    C=np.zeros(4)
    for i in range(len(w_k)):
        for j in range(len(w_k)):
            J,B,N=JBN(x_k[i],y_k[j],x1,y1,x2,y2,x3,y3,x4,y4)
            C=C+rho*cp*np.linalg.det(J)*np.dot(np.transpose(N),N)
    return C

def MAT_RIGI_THER(x1,y1,x2,y2,x3,y3,x4,y4): #Matrice de conductivité ou rigidité thermique 
    w_k=[1,1,1,1]
    x_k=[-ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3),ma.sqrt(1/3)]
    y_k=[-ma.sqrt(1/3),ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3)]
    K=np.zeros(4)
    for i in range(len(w_k)):
        for j in range(len(w_k)):
            J,B,N=JBN(x_k[i],y_k[j],x1,y1,x2,y2,x3,y3,x4,y4)
            K=K+lamda*np.linalg.det(J)*np.dot(np.transpose(B),B)
    return K
    
def MAT_RIGI_MECA(x1,y1,x2,y2,x3,y3,x4,y4): #Matrice de rigidité mécanique 
    w_k=[1,1,1,1]
    x_k=[-ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3),ma.sqrt(1/3)]
    y_k=[-ma.sqrt(1/3),ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3)]
    K=np.zeros(8)
    Q=np.zeros((3,4))
    D=(E/(1-nu**2))*np.array([[1,nu,0],[nu,1,0],[0,0,0.5*(1-nu)]])
    for i in range(len(w_k)):
        for j in range(len(w_k)):
            J,B,N=JBN(x_k[i],y_k[j],x1,y1,x2,y2,x3,y3,x4,y4)
            Q[0,0]=(1/np.linalg.det(J))*J[1,1]
            Q[0,1]=-(1/np.linalg.det(J))*J[0,1]
            Q[1,2]=-(1/np.linalg.det(J))*J[1,0]
            Q[1,3]=(1/np.linalg.det(J))*J[0,0]
            Q[2,0]=-(1/np.linalg.det(J))*J[1,0]
            Q[2,1]=(1/np.linalg.det(J))*J[0,0]
            Q[2,2]=(1/np.linalg.det(J))*J[1,1]
            Q[2,3]=-(1/np.linalg.det(J))*J[0,1]
            B_tilde=np.array([[B[0,0],0,B[0,1],0,B[0,2],0,B[0,3],0],[B[1,0],0,B[1,1],0,B[1,2],0,B[1,3],0],[0,B[0,0],0,B[0,1],0,B[0,2],0,B[0,3]],[0,B[1,0],0,B[1,1],0,B[1,2],0,B[1,3]]])
            B_calc=np.dot(Q,B_tilde)
            K=K+np.linalg.det(J)*np.dot(np.dot(np.transpose(B_calc),D),B_calc)
    return K

##Calcul des différents vecteurs élémentaires

#Pour un élément à quatre noeuds

def VEC_FLUX_NOD(x1,y1,x2,y2,x3,y3,x4,y4): #Vecteur des flux nodaux 
    w_k=[1,1,1,1]
    x_k=[-ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3),ma.sqrt(1/3)]
    y_k=[-ma.sqrt(1/3),ma.sqrt(1/3),-ma.sqrt(1/3),ma.sqrt(1/3)]
    F=np.zeros((4,1))
    for i in range(len(w_k)):
        for j in range(len(w_k)):
            J,B,N=JBN(x_k[i],y_k[j],x1,y1,x2,y2,x3,y3,x4,y4)
            F=F+q*np.linalg.det(J)*np.transpose(N)
    return F
    
##Pour nr éléments dans l'épaisseur et na éléments circonférentiel à 4 noeuds (2D)

#Table de coordonnées

r=[Rint+i*(Rext-Rint)/(nr) for i in range(nr+1)]
theta=[i*(2*ma.pi)/(na) for i in range(na+1)]

R=[]
for i in range(na+1):
    R.extend(r)

THETA=[]
for i in range(nr+1):
    THETA.extend(theta)

def tri_insertion(L): #Algorithme de tri   
    for i,v in enumerate(L):
        j=i
        while 0<j and v<L[j-1]:
            L[j]=L[j-1]
            j=j-1
        L[j]=v
    return L

def sign(n): #Fonction signe
    if n<0:
        return -1
    else:
        return 1

Theta=tri_insertion(THETA)
coord = np.zeros(((nr+1)*(na+1),2))
for i in range(len(R)):
    coord[i,0]=R[i]*ma.cos(-Theta[i])
    coord[i,1]=R[i]*ma.sin(-Theta[i])

for i in range(nr+1):
    coord=np.delete(coord,-1,axis=0)

#Table de connectivité

connect = np.zeros((nr*na,4))
L=[j for j in range(1,nr*na)]
for k in range(len(L)+1):
    connect[k,0]=ma.floor(k/nr)
for i in range (nr*na):
    connect[i,0]=connect[i,0]+1+i
    connect[i,1]=connect[i,0]+1
    connect[i,2]=connect[i,1]+(nr+1)
    connect[i,3]=connect[i,2]-1

for i in range(len(connect)):
    for j in range(4):
        if connect[i,j]-(nr+1)*na<0.5:
            connect[i,j]=connect[i,j]
        else:
            connect[i,j]=connect[i,j]-(nr+1)*na

#Table de localisation

location=np.zeros((nr*na,8))
for i in range(na*nr):
    for j in range(4):
        location[i,2*j]=2*connect[i,j]-1
        location[i,2*j+1]=2*connect[i,j]

#Noeuds sur le rayon extérieur

N_ext=[(nr+1)*(i+1) for i in range(na)]

#Noeuds sur le rayon intérieur

N_int=[(nr+1)*(i+1)-nr for i in range(na)]

#Assemblage

def MAT_CAPA_THER_ASSE_2D(): #Assemblage de la matrice de capacité thermique
    C=np.zeros(((nr+1)*(na),(nr+1)*(na)))
    for i in range(len(connect)):
        M=MAT_CAPA_THER(int(coord[int(connect[i,0])-1,0]),int(coord[int(connect[i,0])-1,1]),int(coord[int(connect[i,1])-1,0]),int(coord[int(connect[i,1])-1,1]),int(coord[int(connect[i,2])-1,0]),int(coord[int(connect[i,2])-1,1]),int(coord[int(connect[i,3])-1,0]),int(coord[int(connect[i,3])-1,1]))
        C[int(connect[i,0])-1,int(connect[i,0])-1]+=M[0,0]
        C[int(connect[i,0])-1,int(connect[i,1])-1]+=M[0,1]
        C[int(connect[i,0])-1,int(connect[i,2])-1]+=M[0,2]
        C[int(connect[i,0])-1,int(connect[i,3])-1]+=M[0,3]
        C[int(connect[i,1])-1,int(connect[i,0])-1]+=M[1,0]
        C[int(connect[i,1])-1,int(connect[i,1])-1]+=M[1,1]
        C[int(connect[i,1])-1,int(connect[i,2])-1]+=M[1,2]
        C[int(connect[i,1])-1,int(connect[i,3])-1]+=M[1,3]
        C[int(connect[i,2])-1,int(connect[i,0])-1]+=M[2,0]
        C[int(connect[i,2])-1,int(connect[i,1])-1]+=M[2,1]
        C[int(connect[i,2])-1,int(connect[i,2])-1]+=M[2,2]
        C[int(connect[i,2])-1,int(connect[i,3])-1]+=M[2,3]
        C[int(connect[i,3])-1,int(connect[i,0])-1]+=M[3,0]
        C[int(connect[i,3])-1,int(connect[i,1])-1]+=M[3,1]
        C[int(connect[i,3])-1,int(connect[i,2])-1]+=M[3,2]
        C[int(connect[i,3])-1,int(connect[i,3])-1]+=M[3,3]
    return(C)
    
def MAT_RIGI_THER_ASSE_2D(): #Assemblage de la matrice de conductivité thermique
    C=np.zeros(((nr+1)*(na),(nr+1)*(na)))
    for i in range(len(connect)):
        M=MAT_RIGI_THER(int(coord[int(connect[i,0])-1,0]),int(coord[int(connect[i,0])-1,1]),int(coord[int(connect[i,1])-1,0]),int(coord[int(connect[i,1])-1,1]),int(coord[int(connect[i,2])-1,0]),int(coord[int(connect[i,2])-1,1]),int(coord[int(connect[i,3])-1,0]),int(coord[int(connect[i,3])-1,1]))
        C[int(connect[i,0])-1,int(connect[i,0])-1]+=M[0,0]
        C[int(connect[i,0])-1,int(connect[i,1])-1]+=M[0,1]
        C[int(connect[i,0])-1,int(connect[i,2])-1]+=M[0,2]
        C[int(connect[i,0])-1,int(connect[i,3])-1]+=M[0,3]
        C[int(connect[i,1])-1,int(connect[i,0])-1]+=M[1,0]
        C[int(connect[i,1])-1,int(connect[i,1])-1]+=M[1,1]
        C[int(connect[i,1])-1,int(connect[i,2])-1]+=M[1,2]
        C[int(connect[i,1])-1,int(connect[i,3])-1]+=M[1,3]
        C[int(connect[i,2])-1,int(connect[i,0])-1]+=M[2,0]
        C[int(connect[i,2])-1,int(connect[i,1])-1]+=M[2,1]
        C[int(connect[i,2])-1,int(connect[i,2])-1]+=M[2,2]
        C[int(connect[i,2])-1,int(connect[i,3])-1]+=M[2,3]
        C[int(connect[i,3])-1,int(connect[i,0])-1]+=M[3,0]
        C[int(connect[i,3])-1,int(connect[i,1])-1]+=M[3,1]
        C[int(connect[i,3])-1,int(connect[i,2])-1]+=M[3,2]
        C[int(connect[i,3])-1,int(connect[i,3])-1]+=M[3,3]
    return(C)

def VEC_FLUX_NOD_ASSE_2D(): #Assemblage du vecteur des flux nodaux
    F=np.zeros(((nr+1)*(na),1))
    for i in range(len(connect)):
        V=VEC_FLUX_NOD(int(coord[int(connect[i,0])-1,0]),int(coord[int(connect[i,0])-1,1]),int(coord[int(connect[i,1])-1,0]),int(coord[int(connect[i,1])-1,1]),int(coord[int(connect[i,2])-1,0]),int(coord[int(connect[i,2])-1,1]),int(coord[int(connect[i,3])-1,0]),int(coord[int(connect[i,3])-1,1]))
        F[int(connect[i,0])-1,0]+=V[0,0]
        F[int(connect[i,1])-1,0]+=V[1,0]
        F[int(connect[i,2])-1,0]+=V[2,0]
        F[int(connect[i,3])-1,0]+=V[3,0]
    return(F)

def VEC_FORCE_NOD_ASSE_2D(): #Assemblage du vecteur des forces nodales
    F=np.zeros(((nr+1)*(2*na),1))
    for k in range(len(N_int)):
        F[2*(nr+1)*k,0]=F_l*ma.cos(Theta[N_int[k]-1])
        F[2*(nr+1)*k+1,0]=F_l*ma.sin(Theta[k*(nr+1)])
    return F

##Calcul du champ de température 

#Champ de température à l'intant initial

# Create an new Excel file and add a worksheet

workbook = xlsxwriter.Workbook('T_air.xlsx')
worksheet = workbook.add_worksheet('CHAMP_TEMP')

T0 = np.ones(((nr+1)*na, 1))
T0=To*T0

def CALC_CHAMP_TEMP(T0,C,K,F):
    TEMP=T0
    T_t=[To+i*(Tf-To)/N_time for i in range(N_time+1)]
    for j in N_ext:
            TEMP[j-1,0]=Tf
    Mat_inv=np.linalg.inv((1/delta_t)*C+K)
    for i in range(N_time):
        T=np.dot(Mat_inv,F+np.dot((1/delta_t)*C,np.reshape(TEMP[:,i],(na*(nr+1),1))))
        for k in N_ext:
            T_change=T
            T_change[k-1,:]=Tf
        TEMP=np.hstack((TEMP,T_change))
    return TEMP
    

CHAMP_TEMP=CALC_CHAMP_TEMP(T0,MAT_CAPA_THER_ASSE_2D(),MAT_RIGI_THER_ASSE_2D(),VEC_FLUX_NOD_ASSE_2D())
l,c=CHAMP_TEMP.shape
for j in range(c):
    for i in range(l):
        worksheet.write(i, j, CHAMP_TEMP[i,j])

#Ecriture du classeur sur le disque

workbook.close()

##Tracé des différents graphes

#Visualisation du maillage

plt.figure()
a=plt.subplot(121)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
a.plot(x1,x2,color="r")
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
a.plot(x1,x2,label="Air scheme",color="r")
plt.scatter(coord[:,0],coord[:,1],label="Mesh nodes",s=5)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("Mesh ("+str((nr+1)*na)+" Nodes and "+str(na*nr)+" Elements)",fontsize=16)
plt.xlabel("$x (mm)$",style='italic',fontsize=16)
plt.ylabel("$y (mm)$",style='italic',fontsize=16)
plt.legend(loc='upper left')
plt.tick_params(labelsize=16)

b=plt.subplot(122)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
b.plot(x1,x2,color="r")
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
b.plot(x1,x2,label="Air scheme",color="r")
plt.scatter(coord[:,0],coord[:,1],label="Mesh nodes",s=15)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("Mesh detail",fontsize=16)
plt.xlabel("$x (mm)$",style='italic',fontsize=16)
plt.ylabel("$y (mm)$",style='italic',fontsize=16)
plt.axis([-Rext,-Rint,-Rext/6,Rext/6])
plt.legend(loc='upper left')
plt.tick_params(labelsize=16)

plt.show()

#Tracé de la température en fonction du temps

plt.figure()
Temperature=CHAMP_TEMP[0,:]
Temperature1=CHAMP_TEMP[ma.floor(nr/2),:]
Temperature2=CHAMP_TEMP[nr,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,Temperature,Temps,Temperature1,Temps,Temperature2)
plt.title("Temperature variation as a function of time",fontsize=16)
plt.xlabel("$t (s)$",style='italic',fontsize=16)
plt.ylabel("T ($^\circ$K)",style='italic',fontsize=16)
plt.legend(["Temperature evolution at the inner radius","Temperature evolution at the mean radius","Temperature evolution at the outer radius"])
plt.tick_params(labelsize=12)
plt.show()

#Cartographie  du champ de température

plt.figure()
plt.suptitle('Temperature field at different time steps ($^\circ$K)',fontsize=16)
a=plt.subplot(231)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
ax=plt.scatter(coord[:,0],coord[:,1],c=CHAMP_TEMP[:,0],cmap='jet',s=50,label="Temperature in the nodes")
plt.colorbar(ax)
#ax.set_clim(min(CHAMP_TEMP[:,0]),max(CHAMP_TEMP[:,0]))
ax.set_clim(To,Tf)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("t="+str(round(0,2))+"s")
plt.ylabel("$y (mm)$",style='italic',fontsize=16)
plt.tick_params(labelsize=12)

b=plt.subplot(236)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
ax=plt.scatter(coord[:,0],coord[:,1],c=CHAMP_TEMP[:,-1],cmap='jet',s=50,label="Temperature in the nodes")
plt.colorbar(ax)
#ax.set_clim(min(CHAMP_TEMP[:,-1]),max(CHAMP_TEMP[:,-1]))
ax.set_clim(To,Tf)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("t="+str(round(time,2))+"s") 
plt.xlabel("$x (mm)$",style='italic',fontsize=16)
plt.tick_params(labelsize=12)

c=plt.subplot(232)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
ax=plt.scatter(coord[:,0],coord[:,1],c=CHAMP_TEMP[:,ma.floor(0.05*N_time)],cmap='jet',s=50,label="Temperature in the nodes")
plt.colorbar(ax)
#ax.set_clim(min(CHAMP_TEMP[:,ma.floor(N_time/3)]),max(CHAMP_TEMP[:,ma.floor(N_time/3)]))
ax.set_clim(To,Tf)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("t="+str(round(0.05*time,2))+"s")
plt.tick_params(labelsize=12)

c=plt.subplot(233)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
ax=plt.scatter(coord[:,0],coord[:,1],c=CHAMP_TEMP[:,ma.floor(0.1*N_time)],cmap='jet',s=50,label="Température aux noeuds")
plt.colorbar(ax)
#ax.set_clim(min(CHAMP_TEMP[:,ma.floor(2*N_time/3)]),max(CHAMP_TEMP[:,ma.floor(2*N_time/3)]))
ax.set_clim(To,Tf)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("t="+str(round(0.1*time,2))+"s")
plt.tick_params(labelsize=12)

c=plt.subplot(234)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
ax=plt.scatter(coord[:,0],coord[:,1],c=CHAMP_TEMP[:,ma.floor(0.15*N_time)],cmap='jet',s=50,label="Température aux noeuds")
plt.colorbar(ax)
#ax.set_clim(min(CHAMP_TEMP[:,ma.floor(2*N_time/3)]),max(CHAMP_TEMP[:,ma.floor(2*N_time/3)]))
ax.set_clim(To,Tf)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("t="+str(round(0.15*time,2))+"s")
plt.xlabel("$x (mm)$",style='italic',fontsize=16)
plt.ylabel("$y (mm)$",style='italic',fontsize=16)
plt.tick_params(labelsize=12)

c=plt.subplot(235)
theta0 = np.linspace(0, 2*np.pi, na)
rint = (Rint)
x1 = rint*np.cos(theta0)
x2 = rint*np.sin(theta0)
rext = (Rext)
x1 = rext*np.cos(theta0)
x2 = rext*np.sin(theta0)
ax=plt.scatter(coord[:,0],coord[:,1],c=CHAMP_TEMP[:,ma.floor(0.3*N_time)],cmap='jet',s=50,label="Température aux noeuds")
plt.colorbar(ax)
#ax.set_clim(min(CHAMP_TEMP[:,ma.floor(2*N_time/3)]),max(CHAMP_TEMP[:,ma.floor(2*N_time/3)]))
ax.set_clim(To,Tf)
plt.axis('equal')
plt.axis([-1.2*Rext, 1.2*Rext, -1.2*Rext, 1.2*Rext])
plt.title("t="+str(round(0.3*time,2))+"s")
plt.xlabel("$x (mm)$",style='italic',fontsize=16)
plt.tick_params(labelsize=12)

plt.show()

#Evolution de la température en espace et en temps

X=[]
X1=coord[:,0].tolist()
for i in range(time+1):
    X.extend(X1)
    
Y=[]
Y1=coord[:,1].tolist()
for i in range(time+1):
    Y.extend(Y1)

t=[]
z=[i for i in range(time+1)]
for i in range((nr+1)*na):
    t.extend(z)
T=tri_insertion(t)

ColorTemp=[]
for i in range(time+1):
    ColorTemp.extend(CHAMP_TEMP[:,int(i/delta_t)].tolist())
    
fig=go.Figure(data=[go.Scatter3d(x=X,y=Y,z=T,mode='markers',marker=dict(size=5,color=ColorTemp,colorbar=dict(
            title="Temperature (°K)"
        ),colorscale='Jet',opacity=0.8))])

fig.update_layout(scene = dict(
        xaxis = dict(title='x (mm)',nticks=4, range=[-1.5*Rext,1.5*Rext],),
                     yaxis = dict(title='y (mm)',nticks=4, range=[0,3*Rext],),
                     zaxis = dict(title='t (s)',nticks=4, range=[0,time],),),margin=dict(l=0,r=0,b=0,t=0))
fig.update_xaxes(title_text='x (mm)')
fig.update_yaxes(title_text='y (mm)')
fig.update_xaxes(title_text='t (s)')

py.offline.plot(fig,filename='Champ_de_temperature_Air.html')

#Ecriture des données dans un fichier Excel

path = r"D:\EPF 5A\Projet Recherche\Champ_Temperature.xls"
 
# On créer un "classeur"
classeur = Workbook()
# On ajoute une feuille au classeur
feuille = classeur.add_sheet("CHAMP_TEMP")
for i in range(len(X1)):
    feuille.write(i, 0, X1[i])
for i in range(len(Y1)):
    feuille.write(i, 1, Y1[i])
for j in z:
    for i in range(len(X1)):
        feuille.write(i, j+2, CHAMP_TEMP[i,int(j/delta_t)])
# Ecriture du classeur sur le disque
classeur.save(path)