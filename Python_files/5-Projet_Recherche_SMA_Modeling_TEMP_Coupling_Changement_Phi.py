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
import time
from mpl_toolkits.mplot3d import Axes3D
import xlsxwriter, pandas
from scipy.integrate import quad
from scipy.optimize import bisect

# ouverture du fichier Excel 

wb = xlrd.open_workbook('SIGMA_MECA.xlsx')
 
# lecture des données dans chaque feuille

RAYON= wb.sheet_by_name(u'SIGMA_MECA')

nb_ligne=RAYON.nrows
nb_colonne=RAYON.ncols

Rayon=np.zeros((nb_ligne,nb_colonne))
for i in range(nb_colonne):
    Rayon[:,i]=RAYON.col_values(i)


CONTRMECA=np.zeros((7,1))
for i in range(7):
    CONTRMECA[i,0]=Rayon[i,0]

#Paramètres géométriques

Rint=Rayon[0,1] #Rayon intérieur
Rext=Rayon[-1,1] #Rayon extérieur 

#Paramètres temporels et discrétisation 

time=50 #durée
N_time=500 #STEP de calcul
delta_t=time/N_time #pas de temps
nr=6 #Nombre d'éléments dans l'épaisseur en 2D
na=120 #Nombre d'éléments circonférentiel en 2D

#Paramètres matériaux et physiques

rho=6.45/1000000 #Masse volumique
lamda=18/1000 #Conductivité thermique
cp=322 #Capacité thermique
alpha=11/10**6 #Dilatation thermique 
alpha_A=11/10**6 #Dilatation thermique Austenite
alpha_M=6.6/10**6 #Dilatation thermique Martensite
q=0 #Flux de chaleur
E_M=28e3 #Module d'Young Martensite
E_A=75e3 #Module d'Young Austenite
nu=0.3 #Coefficient de Poisson
To=20+273 #Température homogène à l'intant initial
Tf=70+273 #Température finale de chauffe 
P=0.85 #Pression
F_l=P*2*ma.pi*Rint/na #Force linéique appliqué aux noeuds

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

##Fonctions pour calculer la contrainte vraie 

#Calcul du champ de température 

# ouverture du fichier Excel 

wb = xlrd.open_workbook('T_bague.xlsx')
 
# lecture des données dans chaque feuille

TEMP_BAGUE= wb.sheet_by_name(u'CHAMP_TEMP')

nb_ligne=TEMP_BAGUE.nrows
nb_colonne=TEMP_BAGUE.ncols

T_BAGUE=np.zeros((nb_ligne,nb_colonne))
for i in range(nb_colonne):
    T_BAGUE[:,i]=TEMP_BAGUE.col_values(i)
   
CHAMP_TEMP=T_BAGUE[:nr+1,:]

#Calcul de la déformation thermique

def CALC_DEFO_THER(CHAMP_TEMP):
    Epsilon=np.zeros((nr+1,1))
    for i in range(N_time):
        DIFF_TEMP=np.reshape(CHAMP_TEMP[:,i+1]-CHAMP_TEMP[:,i],((nr+1),1))
        Epsilon=np.hstack((Epsilon,DIFF_TEMP))
    return alpha*Epsilon

DEFO_THER=CALC_DEFO_THER(CHAMP_TEMP)

#Calcul du déplacement thermique

def CALC_DEPL_THER(DEFO_THER):
    DEPL=np.zeros((((nr+1)),N_time+1))
    r=np.zeros((((nr+1)),N_time+1))
    for j in range((nr+1)):
        r[j,0]=ma.sqrt(coord[j,0]**2+coord[j,1]**2)
        for i in range(N_time):
            r[0,i+1]=r[0,0]
    for i in range(N_time):
        for k in range((nr)):
            DEPL[k+1,i]=DEPL[k,i]+DEFO_THER[k+1,i]*(r[k+1,i]-r[k,i])
            r[k+1,i+1]=r[k+1,i]+DEPL[k+1,i]
    return DEPL,r

DEPL_THER,RAYON_THER=CALC_DEPL_THER(DEFO_THER)

#Coefficients de Lamé

# l=(E*nu)/((1+nu)*(1-2*nu)) #Coefficient de lamé
# m=(E)/(2*(1+nu)) #Coefficient de lamé

#Fonction de densité de probabilité

def DDP(x,sigma,mean):
    return ma.exp(-(x-mean)**2/(2*sigma**2))/(sigma*ma.sqrt(2*ma.pi))

#Fonctions de répartition

meanE=45+273
sigmaE=0.605

def phiE(x):
    return E_M+(E_A-E_M)*(1.0 + ma.erf((x-meanE)/(sigmaE*ma.sqrt(2.0))))/2.0
    
meanP=45+273
sigmaP=0.605

def phiP(x):
    return 0.001*(1.0 + ma.erf((x-meanP)/(sigmaP*ma.sqrt(2.0))))/2.0
    
meanP=45+273
sigmaP=0.605

def phi(x):
    return (1.0 + ma.erf((x-meanP)/(sigmaP*ma.sqrt(2.0))))/2.0
    
meanE=45+273
sigmaE=0.605

def phi_alpha(x):
    return alpha_M+(alpha_A-alpha_M)*(1.0 + ma.erf((x-meanE)/(sigmaE*ma.sqrt(2.0))))/2.0

#Calcul du couplage faible et fort imbriqué

CONTR=np.zeros((7,N_time+1))
for i in range(7):
    CONTR[i,0]=CONTRMECA[i,0]

DEFO_PHASE=np.zeros((7,N_time+1))
for i in range(nr+1):
    for j in range(N_time+1):
        DEFO_PHASE[i,j]=phiP(CHAMP_TEMP[i,j])
        
ALPHA=np.zeros((7,N_time+1))
for i in range(nr+1):
    for j in range(N_time+1):
        ALPHA[i,j]=phi_alpha(CHAMP_TEMP[i,j])

def CALC_CPLG(DEFO_THER,CONTR):
    CONTR_weak=CONTR
    CONTR_strong=np.zeros((7,N_time+1))
    T=np.zeros((7,N_time+1))
    #Gamma=l+2*m
    for k in range(7):
        for i in range(N_time):
            CONTR_weak[:,i+1]=CONTR_weak[:,i]+(1/alpha)*DEFO_THER[:,i+1]*ALPHA[k,i+1]
    CONTR_weak=CONTR_weak+DEFO_PHASE*(1)
    for k in range(7):
        for j in range(N_time+1):
            T[:,j]=((1/(ALPHA[k,j])))*CONTR_weak[:,j]
            CONTR_strong[:,j]=CONTR_weak[:,j]+alpha*T[:,j]
    return CONTR_weak,CONTR_strong,T

DEFO_weak,DEFO_strong,T=CALC_CPLG(DEFO_THER,CONTR)

E_PHASE=np.zeros((7,N_time+1))
for i in range(nr+1):
    for j in range(N_time+1):
        E_PHASE[i,j]=phiE(CHAMP_TEMP[i,j])

CONTR_weak=np.zeros((7,N_time+1))
for i in range(nr+1):
    for j in range(N_time+1):
        CONTR_weak[i,j]=E_PHASE[i,j]*DEFO_weak[i,j]*((1+(nu/(1-2*nu)))/((1+nu)))
        
CONTR_strong=np.zeros((7,N_time+1))
for i in range(nr+1):
    for j in range(N_time+1):
        CONTR_strong[i,j]=E_PHASE[i,j]*DEFO_strong[i,j]*((1+(nu/(1-2*nu)))/((1+nu)))

##Tracé des différents graphes

#Evolution du rayon intérieur en fonction du temps et selon le chargement

plt.figure()
Rayon_ther=RAYON_THER[6,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,Rayon_ther,label='Outer radius evolution with a thermal charge')
plt.title("Outer radius evolution as a function of time")
plt.xlabel("$t (s)$")
plt.ylabel("$R_{ext} (mm)$",style='italic')
#plt.legend(["Evolution du rayon extérieur avec un chargement thermique"])
plt.show()

workbook = xlsxwriter.Workbook('T_ext.xlsx')
worksheet = workbook.add_worksheet('CHAMP_TEMP')

l=len(Rayon_ther)
for i in range(l):
        worksheet.write(i, 0, Rayon_ther[i])

#Ecriture du classeur sur le disque

workbook.close()

#Evolution du déplacement nodal en fonction du temps

plt.figure()

plt.suptitle("Evolution of nodal displacement of outer radius as a function of time")
b=plt.subplot(311)
Deplacement=DEPL_THER[6,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,Deplacement)
plt.ylabel("$u_{R_{ext}} (mm)$",style='italic')
plt.legend(["Evolution of thermal nodal displacement as a function of time"])

c=plt.subplot(312)
Temperature=CHAMP_TEMP[6,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,Temperature)
plt.xlabel("$t (s)$")
plt.ylabel("T ($^\circ$K)",style='italic')
plt.legend(["Temperature evolution in the outer radius"])

c=plt.subplot(313)
Temperature=CHAMP_TEMP[6,:]
Temps=[delta_t*i for i in range(N_time+1)]

Derive=[]
for i in range(len(Temps)-1):
    Derive.append((Temperature[i+1]-Temperature[i])/(Temps[i+1]-Temps[i]))

plt.plot(np.delete(Temps,-1,axis=0),Derive)
plt.xlabel("$t (s)$")
plt.ylabel("dT/dt ($^\circ$K/s)",style='italic')
plt.legend(["Evolution of temperature variation in the outer radius"])

plt.show()

#Evolution de la contrainte 2D

fig = plotly.tools.make_subplots(rows=2, cols=1,subplot_titles=("Stress as a function of radius and time for a weak coupling", "Stress as a function of radius and time for a strong coupling"))

trace=go.Contour(z=CONTR_weak,colorscale="Jet",colorbar=dict(len=0.5, y=0.80),zmin=np.min(np.min(CONTR_weak)),zmax=np.max(np.max(CONTR_strong)))
trace1=go.Contour(z=CONTR_strong,colorscale="Jet",colorbar=dict(len=0.5, y=0.20),zmin=np.min(np.min(CONTR_strong)),zmax=np.max(np.max(CONTR_strong)))

fig.append_trace(trace, 1, 1)
fig.append_trace(trace1, 2, 1)

# fig.update_layout(scene = dict(
#         xaxis = dict(title='t (s)',nticks=4,),
#                      yaxis = dict(title='r (mm)',nticks=4,),
#                      zaxis = dict(title='sigma (MPa)',nticks=4,),),margin=dict(l=0,r=0,b=0,t=0))
# fig.update_xaxes(title_text='t (s)')
# fig.update_yaxes(title_text='y (mm)')
# fig.update_xaxes(title_text='sigma (MPa)')
fig.update_xaxes(title_text="Temps (ds)", row=1, col=1)
fig.update_xaxes(title_text="Temps (ds)", row=2, col=1)
fig.update_yaxes(title_text="Noeud i", row=1, col=1)
fig.update_yaxes(title_text="Noeud i", row=2, col=1)

#R(i)=R_int+1/6(R_ext-R_int)i

py.offline.plot(fig,filename='Champ_de_Contrainte_Changement_de_phase_2D.html')

#Evolution de la contrainte 3D

fig0=go.Figure(data=[go.Surface(z=CONTR_weak,colorscale="Jet"),go.Surface(z=CONTR_strong,colorscale="Jet")])

fig0.update_layout(scene = dict(
        xaxis = dict(title='t (s)',nticks=4,),
                     yaxis = dict(title='r (mm)',nticks=4,),
                     zaxis = dict(title='sigma (MPa)',nticks=4,),),margin=dict(l=0,r=0,b=0,t=0))
fig0.update_xaxes(title_text='t (s)')
fig0.update_yaxes(title_text='y (mm)')
fig0.update_xaxes(title_text='sigma (MPa)')

py.offline.plot(fig0,filename='Champ_de_Contrainte_Changement_de_phase_3D.html')

#Evolution de la déformation 2D

fig2 = plotly.tools.make_subplots(rows=2, cols=1,subplot_titles=("Deformation as a function of radius and time for a weak coupling", "Deformation as a function of radius and time for a strong coupling"))

trace2=go.Contour(z=DEFO_weak,line_smoothing=1,colorscale="Jet",colorbar=dict(len=0.5, y=0.80,),zmin=np.min(np.min(DEFO_weak)),zmax=np.max(np.max(DEFO_strong)))
trace12=go.Contour(z=DEFO_strong,line_smoothing=1,colorscale="Jet",colorbar=dict(len=0.5, y=0.20),zmin=np.min(np.min(DEFO_strong)),zmax=np.max(np.max(DEFO_strong)))

fig2.append_trace(trace2, 1, 1)
fig2.append_trace(trace12, 2, 1)

# fig.update_layout(scene = dict(
#         xaxis = dict(title='t (s)',nticks=4,),
#                      yaxis = dict(title='r (mm)',nticks=4,),
#                      zaxis = dict(title='sigma (MPa)',nticks=4,),),margin=dict(l=0,r=0,b=0,t=0))
# fig.update_xaxes(title_text='t (s)')
# fig.update_yaxes(title_text='y (mm)')
# fig.update_xaxes(title_text='sigma (MPa)')
fig2.update_xaxes(title_text="Temps (ds)", row=1, col=1)
fig2.update_xaxes(title_text="Temps (ds)", row=2, col=1)
fig2.update_yaxes(title_text="Noeud i", row=1, col=1)
fig2.update_yaxes(title_text="Noeud i", row=2, col=1)

#R(i)=R_int+1/6(R_ext-R_int)i

py.offline.plot(fig2,filename='Champ_de_Deformation_Changement_de_phase_2D.html')

#Evolution de la déformation 3D

fig02=go.Figure(data=[go.Surface(z=DEFO_weak,colorscale="Jet"),go.Surface(z=DEFO_strong,colorscale="Jet")])

fig02.update_layout(scene = dict(
        xaxis = dict(title='t (s)',nticks=4,),
                     yaxis = dict(title='r (mm)',nticks=4,),
                     zaxis = dict(title='sigma (MPa)',nticks=4,),),margin=dict(l=0,r=0,b=0,t=0))
fig02.update_xaxes(title_text='t (s)')
fig02.update_yaxes(title_text='y (mm)')
fig02.update_xaxes(title_text='epsilon (mm/mm)')

py.offline.plot(fig02,filename='Champ_de_Deformation_Changement_de_phase_3D.html')

#Fonction de répartition

plt.figure()

XT=np.linspace(20+273,70+273,N_time)
Repart=[phiE(x) for x in XT]
plt.plot(XT,Repart)
plt.title('Evolution of Young Modulus as a function of temperature')
plt.xlabel("T ($^\circ$K)",style='italic',fontsize=14)
plt.ylabel("E (MPa)",style='italic',fontsize=14)
plt.tick_params(labelsize=14)
# plt.legend()
plt.show()

plt.figure()

XT=np.linspace(20+273,70+273,N_time)
RepartP=[phiP(x) for x in XT]
plt.plot(XT,RepartP)
plt.title('Evolution of the deformation due to phase change as a function of temperature')
plt.xlabel("T ($^\circ$K)",style='italic',fontsize=14)
plt.ylabel("$\epsilon_{ \phi}$",style='italic',fontsize=16)
plt.tick_params(labelsize=14)
#plt.legend(["Evolution of the deformation due to phase change as a function of temperature"])
plt.show()

plt.figure()

XT=np.linspace(20+273,70+273,N_time)
RepartPt=[100*phi(x) for x in XT]
plt.plot(XT,RepartPt)
plt.title('Percentage of austenite as a function of temperature')
plt.xlabel("T ($^\circ$K)",style='italic',fontsize=14)
plt.ylabel("% Austenite",style='italic',fontsize=14)
plt.tick_params(labelsize=14)
#plt.legend(["Percentage of austénite as a function of temperature"])
plt.show()

plt.figure()

XT=np.linspace(20+273,70+273,N_time)
RepartPa=[phi_alpha(x) for x in XT]
plt.plot(XT,RepartPa)
plt.title('Evolution of the coefficient of thermal expansion as a function of temperature')
plt.xlabel("T ($^\circ$K)",style='italic',fontsize=14)
plt.ylabel(r'$\alpha$(T)',style='italic',fontsize=16)
plt.tick_params(labelsize=14)
#plt.legend(["Evolution de la dilatation thermique en fonction de la température"])
plt.show()

#Projection temporelle

plt.figure()
plt.title("Stress evolution as a function of radius for a strong coupling")
plt.plot(Rayon[:,1],CONTR_strong[:,0],label='t=0s')
plt.plot(Rayon[:,1],CONTR_strong[:,ma.floor(0.2*N_time)],label='t=10s')
plt.plot(Rayon[:,1],CONTR_strong[:,ma.floor(0.4*N_time)],label='t=20s')
plt.plot(Rayon[:,1],CONTR_strong[:,ma.floor(0.6*N_time)],label='t=30s')
plt.plot(Rayon[:,1],CONTR_strong[:,ma.floor(0.8*N_time)],label='t=40s')
plt.plot(Rayon[:,1],CONTR_strong[:,ma.floor(1*N_time)],label='t=50s')
plt.legend(loc='best')
plt.tick_params(labelsize=14)
plt.xlabel("r (mm)",style='italic',fontsize=14)
plt.ylabel("${\sigma_{rr}} (MPa)$",style='italic',fontsize=16)
plt.show()

#Tracé de la deformation et contrainte en fonction du temps

plt.figure()

plt.suptitle("Deformation evolution (left) and stress (right) for a strong coupling with phase change")

plt.subplot(121)
defo=DEFO_strong[0,:]
defo1=DEFO_strong[ma.floor(nr/2),:]
defo2=DEFO_strong[nr,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,defo,Temps,defo1,Temps,defo2)
plt.title("Strain evolution as a function of time")
plt.xlabel("$t (s)$",style='italic',fontsize=14)
plt.ylabel("$\epsilon_{rr} (mm/mm)$",style='italic',fontsize=16)
plt.legend(["Strain evolution at the inner radius","Strain evolution at the mean radius","Strain evolution at the outer radius"])
plt.tick_params(labelsize=14)

b=plt.subplot(122)
defo0=CONTR_strong[0,:]
defo10=CONTR_strong[ma.floor(nr/2),:]
defo20=CONTR_strong[nr,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,defo0,Temps,defo10,Temps,defo20)
plt.title("Stress evolution as a function of time")
plt.xlabel("$t (s)$",style='italic',fontsize=14)
plt.ylabel("$\sigma_{rr} (MPa)$",style='italic',fontsize=16)
plt.legend(["Stress evolution at the inner radius","Stress evolution at the mean radius","Stress evolution at the outer radius"])
plt.tick_params(labelsize=14)

plt.show()