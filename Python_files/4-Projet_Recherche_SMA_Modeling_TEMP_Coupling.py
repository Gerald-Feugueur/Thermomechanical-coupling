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
import plotly.tools as tls

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
q=0 #Flux de chaleur
E=28e3 #Module d'Young
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

l=(E*nu)/((1+nu)*(1-2*nu)) #Coefficient de lamé
m=(E)/(2*(1+nu)) #Coefficient de lamé
Gamma=l+2*m #Module d'Young modifié

#Calcul du couplage faible et fort imbriqué

CONTR=np.zeros((7,N_time+1))
for i in range(7):
    CONTR[i,0]=CONTRMECA[i,0]
    
def CALC_CPLG(DEFO_THER,CONTR):
    CONTR_weak=CONTR
    CONTR_strong=np.zeros((7,N_time+1))
    T=np.zeros((7,N_time+1))
    for i in range(N_time):
        CONTR_weak[:,i+1]=CONTR_weak[:,i]+DEFO_THER[:,i+1]
    for j in range(N_time+1):
        T[:,j]=((1/(alpha)))*CONTR_weak[:,j]
        CONTR_strong[:,j]=CONTR_weak[:,j]+alpha*T[:,j]
    return Gamma*CONTR_weak,Gamma*CONTR_strong,T

CONTR_weak,CONTR_strong,T=CALC_CPLG(DEFO_THER,CONTR)

##Tracé des différents graphes

#Evolution du rayon exterieur en fonction du temps et selon le chargement

plt.figure()
Rayon_ther=RAYON_THER[6,:]
plt.title('Outer radius variation by temperature change',fontsize=16)
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,Rayon_ther,label='Outer radius variation by temperature change')
#plt.title("Evolution du rayon extérieur en fonction du temps")
plt.xlabel("$t (s)$",style='italic',fontsize=16)
plt.ylabel("$R_{ext} (mm)$",style='italic',fontsize=16)
plt.legend(["Outer radius variation by temperature change"])
plt.tick_params(labelsize=16)
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

plt.suptitle("Evolution of outer node displacement as a function of time",fontsize=16)
b=plt.subplot(311)
Deplacement=DEPL_THER[6,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,1000*Deplacement)
plt.ylabel("$u_{R_{ext}} (\mu m)$",style='italic',fontsize=16)
plt.legend(["Evolution of thermal node displacement as a function of time"])
plt.tick_params(labelsize=16)

c=plt.subplot(312)
Temperature=CHAMP_TEMP[6,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,Temperature)
plt.ylabel("T ($^\circ$K)",style='italic',fontsize=16)
plt.legend(["Temperature evolution in the outer radius"])
plt.tick_params(labelsize=16)

c=plt.subplot(313)
Temperature=CHAMP_TEMP[6,:]
Temps=[delta_t*i for i in range(N_time+1)]

Derive=[]
for i in range(len(Temps)-1):
    Derive.append((Temperature[i+1]-Temperature[i])/(Temps[i+1]-Temps[i]))

plt.plot(np.delete(Temps,-1,axis=0),Derive)
plt.xlabel("$t (s)$",style='italic',fontsize=16)
plt.ylabel("dT/dt ($^\circ$K/s)",style='italic',fontsize=16)
plt.legend(["Temperature evolution in the outer radius"])
plt.tick_params(labelsize=16)

plt.show()

#Evolution de la contrainte 2D

fig = plotly.tools.make_subplots(rows=2, cols=1,subplot_titles=("Stress as a function of radius and time for a weak coupling", "Stress as a function of radius and time for a strong coupling"))

trace=go.Contour(z=CONTR_weak,line_smoothing=1,colorscale="Jet",colorbar=dict(len=0.5, y=0.80,),zmin=np.min(np.min(CONTR_weak)),zmax=np.max(np.max(CONTR_strong)))
trace1=go.Contour(z=CONTR_strong,line_smoothing=1,colorscale="Jet",colorbar=dict(len=0.5, y=0.20),zmin=np.min(np.min(CONTR_strong)),zmax=np.max(np.max(CONTR_strong)))

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

py.offline.plot(fig,filename='Champ_de_Contrainte_2D.html')

#Evolution de la contrainte 3D

fig0=go.Figure(data=[go.Surface(z=CONTR_weak,colorscale="Jet"),go.Surface(z=CONTR_strong,colorscale="Jet")])

fig0.update_layout(scene = dict(
        xaxis = dict(title='t (s)',nticks=4,),
                     yaxis = dict(title='r (mm)',nticks=4,),
                     zaxis = dict(title='sigma (MPa)',nticks=4,),),margin=dict(l=0,r=0,b=0,t=0))
fig0.update_xaxes(title_text='t (s)')
fig0.update_yaxes(title_text='y (mm)')
fig0.update_xaxes(title_text='sigma (MPa)')

py.offline.plot(fig0,filename='Champ_de_Contrainte_3D.html')

#Evolution de la déformation 2D

fig2 = plotly.tools.make_subplots(rows=2, cols=1,subplot_titles=("Deformation as a function of radius and time for a weak coupling", "Deformation as a function of radius and time for a strong coupling"))

trace2=go.Contour(z=(1/Gamma)*CONTR_weak,line_smoothing=1,colorscale="Jet",colorbar=dict(len=0.5, y=0.80,),zmin=np.min(np.min((1/Gamma)*CONTR_weak)),zmax=np.max(np.max((1/Gamma)*CONTR_strong)))
trace12=go.Contour(z=(1/Gamma)*CONTR_strong,line_smoothing=1,colorscale="Jet",colorbar=dict(len=0.5, y=0.20),zmin=np.min(np.min((1/Gamma)*CONTR_strong)),zmax=np.max(np.max((1/Gamma)*CONTR_strong)))

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

py.offline.plot(fig2,filename='Champ_de_Deformation_2D.html')

#Evolution de la déformation 3D

fig02=go.Figure(data=[go.Surface(z=(1/Gamma)*CONTR_weak,colorscale="Jet"),go.Surface(z=(1/Gamma)*CONTR_strong,colorscale="Jet")])

fig02.update_layout(scene = dict(
        xaxis = dict(title='t (s)',nticks=4,),
                     yaxis = dict(title='r (mm)',nticks=4,),
                     zaxis = dict(title='sigma (MPa)',nticks=4,),),margin=dict(l=0,r=0,b=0,t=0))
fig02.update_xaxes(title_text='t (s)')
fig02.update_yaxes(title_text='y (mm)')
fig02.update_xaxes(title_text='epsilon (mm/mm)')

py.offline.plot(fig02,filename='Champ_de_Deformation_3D.html')

#Tracé de la deformation et contrainte en fonction du temps

plt.figure()

plt.suptitle("Strain for a weak coupling (left) and a strong coupling (right)",fontsize=16)

plt.subplot(121)
defo=(1/Gamma)*CONTR_weak[0,:]
defo1=(1/Gamma)*CONTR_weak[ma.floor(nr/2),:]
defo2=(1/Gamma)*CONTR_weak[nr,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,100*defo,Temps,100*defo1,Temps,100*defo2)
plt.title("Strain evolution as a function of time",fontsize=16)
plt.ylabel("$\epsilon_{rr} $"+"(%)",style='italic',fontsize=16)
plt.xlabel("$t (s)$",style='italic',fontsize=16)
plt.legend(["Strain evolution at the inner radius","Strain evolution at the mean radius","Strain evolution at the outer radius"])
plt.tick_params(labelsize=16)

plt.subplot(122)
defo0=(1/Gamma)*CONTR_strong[0,:]
defo10=(1/Gamma)*CONTR_strong[ma.floor(nr/2),:]
defo20=(1/Gamma)*CONTR_strong[nr,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,100*defo0,Temps,100*defo10,Temps,100*defo20)
plt.title("Strain evolution as a function of time",fontsize=16)
plt.xlabel("$t (s)$",style='italic',fontsize=16)
plt.ylabel("$\epsilon_{rr} $"+"(%)",style='italic',fontsize=16)
plt.legend(["Strain evolution at the inner radius","Strain evolution at the mean radius","Strain evolution at the outer radius"])
plt.tick_params(labelsize=16)

plt.show()

plt.figure()

plt.suptitle("Stress for a weak coupling (left) and a strong coupling (right)",fontsize=16)

plt.subplot(121)
defoa=CONTR_weak[0,:]
defo1a=CONTR_weak[ma.floor(nr/2),:]
defo2a=CONTR_weak[nr,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,defoa,Temps,defo1a,Temps,defo2a)
plt.title("Stress evolution as a function of time",fontsize=16)
plt.xlabel("$t (s)$",style='italic',fontsize=16)
plt.ylabel("$\sigma_{rr} (MPa)$",style='italic',fontsize=16)
plt.legend(["Stress evolution at the inner radius","Stress evolution at the mean radius","Stress evolution at the outer radius"])
plt.tick_params(labelsize=16)

plt.subplot(122)
defo0b=CONTR_strong[0,:]
defo10b=CONTR_strong[ma.floor(nr/2),:]
defo20b=CONTR_strong[nr,:]
Temps=[delta_t*i for i in range(N_time+1)]
plt.plot(Temps,defo0b,Temps,defo10b,Temps,defo20b)
plt.title("Stress evolution as a function of time",fontsize=16)
plt.xlabel("$t (s)$",style='italic',fontsize=16)
plt.ylabel("$\sigma_{rr} (MPa)$",style='italic',fontsize=16)
plt.legend(["Stress evolution at the inner radius","Stress evolution at the mean radius","Stress evolution at the outer radius"])
plt.tick_params(labelsize=16)

plt.show()