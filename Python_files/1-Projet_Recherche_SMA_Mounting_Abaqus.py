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

#Paramètres géométriques

Rint=22 #Rayon intérieur
Rext=30 #Rayon extérieur 

#Paramètres temporels et discrétisation 

time=50 #durée
N_time=500 #STEP de calcul
delta_t=time/N_time #pas de temps
nr=6 #Nombre d'éléments dans l'épaisseur en 2D
na=120 #Nombre d'éléments circonférentiel en 2D

#Paramètres matériaux et physiques

rho=6.45/1000000 #Masse volumique
E=28e3 #Module d'Young nitinol ambiant
nu=0.3 #Coefficient de Poisson
To=20+273 #Température homogène à l'intant initial
Tf=70+273 #Température finale de chauffe 
P=2.3 #Pression
F_l=P*2*ma.pi*Rint/na #Force linéique appliqué à chaque noeud

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

##Calcul de la contrainte mécanique

# ouverture du fichier Excel 

wb = xlrd.open_workbook('SIGMA_MECA_Abaqus.xlsx')
 
# lecture des données dans chaque feuille

RAYON= wb.sheet_by_name(u'xyToExcel')

nb_ligne=RAYON.nrows
nb_colonne=RAYON.ncols

Rayon=np.zeros((nb_ligne,nb_colonne))
for i in range(nb_colonne):
    Rayon[:,i]=RAYON.col_values(i)
#Réduction de modèle (ligne d'intérêt theta=0)

Reduc=np.zeros((7,3))
for i in range(7):
    Reduc[i,0]=abs(Rayon[i,0])
    Reduc[i,1]=Rayon[i,2]+R[i]
    Reduc[i,2]=R[i]

# Create an new Excel file and add a worksheet

workbook = xlsxwriter.Workbook('SIGMA_MECA.xlsx')
worksheet1 = workbook.add_worksheet('SIGMA_MECA')

l1,c1=Reduc.shape
for i in range(l1):
    for j in range(c1):
        worksheet1.write(i, j, Reduc[i,j])
        
#Ecriture du classeur sur le disque

workbook.close()


