# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:46:03 2021

@author: leo

=============================================================================
DISCLAIMER



=============================================================================
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import Audio
import matplotlib.pyplot as pyp

#Définition de nos variables
alpha = 0.1
tau = 200
mu = 6*10**(-4)
L = 0.5
sh = 0
sg = 0
c = np.sqrt(tau/mu)
s = 0
T = (2*L)/c
xmin = 0
xmax = L
nbx = int(L*1000)
nbxComparaison = 100
dt = T/100

#%%Tracer de g, h et f pour une valeur de T

print("0:t=0  1:t=T/4  2:t=T/2  3:t=3*T/4  4:t=T")
n= int(input("Quelle valeur de t ?"))

#Créer un vecteur avec nos valeurs de x
x = np.linspace(xmin, xmax, num=nbx)
#Créer un vecteur avec nos valeurs de t
Tt = np.array([0, T/4, T/2, 3*T/4, T])

#Définition de nos fonction g, h et f
def solexacte_g(t,x):
    tp = (t+T/2)%T-T/2
    sg = -(alpha/(4*c))*((L-abs(x+c*tp))**2)
    return sg

def solexacte_h(t,x):
    tp = (t+T/2)%T-T/2
    sh = (alpha/(4*c))*((L-abs(x-c*tp))**2)
    return sh

def solexacte(t,x):
    tp = (t+T/2)%T-T/2
    s = (alpha/(4*c))*((L-abs(x-c*tp))**2) - (alpha/(4*c))*((L-abs(x+c*tp))**2)
    return s

def plotT(t):
    plt.plot(x,solexacte(t,x),"k-",label="f")
    plt.plot(x,solexacte_g(t,x),"b-",label="g")
    plt.plot(x,solexacte_h(t,x),"r-",label="h")

plt.figure(1)
plotT(Tt[n])
plt.title("Tracé de h, g et f pour t=T")
plt.legend()
plt.xlabel("Spatiale (x)")
plt.ylabel("f approximation de la corde")

#%%

def animate(i): 
    t = i * dt
    y = solexacte(t,x)
    sg = solexacte_g(t,x)
    sh = solexacte_h(t,x)
    line1.set_data(x, y)
    line2.set_data(x,sg)
    line3.set_data(x,sh)
    return line1, line2, line3, 

fig = plt.figure()
line1, = plt.plot([],[],"k-",label="f")
line2, = plt.plot([],[],"b-",label="g")
line3, = plt.plot([],[],"r-",label="h")
plt.xlim (0 , L)
plt.ylim(-(alpha)/(16*c),alpha/(16*c))
plt.title("Tracé de h, g et f sur T")
plt.legend()
plt.xlabel("Spatiale (x)")
plt.ylabel("f approximation de la corde")

ani = animation.FuncAnimation(fig, animate, frames=100, blit=True, interval=20, repeat=True)

#%%

def solapprochee(t,x,p):
    Up = 0
    for i in range(1, p+1):
        fn = (i*c)/(2*L)
        kn = (i*np.pi)/L
        U = (2*alpha*L*L)/(i*i*np.pi**2*c)*np.sin(kn*x)*np.sin(2*np.pi*fn*t)
        Up = Up + U
    return Up

plt.figure(3)
plt.plot(x, solapprochee(T/3,x,1),"y-",label="Up pour p=1")
plt.plot(x, solapprochee(T/3,x,3),"b-",label="Up pour p=3")
plt.plot(x, solapprochee(T/3,x,10),"r-",label="Up pour p=10")
plt.plot(x, solapprochee(T/3,x,100),"c-",label="Up pour p=100")
plt.plot(x, solexacte(T/3,x),"k-", label="Solexacte")
plt.title("Tracé de Up pour p=1,3,10,100")
plt.legend()
plt.xlabel("Spatiale (x)")
plt.ylabel("approximation de Up(t,x)")

#%%

Majoration=np.zeros(100)

#Créer un vecteur avec nos valeurs de xComparaison
xComparaison = np.linspace(xmin, xmax, num=nbxComparaison)
Erreur = np.zeros(100)
pp = np.linspace(1, 100, 100)

for p in range(1, 100+1):
    comparaison = abs(solapprochee(T/3, xComparaison, p) - solexacte(T/3,xComparaison))
    Erreur[p-1] = np.max(comparaison)
    Majoration[p-1] = (2*alpha*L*L)/(c*np.pi**2*p)

plt.figure(4)
plt.loglog(pp, Erreur,label="Erreur")
plt.loglog(pp, Majoration,"r-",label="Majoration")
plt.title("Tracé de l'érreur et de la majoration")
plt.legend()
plt.xlabel("Spatiale (x)")
plt.ylabel("Erreur commise en fonction de p")

#%%

xBar = np.linspace(1,30,30)

def Dn(N):
    dn = (2*alpha*L)/(N*np.pi*c)*np.cos(N*np.pi)
    return dn

plt.figure(5)
plt.bar(xBar,Dn(xBar))
plt.title("Diagramme des 30 premières valeurs de dn")
plt.xlabel("numéro de la valeur")
plt.ylabel("Valeur de dn")

#%%Affiche le son sur la période T = 1/fn

def Vp(t,P):
    vp = 0
    for n in range(1,P+1):
        vp = 2*alpha*L/(n*np.pi*c)*(-1)**n*np.sin((n*np.pi*c)/L*t) + vp
    return vp

t = np.linspace(0, T, 100)

plt.figure(6)
plt.plot(t,Vp(t,30))
plt.title("Représentation du son non amorti sur T")
plt.xlabel("temps")
plt.ylabel("Valeur du son")

#%%Affiche le son sur la période 2sec

fEch = 44100 #La fréquence d'échantillonnage
t = np.linspace(0, 2, 2*fEch+1)

plt.figure(7)
plt.plot(t,Vp(t,30))
plt.title("Représentation du son non amorti")
plt.xlabel("temps")
plt.ylabel("Valeur du son")

#%%Gérérer le son non amorti

# script audio.py
# (C) Fabrice Sincère ; Jean-Claude Meilland
import wave
import math
import binascii

print("Signal non amorti")
print("Création d'un fichier audio au format WAV (PCM 8 bits stéréo 44100 Hz)")
print("Son de forme sinusoïdale\n")

NomFichier = 'son_non_amorti.wav'
Monson = wave.open(NomFichier,'w') # instanciation de l'objet Monson

nbCanal = 1    # stéreo
nbOctet = 1    # taille d'un échantillon : 1 octet = 8 bits
fech = 44100  # fréquence d'échantillonnage
duree = float(input('Durée (en secondes) ? '))

nbEchantillon = int(duree*fech)
print("Nombre d'échantillons :",nbEchantillon)

parametres = (nbCanal,nbOctet,fech,nbEchantillon,'NONE','not compressed')# tuple
Monson.setparams(parametres)    # création de l'en-tête (44 octets)

# niveau max dans l'onde positive : +1 -> 255 (0xFF)
# niveau max dans l'onde négative : -1 ->   0 (0x00)
# niveau sonore nul :                0 -> 127.5 (0x80 en valeur arrondi)

print('Veuillez patienter...')
for i in range(0,nbEchantillon):
    # canal
    # 127.5 + 0.5 pour arrondir à l'entier le plus proche
    val = wave.struct.pack('B',int(128.0 + 127.5*10000*Vp(i/fEch,30)))
    Monson.writeframes(val) # écriture frame

Monson.close()

Fichier = open(NomFichier,'rb')
data = Fichier.read()
tailleFichier = len(data)
print('\nTaille du fichier',NomFichier, ':', tailleFichier,'octets')
print("Lecture du contenu de l'en-tête (44 octets) :")
print(binascii.hexlify(data[0:44]))
print("Nombre d'octets de données :",tailleFichier - 44)
print("")
Fichier.close()

#%%

l1 = 0.25
l2 = 0.75

def envelop(Interval,N):
    env = np.ones(N)
    
    #profil de montée
    x=np.linspace(0,0.5,int(Interval[0]*N))
    for k in range(0,len(x)-1):
        env[k]=16*((1-x[k])*x[k])**2*(x[k]/(2*l1))
    
    #profil de descente
    x=np.linspace(0.5,1,int(N-Interval[1]*N))
    for k in range(N+1-len(x),N+1):
        env[k-1]=16*((1-x[k-(N+1-len(x))])*x[k-(N+1-len(x))])**2*((1-x[k-(N+1-len(x))])/(2*(1-l2)))
    return env

attenuation = envelop([l1,l2],2*fEch+1)

plt.figure(8)
plt.plot(t,attenuation)
plt.title("Représentation de l'enveloppe")
plt.xlabel("temps")
plt.ylabel("Valeur du son")

#%%Affiche le son amorti sur la période 2sec

SonAtenue = Vp(t,30)*attenuation

plt.figure(9)
plt.plot(t,SonAtenue)
plt.title("Représentation du son amorti")
plt.xlabel("temps")
plt.ylabel("Valeur du son")

#%%Gérérer le son amorti

# script audio.py
# (C) Fabrice Sincère ; Jean-Claude Meilland
import wave
import math
import binascii

print("Signal amorti")
print("Création d'un fichier audio au format WAV (PCM 8 bits stéréo 44100 Hz)")
print("Son de forme sinusoïdale\n")

NomFichier = 'son_amorti.wav'
Monson = wave.open(NomFichier,'w') # instanciation de l'objet Monson

print("Nombre d'échantillons :",nbEchantillon)

parametres = (nbCanal,nbOctet,fech,nbEchantillon,'NONE','not compressed')# tuple
Monson.setparams(parametres)    # création de l'en-tête (44 octets)

# niveau max dans l'onde positive : +1 -> 255 (0xFF)
# niveau max dans l'onde négative : -1 ->   0 (0x00)
# niveau sonore nul :                0 -> 127.5 (0x80 en valeur arrondi)

print('Veuillez patienter...')
for i in range(0,nbEchantillon):
    # canal
    # 127.5 + 0.5 pour arrondir à l'entier le plus proche
    val = wave.struct.pack('B',int(128.0 + 127.5*10000*attenuation[i]*Vp(i/fEch,30)))
    Monson.writeframes(val) # écriture frame

Monson.close()

Fichier = open(NomFichier,'rb')
data = Fichier.read()
tailleFichier = len(data)
print('\nTaille du fichier',NomFichier, ':', tailleFichier,'octets')
print("Lecture du contenu de l'en-tête (44 octets) :")
print(binascii.hexlify(data[0:44]))
print("Nombre d'octets de données :",tailleFichier - 44)
Fichier.close()
print("fini")