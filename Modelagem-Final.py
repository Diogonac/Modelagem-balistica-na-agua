"""
Insper
Author: Diogo Nobre de Araujo Cintra
Graduated in Mechanical Engineering
Email: diogonac@al.insper.edu.br

"""
from scipy.integrate import odeint
import numpy as np
import math
import matplotlib.pyplot as plt

t_Ar=np.arange(0,50,1e-5)#Lista de tempo para o Ar
t_Agua=np.arange(0,0.5,1e-5)#Lista de tempo para a Agua
r = 0.00278 #Raio do projetil 2,.78mm 
A=math.pi*r**2#Area de secao transversal do projeti m^2
d = 7850#Densidade do material do projetil kg/m^3
g = 10#Aceleracao da gravidade m/s^2
m = 0.0041#Massa do projetil m
volume = (4/3)*math.pi*r**3#Volume do projetil m^3
p_Ar = 1.2928#Densidade do Ar
visc_Ar = 349.9#Viscosidade cinematica do Ar
p_Agua = 997#Densidade da Agua
visc_Agua = 6.48e3#Viscosidade cinematica da Agua

tempo = [0.49,49.9]#Primeiro na agua depois no ar
dist_X =[11.7,345.4]
meio = [6.48e3,349.9]

Tempo = [0.49,0.880]
Distancia = [11.7,0.720]

def EquacoesDiferenciais_Ar(listaSolucoes, t_Ar): 
    X=listaSolucoes[0]
    Y=listaSolucoes[1]
    Vx=listaSolucoes[2]
    Vy=listaSolucoes[3] 
    
    dXdt=Vx
    dYdt=Vy
    
    V = math.sqrt(Vx**2+Vy**2)#Calculo da velociade 
    Re = V*visc_Ar#Calculo do numero de Reynolds (Re)
    
    #Aproximacao das equacaos que descreve o comportamento do Cd em funcao do numero de Reynolds (Re)
    if V <57:
        #Cd = -0.0007*Re+15.001
        Cd = 10**((-0.7503*(math.log10(Re)+1)+3.9285)-1)

    elif 57<=V<572:
        #Cd = 0.47
        Cd = 10**(-0.39265)
    
    elif 572<=V<2143:
        #Cd = -6e-7*Re+0.5164
        Cd = 10**(((-1.3846*(math.log10(Re)+1))+9.3732)-1)

    #Condicao para parar no chao   
    if Y <= 0:
        dXdt=0
        dYdt=0
        dVxdt=0
        dVydt=0 
        
    dVxdt= (-1/2)*(1/m)*p_Ar *Cd*A*V*Vx
    dVydt= (1/2)*(1/m)*p_Ar*Cd*A*V*-Vy +(1/m)*p_Ar*volume*g -g
    
    return dXdt, dYdt, dVxdt, dVydt,

def EquacoesDiferenciais_Agua(listaSolucoes, t_Agua):
    X=listaSolucoes[0]
    Y=listaSolucoes[1]
    Vx=listaSolucoes[2]
    Vy=listaSolucoes[3]    
    
    dXdt=Vx
    dYdt=Vy
     
    V = math.sqrt(Vx**2+Vy**2)#Calculo da velociade 
    Re = V*visc_Agua#Calculo do numero de Reynolds (Re)
    
    #Aproximacao das equacaos que descreve o comportamento do Cd em funcao do numero de Reynolds (Re)    
    if V <3:
        #Cd = (-0.0007*Re+15.001)*10
        Cd = 10**(((-0.7503*(math.log10(Re)+1))+3.9285)-1)
        
    elif 3<=V<31:
        #Cd = (0.4)*10
        #Cd = 10**(-0.39265)
        Cd = 0.4
        
    elif 31<=V<116:
        #Cd = (-6e-7*Re+0.5164)*10
        Cd = 10**(((-1.3846*(math.log10(Re)+1))+9.3732)-1)
        
    elif 116<V:
       #Cd =  (1e-8*Re+0.0703)*10
       Cd = 10**(((0.4125*(math.log10(Re)+1))-2.9203)-1)

    #Condicao para parar no chao       
    if Y <= 0:
        dXdt=0
        dYdt=0
        dVxdt=0
        dVydt=0    
 
    dVxdt= (-1/2)*(1/m)*p_Agua*Cd*A*V*Vx  
    dVydt= (1/2)*(1/m)*p_Agua*Cd*A*V*-Vy +(1/m)*p_Agua*volume*g -g
    
    return dXdt, dYdt, dVxdt, dVydt,

Vo=911 #Velocidade inicial do projetil m/s
angulo=math.radians(0)#Angulo inicial do projetil rad
cosseno = math.cos(angulo)#Cosseno inicial do projetil rad
seno = math.sin(angulo)#Seno inicial do projetil rad

S0=[0, 1, cosseno*Vo, seno*Vo]#Condicao inicial

S_Ar=odeint(EquacoesDiferenciais_Ar, S0, t_Ar)
S_Agua=odeint(EquacoesDiferenciais_Agua, S0,t_Agua)



plt.plot(S_Ar[:,0], S_Ar[:,1], 'b',label=("Ar"))
plt.title('Alcance do lançamento horizontal de uma projétil de 2,78 mm de raio')
plt.xlabel("Deslocamento horizontal (m)")
plt.ylabel("Deslocamento vertical (m)")
plt.legend()
plt.grid(True)
plt.show()

plt.plot(t_Ar, S_Ar[:,2], 'b',label=("Ar"))
plt.title('Alcance do lançamento horizontal de uma projétil de 2,78 mm de raio')
plt.xlabel("Tempo (s)")
plt.ylabel("Velociade horizontal (m/s)")
plt.legend()
plt.grid(True)
plt.show()

plt.plot(S_Agua[:,0], S_Agua[:,1],'r', label=("Água"))
plt.title('Alcance do lançamento horizontal de uma projétil de 2,78 mm de raio')
plt.xlabel("Deslocamento horizontal (m)")
plt.ylabel("Deslocamento vertical (m)")
plt.legend() 
plt.grid(True)
plt.show()

plt.plot(t_Agua, S_Agua[:,2], 'r',label=("Água"))
plt.title('Alcance do lançamento horizontal de uma projétil de 2,78 mm de raio')
plt.xlabel("Tempo (s)")
plt.ylabel("Velociade horizontal (m/s)")
plt.legend()
plt.grid(True)
plt.show()

for e in range(len(tempo)):    
    plt.scatter(dist_X[e], tempo[e],label=("Viscosidade Cinemática=%.1f" %(meio[e])))

    
plt.title('Alcance máximo do projetil')
plt.xlabel("Deslocamento horizontal (m)")
plt.ylabel("Tempo (s)")
plt.legend()
plt.grid(True)
plt.show()

for i in range(len(Tempo)):    
    plt.scatter(Distancia[i], Tempo[i],label=("Distância máxima=%.1f" %(Distancia[i])))
    
plt.title('Alcance do lançamento horizontal Esperado x Simulação')
plt.xlabel("Deslocamento horizontal (m)")
plt.ylabel("Tempo (s)")
plt.legend()
plt.grid(True)
plt.show()

for i in range(len(Tempo)):    
    plt.scatter(Distancia[i], Tempo[i],label=("Distância máxima=%.1f" %(Distancia[i])))
    
plt.title('Alcance do lançamento horizontal Esperado x Simulação')
plt.xlabel("Deslocamento horizontal (m)")
plt.ylabel("Tempo (s)")
plt.legend()
plt.grid(True)
plt.show()