from matplotlib import pyplot as plt
import numpy as np
import numpy.linalg as alg
import math

## Exemple Warren

#liste des coordonées des points
Points=[(0,0),(1,0),(2,0),(3,0),(0.5,math.sqrt(3)/2),(1.5,math.sqrt(3)/2),(2.5,math.sqrt(3)/2)] 

#liste des liaisons données par les positions dans Points
Liaisons=[(0,1),(1,2),(2,3),(0,4),(4,1),(1,5),(5,2),(2,6),(6,3),(4,5),(5,6)]       

#liste de listes des forces sur chaque point,n°du point+(array)
Forces_ext=[(1,[0.,-50000.]),(2,[0,-50000.]),(0,[0.,50000.]),(3,[0,50000.])]     

#liste contenant les points de laisons au bati
Bati=[0,3]           

#reste vide au debut, puis on ajoute la norme de la force pour la liaison entre les 2 points
F_li=[]           

## Exemple Warren inversé

Points=[(0,0),(1,0),(2,0),(3,0),(0.5,-math.sqrt(3)/2),(1.5,-math.sqrt(3)/2),(2.5,-math.sqrt(3)/2)]

Liaisons=[(0,1),(1,2),(2,3),(0,4),(4,1),(1,5),(5,2),(2,6),(6,3),(4,5),(5,6)]

Forces_ext=[(1,[0.,-50000.]),(2,[0,-50000.]),(0,[0.,50000.]),(3,[0,50000.])]

Bati=[0,3] 

F_li=[]

## Exemple Warren avec montants

Points=[(0,0),(0.5,0),(1,0),(1.5,0),(2,0),(2.5,0),(3,0),(0.5,math.sqrt(3)/2),(1,math.sqrt(3)/2),(1.5,math.sqrt(3)/2),(2,math.sqrt(3)/2),(2.5,math.sqrt(3)/2)]

Liaisons=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(0,7),(7,1),(7,2),(7,8),(8,9),(9,10),(10,11),(2,8),(2,9),(9,3),(9,4),(4,10),(4,11),(11,5),(11,6)]

Forces_ext=[(2,[0.,-20000.]),(4,[0,-20000.]),(0,[0.,50000.]),(6,[0.,50000.]),(1,[0.,-20000.]),(3,[0.,-20000.]),(5,[0.,-20000.]),]

Bati=[0,6] 

F_li=[]

## Exemple Howe

Points=[(0,0),(0.5,0),(1,0),(1.5,0),(2,0),(2.5,0),(3,0),(0.5,math.sqrt(3)/2),(1,math.sqrt(3)/2),(1.5,math.sqrt(3)/2),(2,math.sqrt(3)/2),(2.5,math.sqrt(3)/2)]

Liaisons=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(0,7),(7,1),(7,8),(8,9),(9,10),(10,11),(1,8),(8,2),(9,2),(9,3),(4,9),(4,10),(10,5),(11,5),(11,6)]

Forces_ext=[(2,[0.,-20000.]),(4,[0,-20000.]),(0,[0.,50000.]),(6,[0,50000.]),(1,[0.,-20000.]),(3,[0.,-20000.]),(5,[0.,-20000.]),]

Bati=[0,6] 

F_li=[]

## Pont de l'A3

Points=[(0,0),(13.34,0),(26.68,0),(40.02,0),(53.36,0),(66.7,0),(80.04,0),(93.38,0),(106.72,0),(120.06,0),(133.4,0),(6.67,11.55),(20.01,11.55),(33.35,11.55),(46.69,11.55),(60.03,11.55),(73.37,11.55),(86.71,11.55),(100.05,11.55),(113.39,11.55),(126.73,11.55)]

Liaisons=[(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20),(0,11),(11,1),(1,12),(12,2),(2,13),(13,3),(3,14),(14,4),(4,15),(15,5),(5,16),(16,6),(6,17),(17,7),(7,18),(18,8),(8,19),(19,9),(9,20),(20,10)]

Bati=[0,10]

Forces_ext=[(1,[0.,-1200000.]),(2,[0.,-1200000.]),(3,[0.,-1200000.]),(4,[0.,-1200000.]),(5,[0.,-1200000.]),(6,[0.,-1200000.]),(7,[0.,-1200000.]),(8,[0.,-1200000.]),(9,[0.,-1200000.]),(0,[0,5400000]),(10,[0,5400000])]

F_li=[] 


##  
def voisins(point):
    L=[]
    for k in Liaisons:
        x,y=k
        if x==point:
            L.append(y)
        elif y==point:
            L.append(x)
    return L

    
def liste_liaisons_inconnues (F_li,point):
    L=voisins(point)
    L1=L
    for i in range (len(F_li)):
        tuple,F=F_li[i]
        x,y=tuple
        if x==point:
            L1.remove(y)
        elif y==point:
            L1.remove(x)
    return L1
            
def etre_dans_bati(point):
    for k in Bati:
        if k==point:
            return True
    return False

def norme(Force):
    return(math.sqrt(Force[0]**2+Force[1]**2))
        

def BDF(point):
    Lf=[]
    for k in range (len(Forces_ext)):
        p,F=Forces_ext[k]
        if p==point:
            Lf.append(F)
    for i in F_li:
        tuple,F=i
        x,y=tuple
        
        if x==point:
            a1,b1=Points[x]
            a2,b2=Points[y]
            hyp=math.sqrt((a2-a1)**2+(b2-b1)**2)
            ct=(a1-a2)/hyp
            st=(b1-b2)/hyp
            Lf.append([F*ct,F*st])
        elif y==point:
            a1,b1=Points[y]
            a2,b2=Points[x]
            hyp=math.sqrt((a1-a2)**2+(b1-b2)**2)
            ct=(a1-a2)/hyp
            st=(b1-b2)/hyp
            Lf.append([F*ct,F*st])
    return (Lf)
        


        
def PFS_1(Lf,x,y):    #au point x
    a1,b1=Points[x]
    a2,b2=Points[y]
    X=0
    Y=0
    for k in Lf:
        X+=k[0]
        Y+=k[1]
    if a1==a2:
        F=Y
    else:
        hyp=math.sqrt((a1-a2)**2+(b1-b2)**2)
        ct=(a1-a2)/hyp
        F=-X/ct
    F_li.append(((x,y),F))
    return (F)
    
    

def PFS_2(Lf,x,y1,y2):
    a,b=Points[x]
    a1,b1=Points[y1]
    a2,b2=Points[y2]
    X=0
    Y=0
    
    for k in Lf:
        X+=k[0]
        Y+=k[1]
        
    hyp1=math.sqrt((a-a1)**2+(b-b1)**2)
    hyp2=math.sqrt((a-a2)**2+(b-b2)**2)
    
    ct1=(a-a1)/hyp1
    ct2=(a-a2)/hyp2
    st1=(b-b1)/hyp1
    st2=(b-b2)/hyp2
    
    F=np.array([-X,-Y])
    
    Force=alg.solve(np.array([[ct1,ct2],[st1,st2]]),F)
    
    F_li.append(((x,y1),Force[0]))
    F_li.append(((x,y2),Force[1]))
    
    return (Force[0],Force[1])
    
    



    
def Cremona(Points,Forces_ext,Bati,F_li):
    n=len(Liaisons)
    if len(F_li)==n:
        for k in Forces_ext:
            x,F=k
            a,b=0,0
            if etre_dans_bati(x):
                p1,p2=Points[x]
                if F[0]>0:
                    a=6
                elif F[0]<0:
                    a=-6
                if F[1]>0:
                    b=10
                elif F[1]<0:
                    b=-6
                Nf=norme(F)
                plt.annotate(Nf, xy=(p1,p2), xytext=(p1-a,p2-b),arrowprops={'facecolor':'black', 'shrink':0.05} )
            
            else:
                p1,p2=Points[x]
                if F[0]>0:
                    a=6
                elif F[0]<0:
                    a=-10
                if F[1]>0:
                    b=6
                elif F[1]<0:
                    b=-6
                Nf=norme(F)
                plt.annotate(Nf, xy=(p1,p2), xytext=(p1-a,p2-b),arrowprops={'facecolor':'yellow', 'shrink':0.05} )
        plt.xlabel('x (metre)')
        plt.ylabel('y (metre)')
        plt.title("Diagramme de Cremona (forces en Newton)")
        plt.axis('equal')
        plt.show()
    
    else:
        for k in range (len(Points)):
            vois=voisins(k)
            nb=len(vois)
            Linc=liste_liaisons_inconnues (F_li,k)
            
            if len(Linc)==1:
                x,y=k,Linc[0]
                Lf=BDF(x)
                F=PFS_1(Lf,x,y)
                color='c'
                if F<0:
                    color='g'
                elif F>0:
                    color='r'
                a1,b1=Points[x]
                a2,b2=Points[y]
                plt.plot([a1,a2],[b1,b2],color,linewidth=3.5)
                #plt.text((a1+a2)/2,(b1+b2)/2,F)
                
                Cremona(Points,Forces_ext,Bati,F_li)
                
            elif len(Linc)==2 :
                x=k
                Lf=BDF(x)
                y1,y2=Linc[0],Linc[1]
                F1,F2=PFS_2(Lf,x,y1,y2)
                
                color='c'
                if F1<0:
                    color='g'
                elif F1>0:
                    color='r'
                a1,b1=Points[x]
                a2,b2=Points[y1]
                plt.plot([a1,a2],[b1,b2],color,linewidth=3.5)
                #plt.text((a1+a2)/2,(b1+b2)/2,F1)
                
                color='c'
                if F2<0:
                    color='g'
                elif F2>0:
                    color='r'
                c1,d1=Points[x]
                c2,d2=Points[y2]
                plt.plot([c1,c2],[d1,d2],color,linewidth=3.5)
                #plt.text((c1+c2)/2,(d1+d2)/2,F2)
                
                Cremona(Points,Forces_ext,Bati,F_li)
                
                    
                
            
def courbe(F_li):   ## Pour le pont de L'A3
    X1=[]
    X2=[]
    X3=[]
    Y1=[]
    Y2=[]
    Y3=[]
    for k in range (len(F_li)):
        F=F_li[k][1]
        pt1,pt2=F_li[k][0]
        if pt1<=10 and pt2<=10:
            x1,y1=Points[pt1]
            x2,y2=Points[pt2]
            X1.append((x1+x2)/2)
            Y1.append(abs(F))
        elif pt1>10 and pt2>10:
            x1,y1=Points[pt1]
            x2,y2=Points[pt2]
            X2.append((x1+x2)/2)
            Y2.append(abs(F))
        else:
            x1,y1=Points[pt1]
            x2,y2=Points[pt2]
            X3.append((x1+x2)/2)
            Y3.append(abs(F))
    tri_insertion(X1,Y1)    
    tri_insertion(X2,Y2)
    tri_insertion(X3,Y3)    
    plt.plot(X1,Y1,'r',label="poutres inférieures") 
    plt.plot(X2,Y2,'b',label="poutres supérieures")
    plt.plot(X3,Y3,'g',label="diagonales") 
    plt.xlabel('Distance (m)')
    plt.ylabel('Force (N)')
    plt.title("Efforts au sein de la structure")
    plt.legend()     
    plt.show()
  
  
def tri_insertion(X,F):
    n=len(X)
    for i in range (1,n):
        x=X[i]
        f=F[i]
        p=i
        while p>0 and X[p-1]>x:
            X[p]=X[p-1]
            F[p]=F[p-1]
            p=p-1
        X[p]=x
        F[p]=f