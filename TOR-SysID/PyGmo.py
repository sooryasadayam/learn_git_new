#from PyGMO import algorithm, island
#from PyGMO.problem import base
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import ode
class AD(object):
    def __init__(self, left, right, product, n):
        self.l = left
        self.r = right
        self.p = product
        self.n = n
        self.Rmat = np.array([['-' + self.n + 'kf', self.n + 'kb'],
                              ['-' + self.n + 'kf', self.n + 'kb'],
                              [self.n + 'kf', '-' + self.n + 'kb']],dtype=object)
        self.Rcnt = np.array([[self.l + '*' + self.r], [self.p]],dtype=object)
        self.Rxn = np.array([[self.l], [self.r], [self.p]],dtype=object)

class ES(object):
     def __init__(self, enzyme, substrate, product,n):
         self.enz = enzyme
         self.sub = substrate
         self.product = product
         self.n = n
         self.Rmat = np.array([['-'+self.n+'kf', self.n+'kb'],
                           ['-'+self.n+'kf', self.n+'kb' + '+' +self.n + 'kr'],
                           [self.n+'kf', '-'+self.n+'kb' +'+' +'-' +self.n + 'kr'],
                           ['', self.n+'kr']],dtype=object)

         self.Rcnt = np.array([[self.enz + '*'+self.sub],[self.enz+'-'+self.sub]],dtype=object)
         self.Rxn  = np.array([[self.enz],[self.sub],[self.enz+'-'+self.sub],[self.product]],dtype=object)

def join_mat(a,b):
    c = np.array(np.tile('', (a.Rmat.shape[0] + b.Rmat.shape[0], a.Rmat.shape[1] + b.Rmat.shape[1])),
                 dtype=object)
    n = a.Rmat.shape[0]
    m = a.Rmat.shape[1]
    c[0:n, 0:m] = a.Rmat
    c[n:, m:] = b.Rmat
    return c

def join_Rcnt(a,b):
     return np.hstack((a.Rcnt[0:,0],b.Rcnt[0:,0]))

def join_Rxn (a,b):
     return np.hstack((a.Rxn[0:,0],b.Rxn[0:,0]))

class model(object):
    def __init__(self,a,b):
        self.a = a
        self.b = b
        self.Rxn  = np.vstack((self.a.Rxn,self.b.Rxn))
        self.Rcnt = np.vstack((self.a.Rcnt,self.b.Rcnt))
        c = np.array(np.tile('', (self.a.Rmat.shape[0] + self.b.Rmat.shape[0], self.a.Rmat.shape[1] + self.b.Rmat.shape[1])),
                     dtype=object)
        n = self.a.Rmat.shape[0]
        m = self.a.Rmat.shape[1]
        c[0:n, 0:m] = self.a.Rmat
        c[n:, m:] = self.b.Rmat
        self.Rmat = c
    def fuse(self):
        delete = []
        Delete = []
        n=self.a.Rxn.shape[0]
        m=self.Rxn.shape[0]
        for i in range(n,m):
            for j in range(0,self.a.Rxn.shape[0]):
                if (self.Rxn[i] ==self.Rxn[j]):
                    self.Rmat[j,0:]=self.Rmat[j,0:]+self.Rmat[i,0:]
                    delete=delete+[i]
        self.Rxn = np.delete(self.Rxn,delete,axis=0)
        self.Rmat= np.delete(self.Rmat,delete,axis=0)
        for i in range(self.a.Rcnt.shape[0],self.Rcnt.shape[0]):
            for j in range(0,self.a.Rcnt.shape[0]):
                if (self.Rcnt[i] ==self.Rcnt[j]):
                    self.Rmat[0:,j]=self.Rmat[0:,j]+self.Rmat[0:,i]
                    Delete = Delete + i
        self.Rcnt=np.delete(self.Rcnt,Delete,axis=0)
        self.Rmat=np.delete(self.Rmat,Delete,axis=1)

R2=AD('tap42p','pph21','tap42p-pph21','R2')
R4=AD('cdc55','pph21','pp2a1','R4')
R5=AD('tap42p','sit4','tap42p-sit4','R5')
R7=AD('sit4','sap','pp2a2','R7')
R11=AD('tip41','tap42','tip41-tap42','R11')
R12=AD('fpr1','rapamycin','fpr1-rapamycin','R12')
R13=AD('tor1/2','fpr1-rapamycin','tor1/2-fpr1-rapamycin','R13')

R1=ES('tor1/2','tap42','tap42p','R1')
R3=ES('pp2a1','tap42p','tap42','R3')
R6=ES('pp2a2','tap42p','tap42','R6')
R8=ES('tor1/2','tip41','tip41p','R8')
R9=ES('pp2a1','tip41p','tip41','R9')
R10=ES('pp2a2','tip41p','tip41','R10')

c=model(R1,R2)
c.fuse()
for i in range(3,14):
    exec("c=model(c,R%d)" % i)
    c.fuse()
Names_Rxn  = c.Rxn
Names_Rcnt = c.Rcnt
Rxn_dict={Names_Rxn[i,0]:i for i in range(0,len(Names_Rxn))}
Rcnt_dict={}
for i in range(0,len(Names_Rcnt)):
    if('*' in Names_Rcnt[i,0]):Rcnt_dict[Names_Rcnt[i,0]] = Rxn_dict[Names_Rcnt[i,0].split('*')[0]]*Rxn_dict[Names_Rcnt[i,0].split('*')[1]]
    else:
        Rcnt_dict[Names_Rcnt[i,0]]=Rxn_dict[Names_Rcnt[i,0]]
Rmat_dict={}
for i in range(0,c.Rmat.shape[0]):
    for j in range(0,c.Rmat.shape[1]):
        temp=c.Rmat[i,j]
        if ('+' in temp):
            temp1=temp.split('+')[0]
            temp2=temp.split('+')[1]
            if('-' in temp1):temp1=temp1.split('-')[1]
            elif('-' in temp2):temp2=temp2.split('-')[1]
            if ~(temp1 in Rmat_dict.keys()):Rmat_dict[temp1]=0
            elif ~(temp2 in Rmat_dict.keys()):Rmat_dict[temp2]=0
        elif ('-' in temp):
            if ~(temp.split('-')[1] in Rmat_dict.keys()):
                   Rmat_dict[temp.split('-')[1] ] = 0
        elif ~(temp in Rmat_dict.keys()):
             Rmat_dict[temp]=0

#Original Matrix formula -> c.Rxn = c.Rmat*c.Rcnt
#dictionary form         -> Rxn_dict,Rmat_dict, Rcnt_dict
params = Rmat_dict.keys()
params_value = [Rmat_dict[i] for i in params]
MAT = np.zeros_like(c.Rmat)
params_value=np.array(range(0,len(params_value)))   #### This is where the parameter vector is assigned
Rmat_dict=dict(zip(params,params_value))
# formula parsing MAT
for i in range(0,c.Rmat.shape[0]):
    for j in range(0,c.Rmat.shape[1]):
        temp=c.Rmat[i,j]
        if ('+' in temp):
            temp1=temp.split('+')[0]
            temp2=temp.split('+')[1]
            if('-' in temp1):val1=-1*Rmat_dict[temp1.split('-')[1]]
            else:val1=Rmat_dict[temp1]
            if('-' in temp2):val2=-1*Rmat_dict[temp2.split('-')[1]]
            else:val2 = Rmat_dict[temp2]
            MAT[i,j]=val1+val2
        elif ('-' in temp):
            MAT[i,j] = -1*Rmat_dict[temp.split('-')[1]]
        else:
            MAT[i,j] = Rmat_dict[temp]
# Rmat_dict has an unordered parameter mapping
# params is the ordered set of parameter names
# params_value is the ordered set of parameter values
# formula parsing RXN
rxn_values = np.array(range(0,len(params_value)))     #### This is where the S_vec goes
Rxn_dict = dict(zip([i[0] for i in Names_Rxn],rxn_values))
rcnt_values=np.zeros_like(Names_Rcnt)
#formula parsing RCNT
for i in range(0,Names_Rcnt.shape[0]):
    if('*' in Names_Rcnt[i,0]):
        rcnt_values[i,0]=Rxn_dict[Names_Rcnt[i,0].split('*')[0]]*Rxn_dict[Names_Rcnt[i,0].split('*')[1]]
    else:
        rcnt_values[i,0]=Rxn_dict[Names_Rcnt[i,0]]
print(rxn_values)
print()
Rcnt_dict=dict(zip([i[0] for i in Names_Rcnt],rcnt_values))
### dS_dt is being defined
dS_dt=np.dot(MAT,rcnt_values)
print(params)
print(Names_Rcnt)
print(Names_Rxn)
print(c.Rmat)


