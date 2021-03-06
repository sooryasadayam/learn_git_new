from PyGMO import algorithm, island
from PyGMO.problem import base
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import ode
from PyGMO.problem import base
from PyGMO import algorithm, island
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
source = AD('A', 'P', 'AP', 'R1')                #reactant_1,reactant_2, Product
gate = AD('P', 'E', 'PE', 'R2')                  #reactant_1,reactant_2, Product
drain = ES('A', 'E', '~A', 'R3')                 #enzyme,substrate,product
c = model(source, gate)
c.fuse()
c = model(c,drain)
c.fuse()
Names_Rxn = c.Rxn
Names_Rcnt = c.Rcnt
Rmat_dict = {}
#all parameters initiated to zero
for i in range(0, c.Rmat.shape[0]):
        for j in range(0, c.Rmat.shape[1]):
            temp = c.Rmat[i, j]
            if ('+' in temp):
                temp1 = temp.split('+')[0]
                temp2 = temp.split('+')[1]
                if ('-' in temp1):
                    temp1 = temp1.split('-')[1]
                elif ('-' in temp2):
                    temp2 = temp2.split('-')[1]
                if ~(temp1 in Rmat_dict.keys()):
                    Rmat_dict[temp1] = 0
                elif ~(temp2 in Rmat_dict.keys()):
                    Rmat_dict[temp2] = 0
            elif ('-' in temp):
                if ~(temp.split('-')[1] in Rmat_dict.keys()):
                    Rmat_dict[temp.split('-')[1]] = 0
            elif ~(temp in Rmat_dict.keys()):
                Rmat_dict[temp] = 0
    # Original Matrix formula -> c.Rxn = c.Rmat*c.Rcnt
    # dictionary form         -> Rxn_dict,Rmat_dict, Rcnt_dict
#print (Rmat_dict.keys())
def ds_dt(t,S,P):
    global Rmat_dict
    Rmat_dict[""] = P[0]
    Rmat_dict["R1kf"] = P[1]
    Rmat_dict["R1kb"] = P[2]
    Rmat_dict["R2kf"] = P[3]
    Rmat_dict["R2kb"] = P[4]
    Rmat_dict["R3kf"] = P[5]
    Rmat_dict["R3kb"] = P[6]
    Rmat_dict["R3kr"] = P[7]
    global c
    global Names_Rcnt
    global Names_Rxn
    MAT = np.zeros_like(c.Rmat)
    # formula parsing MAT
    for i in range(0, c.Rmat.shape[0]):
        for j in range(0, c.Rmat.shape[1]):
            temp = c.Rmat[i, j]
            if ('+' in temp):
                temp1 = temp.split('+')[0]
                temp2 = temp.split('+')[1]
                if ('-' in temp1):
                    val1 = -1 * Rmat_dict[temp1.split('-')[1]]
                else:
                    val1 = Rmat_dict[temp1]
                if ('-' in temp2):
                    val2 = -1 * Rmat_dict[temp2.split('-')[1]]
                else:
                    val2 = Rmat_dict[temp2]
                MAT[i, j] = val1 + val2
            elif ('-' in temp):
                MAT[i, j] = -1 * Rmat_dict[temp.split('-')[1]]
            else:
                MAT[i, j] = Rmat_dict[temp]
    # Rmat_dict has an unordered parameter mapping
    # params is the ordered set of parameter names
    # params_value is the ordered set of parameter values
    # formula parsing RXN
    rxn_values = np.reshape(S, (7,))  #### This is where the S_vec goes,There's an adhoc Jugad of size changing here
    Rxn_dict = dict(zip([i[0] for i in Names_Rxn], rxn_values))
    rcnt_values = np.zeros_like(Names_Rcnt)
    # formula parsing RCNT
    for i in range(0, Names_Rcnt.shape[0]):
        if ('*' in Names_Rcnt[i, 0]):
            rcnt_values[i, 0] = Rxn_dict[Names_Rcnt[i, 0].split('*')[0]] * Rxn_dict[Names_Rcnt[i, 0].split('*')[1]]
        else:
            rcnt_values[i, 0] = Rxn_dict[Names_Rcnt[i, 0]]
    # print(rxn_values)
    # print(rcnt_values)
    # print(MAT)
    Rcnt_dict = dict(zip([i[0] for i in Names_Rcnt], rcnt_values))
    ### dS_dt is being defined
    dS_dt = np.dot(MAT, rcnt_values)
    return dS_dt
class my_problem(base):
    def __init__(self, dim=7):
        super(my_problem,self).__init__(dim)
        self.set_bounds(0.00, 2.00)
    def _objfun_impl(self, Params):
        global ds_dt
        r = ode(ds_dt)
        r.set_integrator("vode")
        r.set_initial_value(np.ones((7, 1)))
        r.set_f_params(np.hstack((np.array([0]), Params)))  # 7 params
        solution = np.vstack((np.array([0]),np.ones((7,1))))
        T_max = 3.0
        dt = 0.1
        while r.successful() and r.t < T_max - dt:
            solution = np.hstack((solution, np.vstack((np.array([r.t + dt]), r.integrate(r.t + dt)))))
        error = 3 - np.linspace(0.0, 2.9, 30) ** 2 - np.reshape((solution[1, :]), 30)
        sse = np.sum(error ** 2)
        return (sse,)
    def human_readable_extra(self):
        return "\n\t Problem dimension: " + str(self.__dim)
prob = my_problem(dim=7)           # Create a 7-dimensional problem
algo = algorithm.mde_pbx(gen=500)  # 500 generations of mde_pbx algorithm
isl = island(algo, prob, 20)  # Instantiate population with 20 individuals
isl.evolve(1)  # Evolve the island once
isl.join()
print(isl.population.champion)
