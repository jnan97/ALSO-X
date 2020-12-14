#!/usr/bin/env python
# coding: utf-8

# In[5]:


import sys
from gurobipy import *
import numpy as np
import math
import time
import copy
from numpy import linalg as LA


# In[6]:


# Set problem parameter
# random_size = int(input('Input the size of the model: '))
# x_random_size = int(input('Input the size of the model of x: '))
random_size=1000
x_random_size=100


# In[7]:


# Set seed
np.random.seed(45)


# In[8]:


# distribution
a = np.random.randint(1,99, size=(random_size,x_random_size))
c = np.random.randint(-10,0,x_random_size)
bb = np.random.normal(100,0,random_size)


# In[9]:


# Set initial values for epsilon
epsilon = 0.05
violation=math.floor(random_size*epsilon)


# In[ ]:





# In[10]:


def proj2(x_hat,x_random_size,c,t):
    x_hat=x_hat.reshape(x_random_size).tolist()
    hatxx = x_hat
    n=x_random_size
    xsol=[1.0]*n

    iter1=0
    Gap=1e20
    bestlb=0
    bestub=1e20

    maxiter=1000
    len_c=float(sum([abs(c)[i]*abs(c)[i] for i in range(n)]))
    check=max(abs(t)-sum([abs(c)[i]*hatxx[i] for i in range(n)]),0)
    xsol=[hatxx[i]+check/len_c*abs(c)[i] for i in range(n)]
    xsol=[max(min(xsol[i],1.0),0.0) for i in range(n)]
    prevub=1e20

    while iter1<=maxiter and Gap>1e-4:
        iter1=iter1+1
        subg=[0.0]*n
        ub=sum([(xsol[l]-hatxx[l])**2.0 for l in range(n) ])

        check=max(abs(t)-sum([abs(c)[i]*xsol[i] for i in range(n)]),0)
        xsol=[xsol[i]+check/len_c*abs(c)[i] for i in range(n)]
        xsol=[max(min(xsol[i],1.0),0.0) for i in range(n)]

        bestub=min(ub, bestub)
        Gap=(prevub-bestub)/(abs(bestub)+1e-5)
        prevub=bestub
        #print("at iteration %f,bestub =%f, bestlb =%f, gap=%f" %(iter1, bestub,bestlb,Gap))
    return xsol


# In[11]:


def proj_d(x_hat,x_random_size,c,t):
    xsol=copy.deepcopy(x_hat)
    # x project on box constraint, y project on knapsack constraint
    ysol=np.zeros((x_random_size,1))
    p=np.zeros((x_random_size,1))
    q=np.zeros((x_random_size,1))
    maxiter=1000

    iter1=0
    Gap=1e20
    bestlb=0
    bestub=1e20
    prevub=1e20
    ub=1e20

    while iter1<=maxiter and Gap>1e-4:

        iter1+=1
        len_c = float(np.sum(np.square(c)))
        check = max(abs(t)-np.dot(abs(c), xsol+p)[0],0)
        ysol=xsol+p+ check/len_c*abs(c).reshape((x_random_size,1))
        p = xsol+p-ysol
        xsol = np.maximum(np.minimum(ysol+q, np.ones((x_random_size,1))),np.zeros((x_random_size,1)))
        q = ysol + q - xsol 
        ub = np.sum(np.square(xsol-x_hat))
        bestub=min(ub, bestub)
        Gap=(prevub-bestub)/(abs(bestub)+1e-5)
        prevub=bestub
    return xsol


# In[12]:


def g(x,ep):
    return np.dot(ep,np.square(x))
    
def dg(x,ep):
    m=np.dot(x,np.ones((1,ep.shape[0]))).transpose()
    return 2*ep*m


# In[13]:


start=time.time()

support=random_size
delta_1= 1e-2
t_1= -5
t_2=-200
t=(t_1+t_2)/float(2)
x_s=np.zeros((x_random_size,1))
x_s = abs(t)/np.sum(abs(c)) * np.ones((x_random_size,1))
# t_last = t

while t_1-t_2>=delta_1:
#     t_last = t
    t=(t_1+t_2)/float(2)
#     x_s = x_s/t_last*t
    x_s = proj_d(x_s,x_random_size,c,t)
    
    subgrad = np.zeros((random_size,x_random_size))
    obj=100
    iteration = 0 

    Difference=1e20

    maxiter=5000


    while iteration<=maxiter:
#     while iteration<=maxiter and Difference>1e-5:
        
        iteration +=1
        g_val=g(x_s,a)
        support=np.count_nonzero(g_val>100+1e-5)
        if support <= violation:
            x_f=x_s
            break
        if support >= 1.15*violation and iteration>=maxiter/10:
            x_f=x_s
            break
        
        extended_g_val=np.dot(g_val,np.ones((1,a.shape[1])))
        subgrad=dg(x_s,a)

        subgrad[extended_g_val<=100+1e-5]=0

        subgrad_coefficient = np.average(subgrad,axis=0)

        gamma=1/iteration
        x_f=np.zeros((x_random_size,1))

        x_hat = x_s - gamma*subgrad_coefficient.reshape((x_random_size,1))

        if np.dot(x_hat.transpose(),c)==t and np.all(x_hat <= 1) and np.all(x_hat >= 0):
            x_f = x_hat
        else:
#             x_val = Proj(c,t,x_hat)
#             x_f = x_val.values()
#             x_f = np.array(np.matrix(x_f).T)
#            x_f = np.array(proj2(x_hat,x_random_size,c,t)).reshape(x_random_size,1)
            x_f = proj_d(x_hat,x_random_size,c,t)
#        if obj >np.average(np.maximum(g(x_f,a)-100,np.zeros((random_size,1)))):
#             obj = np.average(np.maximum(g(x_f,a)-100,np.zeros((random_size,1))))
#             best_sol_f = x_f
#         ub =  np.average(np.maximum(g(x_f,a)-100,np.zeros((random_size,1))))
#         bestub=min(ub, bestub)
#         Gap=(prevub-bestub)/(abs(bestub)+1e-5)
#         prevub=bestub
#         print("at iteration %f,bestub =%f, bestlb =%f, gap=%f" %(iteration, bestub,bestlb,Gap))

#         Difference = abs(np.sum(x_s-x_f))
#         print ('current Difference, ',Difference)


    #     if Difference > abs(np.average(np.maximum(g(x_s,a)-100,np.zeros((random_size,1))))-np.average(np.maximum(g(x_f,a)-100,np.zeros((random_size,1))))):
    #         Difference = abs(np.average(np.maximum(g(x_s,a)-100,np.zeros((random_size,1))))-np.average(np.maximum(g(x_f,a)-100,np.zeros((random_size,1)))))
    #         best_sol_f = x_f
        x_s = x_f

    #     if GAP>(Upper-Lower)/abs(Upper):
    #         GAP = (Upper-Lower)/abs(Upper)
    #         best_sol_f = x_f
    #         best_iteration = iteration    
    

    g_val=g(x_f,a)
    support=np.count_nonzero(g_val>100+1e-5)
    if support<=violation:
        t_1=t
    else:

        t_2=t
        
    print('current t is : ', t)
    print('current support is: ', support)


modeltime = time.time() - start


# In[14]:


print('G-D')
print(f'The modeltime is {modeltime:.3f} s.')
print(f'The objective s is {t_1:.5f}.')


# In[ ]:





# In[ ]:





# In[16]:


def AP(t,x_s,iteration):
    delta_d = 10000
    #delta_2 = 1e-3*CVAR
    delta_2 = 1e-2
    value = 1000
    subgrad = np.zeros((random_size,x_random_size))
    # subgrad starting 
    iteration_1 = 0 


#     while delta_d >= delta_2 and iteration_1<=10 :
    while delta_d >= delta_2:

        value_1 = value
        iteration_1 +=1
        g_val=g(x_s,a)

        kkkkk=np.sort(g_val,axis=None) 
        g_val[g_val>kkkkk[int(random_size*(1-epsilon))-1]]=0

        value = np.average(g_val)
    #     print('current value  is ', value)

        if  value <= 1e-4:
            delta_d = 0
            break

        delta_d = abs( value - value_1 )
#         print('current difference  is ', delta_d)

        extended_g_val=np.dot(g_val,np.ones((1,a.shape[1])))

        subgrad=dg(x_s,a)

        subgrad[extended_g_val<100+1e-5]=0
        
        subgrad_coefficient = np.average(subgrad,axis=0)

        gamma=1/iteration_1
        x_f=np.zeros((x_random_size,1))

        x_hat = x_s - gamma*subgrad_coefficient.reshape((x_random_size,1))

        if np.dot(x_hat.transpose(),c)==t and np.all(x_hat <= 1) and np.all(x_hat >= 0):
            x_f = x_hat
        else:
    #             x_val = Proj(c,t,x_hat)
    #             x_f = x_val.values()
    #             x_f = np.array(np.matrix(x_f).T)
            x_f = np.array(proj2(x_hat,x_random_size,c,t)).reshape(x_random_size,1)

        x_s = x_f

    g_val=g(x_s,a)
    
    support=np.count_nonzero(g_val>100+1e-5)
    
    return support


# In[17]:


support=random_size
delta_1= 1e-3
t_1= -5
t_2=-200


start=time.time()
x_s=np.zeros((x_random_size,1))
x_s = abs(t)/np.sum(abs(c)) * np.ones((x_random_size,1))
t_last = t

while t_1-t_2>=delta_1:
    
    t_last = t
    t=(t_1+t_2)/float(2)
    x_s = x_s/t_last*t
    
    g_val=g(x_s,a)
    support=np.count_nonzero(g_val>100+1e-5)

    subgrad = np.zeros((random_size,x_random_size))
    obj=100
    iteration = 0 

    Difference=1e20
    maxiter = 5000
 

    while iteration<=maxiter:
    
        iteration +=1
        g_val=g(x_s,a)
        support=np.count_nonzero(g_val>100+1e-5)
        
        if support <= violation:
            x_f=x_s
            break
        if support >= 1.15*violation and iteration>=maxiter/10:
            x_f=x_s
            break


        extended_g_val=np.dot(g_val,np.ones((1,a.shape[1])))
        subgrad=dg(x_s,a)

        subgrad[extended_g_val<=100+1e-5]=0

        subgrad_coefficient = np.average(subgrad,axis=0)

        gamma=1/iteration
        x_f=np.zeros((x_random_size,1))

        x_hat = x_s - gamma*subgrad_coefficient.reshape((x_random_size,1))

        if np.dot(x_hat.transpose(),c)==t and np.all(x_hat <= 1) and np.all(x_hat >= 0):
            x_f = x_hat
        else:

            x_f = np.array(proj2(x_hat,x_random_size,c,t)).reshape(x_random_size,1)


        x_s = x_f



    g_val=g(x_f,a)
    support=np.count_nonzero(g_val>100+1e-5)
    if support<=violation:
        t_1=t
        print('current t is : ', t)
        print('current support is: ', support)
    else:
        print('current t is : ', t)
        print('current iteration  is ', iteration)
        print('current support is: ', support)
        AP_checking=AP(t,x_f,iteration)

        if AP_checking<=violation:
            t_1=t

        else:
            t_2=t

        print('AP current t is : ', t)
        print('AP current support is: ', AP_checking)
    

modeltime = time.time() - start


# In[ ]:





# In[18]:


print('G-D - AP')
print(f'The modeltime is {modeltime:.3f} s.')
print(f'The objective s is {t_1:.5f}.')


# In[ ]:





# In[ ]:


sample_size=1


# In[ ]:


# CVaR Method

def Model_cvar():
    #Model
    print ("Begin to solve CVaR")

    m = Model()
    #m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    s = m.addVars(random_size,lb=-GRB.INFINITY,name="s")
  #  z = m.addVars(random_size,name="z")
    beta = m.addVar(ub=0,lb=-GRB.INFINITY,name="beta")
    m.update()
    
    # Set functions
    
    #obj = 0

    m.setObjective(quicksum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
        
    m.addConstrs(quicksum(a[i,j]*x[j]*x[j] for j in range(x_random_size))  <= bb[i] + s[i]  for i in range(random_size))
    m.addConstrs(s[i] >= beta for i in range(random_size))
    m.addConstr(s.sum()/random_size -beta*(1-epsilon) <= 0)
  
    # Solve the problem
    m.update()
        
    m.optimize()
    

#     Store solutions

    CVAR_s = m.getAttr('x', s)
    CVAR_x = m.getAttr('x', x)
    

    return m.objVal,CVAR_s,CVAR_x


    


# In[ ]:


CVaR=Model_cvar()


# In[26]:





# In[30]:


# Big-M Method
Big_M = []
for i in range(random_size):
   
    Big_M.append(math.ceil(sum(x for x in a[i,:] if x > 0) -bb[i])) 

def Model3():
    #Model
    print ("Begin to solve model 3")

    m = Model()
    m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()
    # s1 = m.addVars(random_size,name="s1")
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    z = m.addVars(random_size,vtype=GRB.BINARY,name="z")
    m.update()
    
    # Set functions
    
    #obj = 0
    g_x_xi = [None] * len(z)
    f_x = 0
    g_x = 0
    
    # Set functions
    for i in range(x_random_size):
      
        f_x += c[i]*x[i]
        m.update()
    
#     for i in range(random_size):
#         g_x_si[i] = g_x +       SSS   
    for i in range(random_size):
     #   g_x_si[i] = g_x + si[i] 
        g_x = 0
        for j in range(x_random_size):
            g_x += a[i,j]*x[j]*x[j]
            m.update()
        g_x_xi[i] = g_x  
    
    # Add objective    
    m.setObjective(f_x, GRB.MINIMIZE)
    m.update()
    # Add constraints
   
    m.addConstrs(g_x_xi[i] <= bb[i] + Big_M[i]*(1-z[i]) for i in range(random_size))

    m.addConstr(z.sum() >= math.ceil(random_size*(1-epsilon)))
  
    # Solve the problem
    m.update()
        
    m.optimize()
    
    # Store solutions
    kk = m.getAttr('x', x)
    
    #return pp,ii,kk
    
    return kk


# In[31]:


start=time.time()
kk=Model3()
obj = 0
for i in range(x_random_size):
    obj += c[i]*kk[i]
modeltime = time.time() - start
print(f'The value of optimal value t is exactly {obj:.5f}.')
print(f'The value of running time is approximately {modeltime:.3f} s.')


# In[ ]:





# In[ ]:




