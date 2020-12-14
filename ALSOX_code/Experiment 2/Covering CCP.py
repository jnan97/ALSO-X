#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from gurobipy import *
import numpy as np
import math
import time
import copy
from numpy import linalg as LA


# In[2]:


# Set parameters for the problem
# random_size = int(input('Input the size of the model: '))
# x_random_size = int(input('Input the size of the model of x: '))
random_size = 400
x_random_size  = 20
sample_size = 1


# In[3]:


# set seed
np.random.seed(45)


# In[4]:


# distribution
a = np.random.randint(1,50, size=(sample_size,random_size,x_random_size,))
c = np.random.randint(1,10,x_random_size)
bb = np.random.normal(40,0,random_size)


# In[5]:


# Set initial values for epsilon
epsilon = 0.05
violation=math.floor(random_size*epsilon)


# In[6]:


# CVaR Method

def Model_cvar():
    #Model
    print ("Begin to solve model 4")

    m = Model()
    m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()

    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    s = m.addVars(random_size,lb=-GRB.INFINITY,name="s")

    beta = m.addVar(lb=-GRB.INFINITY,ub=0,name="beta")
    m.update()
    
    # Set functions
    
    #obj = 0
    g_x_xi = [None] * len(s)

    z = [None] * len(s)
    f_x = 0
    g_x = 0
    g_x_k = [None] * sample_size
    
    # Set functions
    for i in range(x_random_size):
      
        f_x += c[i]*x[i]
        m.update()

        
    for k in range(sample_size):
        g_x_xi = [None] * random_size
        for i in range(random_size):
         #   g_x_si[i] = g_x + si[i] 
            g_x = 0
            
            for j in range(x_random_size):
                g_x += a[k,i,j]*x[j]
                m.update()
            g_x_xi[i] = g_x
        g_x_k[k]=g_x_xi
    
    # Add objective    
    m.setObjective(f_x, GRB.MINIMIZE)
    m.update()
    # Add constraints
    
   
    m.addConstrs(g_x_k[j][i] >= bb[i] - s[i]  for j in range(sample_size) for i in range(random_size))
    m.addConstrs(s[i] >= beta for i in range(random_size))

    
    m.addConstr(s.sum()/random_size -beta*(1-epsilon) <= 0)
  
    # Solve the problem
    m.update()
        
    m.optimize()
    

#     Store solutions

    CVAR_s = m.getAttr('x', s)
    CVAR_x = m.getAttr('x', x)
    

    aaaa=m.objVal
    return m.objVal,CVAR_s,CVAR_x,m.Runtime
    
    


# In[7]:


CVaR=Model_cvar()
CVAR=CVaR[0]
CVAR_s=CVaR[1]
CVAR_x=CVaR[2]
CVAR_time=CVaR[3]


# In[ ]:





# In[8]:


# relaxation
print ("Begin to solve relaxation")

m = Model()

# Create variables
m.update()

x = m.addVars(x_random_size,lb=0,ub=1,name="x")
z = m.addVars(random_size,lb=0,ub=1,name="z")

m.update()

# Set functions

obj = 0
g_x_xi = [None] * random_size
g_x_k = [None] * sample_size
f_x = 0

    
for k in range(sample_size):
    g_x_xi = [None] * random_size
    for i in range(random_size):
         #   g_x_si[i] = g_x + si[i] 
        g_x = 0
            
        for j in range(x_random_size):
            g_x += a[k,i,j]*x[j]
            m.update()
        g_x_xi[i] = g_x
    g_x_k[k]=g_x_xi

# Add objective  
m.update()
m.setObjective(quicksum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
# Add constraints

m.update()
m.addConstrs(g_x_k[k][i] >= bb[i]*z[i] for k in range(sample_size) for i in range(random_size))
m.addConstr(z.sum()>=random_size-violation)

# Solve the problem
#m.params.OutputFlag=0 
m.optimize()
relaxation_runtime=m.Runtime


# Store solutions

x_value_relaxation = m.getAttr('x', x)

value_rel=sum(x_value_relaxation.values()[i]*c[i] for i in range(x_random_size))


# In[ ]:





# In[ ]:





# In[9]:


# alsox
print ("Begin to solve alsox")

start=time.time()

# Set initial values
delta_1=1e-2
t1 =  value_rel
t2 =  (violation+1)*t1
delta_t = t2-t1
t = (t1+t2)/2.0
violation=math.floor(random_size*epsilon)

m = Model()

# Create variables
m.update()
s = m.addVars(random_size,lb=0,name="s")
x = m.addVars(x_random_size,lb=0,ub=1,name="x")

m.update()

# Set functions

obj = 0
g_x_xi = [None] * len(s)
g_x_k = [None] * sample_size
f_x = 0
# g_x = 0


o = [1] * random_size
# Set functions
for i in range(random_size):
    obj += s[i]*o[i]
    m.update()

for i in range(x_random_size):   
    f_x += c[i]*x[i]
    m.update()


    
for k in range(sample_size):
    g_x_xi = [None] * random_size
    for i in range(random_size):
         #   g_x_si[i] = g_x + si[i] 
        g_x = 0
            
        for j in range(x_random_size):
            g_x += a[k,i,j]*x[j]
            m.update()
        g_x_xi[i] = g_x
    g_x_k[k]=g_x_xi

# Add objective    
m.setObjective(obj, GRB.MINIMIZE)
# Add constraints

m.addConstrs(g_x_k[k][i] >= bb[i] - s[i]  for k in range(sample_size) for i in range(random_size))
con1=m.addConstr(f_x <= t)

m.params.OutputFlag=0 
m.optimize()

# Store solutions
ii = m.getAttr('x', s)
kkk = m.getAttr('x', x)
iteration = 0

while delta_t >= delta_1: 
    
    support=0
    for j in ii:
        if ii[j]!=0:
            support +=1
    if support > violation:
        t1=con1.SARHSUp
    else:
        t2=con1.SARHSLow
    t=(t1+t2)/2.0
    delta_t=t2-t1
    con1.RHS = t
    for i in range(random_size):
        s[i].start = ii[i]
    for i in range(x_random_size):
        x[i].start = kkk[i] 
    m.params.OutputFlag=0 
    m.optimize()
    iteration+=1
    ii = m.getAttr('x', s)
    kkk = m.getAttr('x', x)
con1.RHS = t2
m.params.OutputFlag=0 
m.optimize()
ii = m.getAttr('x', s)
kkk = m.getAttr('x', x)
heuristics_modeltime = time.time() - start
heuristics_value=t2


# In[10]:


# Relax & Scale
print ("Begin to solve Relax & Scale")
start=time.time()


# Set initial values
x_bar=[None] * x_random_size
z_bar=[None] * len(s)
delta_1=1e-2
lower =  1
upper =  (violation+1)
delta_error = upper-lower
tau = (lower+upper)/2.0

m = Model()

# Create variables
m.update()

x = m.addVars(x_random_size,lb=0,ub=1,name="x")
z = m.addVars(random_size,lb=0,ub=1,name="z")

m.update()

# Set functions

obj = 0
g_x_xi = [None] * random_size
g_x_k = [None] * sample_size
f_x = 0



    
for k in range(sample_size):
    g_x_xi = [None] * random_size
    for i in range(random_size):
         #   g_x_si[i] = g_x + si[i] 
        g_x = 0
            
        for j in range(x_random_size):
            g_x += a[k,i,j]*x[j]
            m.update()
        g_x_xi[i] = g_x
    g_x_k[k]=g_x_xi

# Add objective  
m.update()
m.setObjective(quicksum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
# Add constraints

m.update()
m.addConstrs(g_x_k[k][i] >= bb[i]*z[i] for k in range(sample_size) for i in range(random_size))
m.addConstr(z.sum()>=random_size-violation)

# Solve the problem
m.params.OutputFlag=0 
m.optimize()

# Store solutions

x_value_relaxation = m.getAttr('x', x)
z_value_relaxation = m.getAttr('x', z)

while delta_error >= delta_1: 
    
    for i in range(x_random_size):
          x_bar[i] = x_value_relaxation.values()[i]*tau
    for i in range(random_size):
  
        z_bar[i] = min(math.floor(z_value_relaxation.values()[i]*tau),1)
    if sum(z_bar)>=random_size-violation:
        upper = tau
    else:
        lower=tau  
    tau = (lower+upper)/2.0
    delta_error = upper-lower
    
modeltime = time.time() - start

for i in range(x_random_size):
      x_bar[i] = x_value_relaxation.values()[i]*upper
for i in range(random_size):

    z_bar[i] = min(math.floor(z_value_relaxation.values()[i]*upper),1)

modeltime = time.time() - start

value_scale=sum(x_bar[i]*c[i] for i in range(x_random_size))


# In[ ]:





# In[ ]:





# In[11]:


def AP(t,ii,kkk):
# def AP(t):   
    o = [1] * random_size
    delta_d = 10000
    #delta_2 = 1e-3*CVAR
    delta_2 = 1e-3
    value = 1000
    
    m = Model()

    # Create variables
    m.update()
    s = m.addVars(random_size,lb=0,name="s")
    x = m.addVars(x_random_size,ub=1,lb=0,name="x")
  
    m.update()
    
    # Set functions
    
    obj = 0
    g_x_xi = [None] * random_size
    g_x_k = [None] * sample_size
    f_x = 0
   # g_x = 0
    # Set functions
    for i in range(random_size):
        obj += s[i]*o[i]
        m.update()

    for i in range(x_random_size):   
        f_x += c[i]*x[i]
        m.update()
    
#     for i in range(random_size):
#      #   g_x_si[i] = g_x + si[i] 
#         g_x = 0
#         for j in range(x_random_size):
#             g_x += a[i,j]*x[j]
#             m.update()
#         g_x_xi[i] = g_x 
        
    for k in range(sample_size):
        g_x_xi = [None] * random_size
        for i in range(random_size):
         #   g_x_si[i] = g_x + si[i] 
            g_x = 0
            
            for j in range(x_random_size):
                g_x += a[k,i,j]*x[j]
                m.update()
            g_x_xi[i] = g_x
        g_x_k[k]=g_x_xi
        
    # Add objective    
    m.setObjective(obj, GRB.MINIMIZE)
    # Add constraints
    m.addConstr(f_x <= t)
    m.addConstrs(g_x_k[k][i] >= bb[i] - s[i]  for k in range(sample_size) for i in range(random_size))
    
#     for i in range(random_size):
#         s[i].start = CVAR_s[i]
#     for i in range(x_random_size):
#         x[i].start = CVAR_x[i] 
    # Solve the problem
    
    m.params.OutputFlag=0 
    m.optimize()
    
    # Store solutions
    i_1 = m.getAttr('x', s)
    k_1 = m.getAttr('x', x)
    
    iteration_ap = 0
    
    while delta_d >= delta_2 and iteration_ap <= 5 :
        
        iteration_ap += 1 
        value_1 = value
        o = [0] * random_size
        in_arr=np.array(i_1.values())
        out_arr = np.argsort(in_arr) 
        for i in range(int(random_size*(1-epsilon))):
            o[out_arr[i]]=1
            
        value = sum(i_1[i]*o[i] for i in range(random_size))
        
        if  value == 0:
            delta_d = 0
            break
        delta_d = abs( value - value_1 )
        m.NumObj=0
        obj = 0
        for i in range(random_size):
            obj += s[i]*o[i]
            m.update()
        m.setObjective(obj, GRB.MINIMIZE)
        for i in range(random_size):
            s[i].start = i_1[i]
        for i in range(x_random_size):
            x[i].start = k_1[i]  
        m.params.OutputFlag=0 
        m.optimize()

        i_1 = m.getAttr('x', s)
        k_1 = m.getAttr('x', x)
        
    s_l_0 = 0
    for i in range(random_size):
        if i_1[i] != 0:
            s_l_0 +=1

    return s_l_0
 
    


# In[12]:


# alsox+
print ("Begin to solve alsox+ ")
start=time.time()


# Set initial values
delta_1=1e-4
t1 =  value_rel
t2 =  (violation+1)*t1
delta_t = t2-t1
t = (t1+t2)/2.0

m = Model()

# Create variables
m.update()
s = m.addVars(random_size,lb=0,name="s")
x = m.addVars(x_random_size,lb=0,ub=1,name="x")

m.update()

# Set functions

obj = 0
g_x_xi = [None] * len(s)
g_x_k = [None] * sample_size
f_x = 0
# g_x = 0


o = [1] * random_size
# Set functions
for i in range(random_size):
    obj += s[i]*o[i]
    m.update()

for i in range(x_random_size):   
    f_x += c[i]*x[i]
    m.update()

# for i in range(random_size):
#  #   g_x_si[i] = g_x + si[i] 
#     g_x = 0
#     for j in range(x_random_size):
#         g_x += a[i,j]*x[j]
#         m.update()
#     g_x_xi[i] = g_x 
for k in range(sample_size):
    g_x_xi = [None] * random_size
    for i in range(random_size):
         #   g_x_si[i] = g_x + si[i] 
        g_x = 0
            
        for j in range(x_random_size):
            g_x += a[k,i,j]*x[j]
            m.update()
        g_x_xi[i] = g_x
    g_x_k[k]=g_x_xi
        

# Add objective    
m.setObjective(obj, GRB.MINIMIZE)
# Add constraints

m.addConstrs(g_x_k[k][i] >= bb[i] - s[i]  for k in range(sample_size) for i in range(random_size))
con1=m.addConstr(f_x <= t)
# for i in range(random_size):
# s[i].start = CVAR_s[i]
# for i in range(x_random_size):
# x[i].start = CVAR_x[i] 
# Solve the problem
m.params.OutputFlag=0 
m.optimize()

# Store solutions
ii = m.getAttr('x', s)
kkk = m.getAttr('x', x)
iteration = 0


while delta_t >= delta_1: 
    
    support=0
    for j in ii:
        if ii[j]!=0:
            support +=1
    if support <= violation:
        t2=con1.SARHSLow
    else:
        t3=con1.SARHSUp
        t4=con1.SARHSLow
        
        AP_checking=AP(t,ii,kkk)
#         AP_checking=AP(t)
        if AP_checking<=violation:
            t2=t4
        else:
            t1=t3
    t=(t1+t2)/2
    delta_t=t2-t1
    con1.RHS = t
    for i in range(random_size):
        s[i].start = ii[i]
    for i in range(x_random_size):
        x[i].start = kkk[i]  
    m.params.OutputFlag=0 
    m.optimize()
    iteration+=1
    ii = m.getAttr('x', s)
    kkk = m.getAttr('x', x)
con1.RHS = t2
m.params.OutputFlag=0 
m.optimize()
ii = m.getAttr('x', s)
kkk = m.getAttr('x', x)
am_modeltime = time.time() - start
am_value=t2


# In[13]:


print('cvar')
# print(f'The modeltime is {modeltime:.3f} s.')
#print(f'The number of iterations is {iteration:.0f}.')
print(f'The modeltime of cvar is {CVAR_time:.3f} s.')
print(f'The objective s is {CVAR:.5f}.')


# In[14]:


print('relaxation')
# print(f'The modeltime is {modeltime:.3f} s.')
#print(f'The number of iterations is {iteration:.0f}.')
print(f'The modeltime of relaxation is {relaxation_runtime:.3f} s.')
print(f'The objective s is {value_rel:.5f}.')


# In[15]:


print('relax and scale')
print(f'The modeltime is {modeltime:.3f} s.')
print(f'The objective s is {value_scale:.5f}.')


# In[16]:


print('alsox')
print(f'The modeltime is {heuristics_modeltime:.3f} s.')
print(f'The objective s is {heuristics_value:.5f}.')


# In[17]:


print('alsox+')
print(f'The modeltime is {am_modeltime:.3f} s.')
print(f'The objective s is {am_value:.5f}.')


# In[ ]:




