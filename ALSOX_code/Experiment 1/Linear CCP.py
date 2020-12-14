#!/usr/bin/env python
# coding: utf-8

# In[8]:


import sys
from gurobipy import *
import numpy as np
import math
import time
from numpy import linalg as LA


# In[9]:


# Set parameters for the problem
# random_size = int(input('Input the size of the model: '))
# x_random_size = int(input('Input the size of the model of x: '))
# sample_size = int(input('Input the constraint size of the model I: '))

random_size = 1000
x_random_size  = 40
sample_size = 1


# In[10]:


# set seed
np.random.seed(45)


# In[11]:


# distribution
a = np.random.randint(1,50, size=(sample_size,random_size,x_random_size,))
c = np.random.randint(-10,-1,x_random_size)
bb = np.random.normal(100,0,random_size)


# In[12]:


# Set initial values for epsilon
epsilon = 0.05
violation=math.floor(random_size*epsilon)


# In[13]:


# CVaR Method

def Model_cvar():
    #Model
    print ("Begin to solve model 4")

    m = Model()
    #m.setParam('TimeLimit', 60*60)
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

    for k in range(sample_size):
        g_x_xi = [None] * random_size
        for i in range(random_size):
            g_x = 0
            
            for j in range(x_random_size):
                g_x += a[k,i,j]*x[j]
                m.update()
            g_x_xi[i] = g_x
        g_x_k[k]=g_x_xi
    # not necessary
    
    # Add objective    
    m.setObjective(quicksum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    # Add constraints
   
    m.addConstrs(g_x_k[j][i] <= bb[i] + s[i]  for j in range(sample_size) for i in range(random_size))
    m.addConstrs(s[i] >= beta for i in range(random_size))
    
    m.addConstr(s.sum()/random_size -beta*(1-epsilon) <= 0)
  
    # Solve the problem
    m.update()
        
    m.optimize()
    

#     Store solutions
    CVAR_s = m.getAttr('x', s)
    CVAR_x = m.getAttr('x', x)
    
    return m.objVal,CVAR_s,CVAR_x
    
    


# In[14]:


CVaR=Model_cvar()
CVAR=CVaR[0]
CVAR_s=CVaR[1]
CVAR_x=CVaR[2]


# In[15]:


# quantitle bounds

def Model_q():

    print ("Begin to solve model quantitle relaxtion")

    m = Model()
    
    m.setParam('TimeLimit', 60*60)
    m.update()
    
    # Create variables
    x = m.addVars(x_random_size,lb=0,ub=1, name="x")
    m.update()
    
    # Set functions
    
    #obj = 0
    g_x_k = [None] * sample_size
    f_x = 0
    g_x = 0
    
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
    m.setObjective(quicksum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    # Add constraints
    

    con1=m.addConstrs(g_x_k[j][0] <= bb[0]  for j in range(sample_size))
    # Solve the problem
    m.update()
    m.params.OutputFlag=0    
    m.optimize()
    kkk = m.getAttr('x', x)
    result = [0]*random_size
    obj = 0
    for i in range(x_random_size):
        obj += c[i]*kkk[i]
    result[0]=obj
    
    
    for i in range(1,random_size):   
        m.remove(con1)
        con1=m.addConstrs( g_x_k[j][i]  <= bb[i] for j in range(sample_size) )
        m.update()
        m.params.OutputFlag=0    
        m.optimize()
        kkk = m.getAttr('x', x)
        obj = 0
        for j in range(x_random_size):
            obj += c[j]*kkk[j]
        result[i]=obj
    return result
        


# In[16]:


result=Model_q()
sorted_array=np.sort(result)
v_q = sorted_array[random_size-violation+1]
print(v_q)


# In[17]:


# alsox
print ("Begin to solve alsox")

start=time.time()
# Set initial values
delta_1=1e-2
t1 =  v_q
t2 =  CVAR
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

m.addConstrs(g_x_k[k][i] <= bb[i] + s[i]  for k in range(sample_size) for i in range(random_size))
con1=m.addConstr(f_x <= t)
for i in range(random_size):
    s[i].start = CVAR_s[i]
for i in range(x_random_size):
    x[i].start = CVAR_x[i] 
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
modeltime = time.time() - start


# In[18]:


print('alsox')
print(f'The modeltime is {modeltime:.3f} s.')
print(f'The number of iterations is {iteration:.0f}.')
print(f'The objective s is {t2:.5f}.')


# In[19]:


def AM(t,ii,kkk):
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
    m.addConstrs(g_x_k[k][i] <= bb[i]+s[i] for i in range(random_size) for k in range(sample_size))
    
    for i in range(random_size):
        s[i].start = CVAR_s[i]
    for i in range(x_random_size):
        x[i].start = CVAR_x[i] 
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
 
    


# In[20]:


# alsox+
print ("Begin to solve alsox + ")



start=time.time()


delta_1 = 1e-3
t1 =  v_q
t2 =  CVAR
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

m.addConstrs(g_x_k[k][i] <= bb[i] + s[i] for i in range(random_size) for k in range(sample_size))
con1=m.addConstr(f_x <= t)
for i in range(random_size):
    s[i].start = CVAR_s[i]
for i in range(x_random_size):
    x[i].start = CVAR_x[i] 
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
        
        AP_checking=AM(t,ii,kkk)
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
modeltime = time.time() - start


# In[21]:


print('alsox+')
print(f'The modeltime is {modeltime:.3f}.')
print(f'The number of iterations is {iteration:.0f}.')
print(f'The objective s is {t2:.5f}.')


# In[22]:






# In[ ]:





# In[ ]:





# In[ ]:




