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
random_size=400
x_random_size=100


# In[3]:


# set seed
np.random.seed(35)


# In[4]:


# distribution
a = np.random.randint(1,50, size=(random_size,x_random_size))
c = np.random.randint(-10,-1,x_random_size) # matrix coefficient for f(x)
bb = np.random.normal(100,0,random_size) 


# In[5]:


# Set initial values for epsilon and theta
epsilon = 0.05
delta_r = 0.01
violation=math.floor(random_size*epsilon)


# In[6]:


# CVaR
def Model_cvar():

    m = Model()
    #m.setParam('TimeLimit', 60*60)
    # Create variables
    m.update()

    x = m.addVars(x_random_size,lb=0,ub=1,name="x")

    s = m.addVars(random_size,lb=-GRB.INFINITY,name="s")

    gamma = m.addVar(lb=-GRB.INFINITY,ub=0,name="gamma")
    lambda_1 = m.addVar(lb=0,name="lambda_1")
    m.update()
    

    # Add objective    
    m.setObjective(quicksum(c[i]*x[i] for i in range(x_random_size)), GRB.MINIMIZE)
    m.update()
    
    # Add constraints
   
    m.addConstrs(quicksum(a[i,j]*x[j] for j in range(x_random_size))  <= bb[i] + s[i]  for i in range(random_size))
 
    m.addConstrs(s[j] >= gamma for j in range(random_size))

    m.addConstr(lambda_1*delta_r + s.sum()/float(random_size)-(1.0-epsilon)*gamma <= 0)
     
    m.addConstr(sum(x[j]*x[j] for j in range(x_random_size))<=lambda_1*lambda_1)
    
    # Solve the problem
    m.update()
    m.params.OutputFlag=0
        
    m.optimize()
    

#     print('Obj: %g' % m.objVal)
#     for v in m.getVars():
#         print('%s %g' % (v.varName, v.x))
#     return m.objVal
    return m.objVal,m.Runtime
    


# In[7]:


Model_cvar()


# In[8]:


# DRCC set
def Model_checking(V_x):
    F= [0] * random_size
    for i in range(random_size):
        F[i] = -(max((bb[i]-sum( a[i,j]*V_x.values()[j] for j in range(x_random_size))),0))/math.sqrt(sum(V_x.values()[j]*V_x.values()[j] for j in range(x_random_size))) # g(x) is \xi*x pointwise
     
        
    if delta_r == 0:

        s_value=[]
        for i in range(random_size):
            s_value.append( V_s[i])
        sorted_array=np.sort(s_value)
        result = sorted_array[j_size-violation-1]

        
    else:
        if (epsilon%(1/random_size))<=1e-10:
            sorted_array=np.sort(F)
            truncated = sorted_array[random_size-violation:]
            result = np.average(truncated) + delta_r / epsilon
        else:
            integer_part = epsilon//(1/random_size)
            fraction_part = epsilon%(1/random_size)
            sorted_array=np.sort(F)
            truncated = sorted_array[random_size-violation:]
            result = delta_r / epsilon + (fraction_part * sorted_array[random_size-violation-1] + sum(truncated)*(1/random_size))/epsilon
    
    return result

  


# In[9]:


# alsox
def alsox():

    start=time.time()
    # Set initial values
    delta_1 = 1e-2
    t1 = -100.0
    t2 = -0
    delta_t = t2-t1
    t = (t1+t2)/2.0
    m = Model()

    # Create variables
    m.update()
    x = m.addVars(x_random_size,lb=0,ub=1,name="x")
    s = m.addVars(random_size,lb=0,name="s")
    lambda_1 = m.addVars(1,lb=0,name="lambda_1")

    m.update()

    m.setObjective(lambda_1[0]*delta_r + s.sum()/float(random_size), GRB.MINIMIZE)
    m.update()

        # Add constraints

    m.addConstrs(quicksum(a[i,j]*x[j] for j in range(x_random_size))  <= bb[i] + s[i]  for i in range(random_size))
 
    m.addConstrs(s[j] >= 0 for j in range(random_size))



    m.update()

    con1=m.addConstr(quicksum(c[i]*x[i] for i in range(x_random_size)) <= t)

    m.addConstr(sum(x[j]*x[j] for j in range(x_random_size))<=lambda_1[0]*lambda_1[0])
  
        # Solve the problem
    m.update()
    m.params.OutputFlag=0

    m.optimize()

    # Store solutions
    V_s = m.getAttr('x', s)
    V_x = m.getAttr('x', x)


    while delta_t >= delta_1: 
        checking = Model_checking(V_x)

        if checking <= 0:
#             t2=con1.SARHSLow
            t2=t
#             print ("t2 is ", t2)
        else:
    #         t1=con1.SARHSUp
            t1=t
#             print ("t1 is ", t1)
        t=(t1+t2)/2.0
        delta_t=t2-t1
        con1.RHS = t
        for i in range(random_size):
            s[i].start = V_s[i]
        for i in range(x_random_size):
            x[i].start = V_x[i]  
        m.optimize()
        V_s = m.getAttr('x', s)
        V_x = m.getAttr('x', x)
        lambda_ = m.getAttr('x', lambda_1)


    modeltime_orginial = time.time() - start
    solution =t2
    return t2,modeltime_orginial


# In[10]:


# alsox()


# In[ ]:





# In[11]:


def Proj(c,t,x_hat):

    m = Model()

    # Create variables
    m.update()

    x = m.addVars(x_random_size,lb=0,ub=1,name="x")

    m.update()
    
    m.addConstr(quicksum( c[i]*x[i] for i in range(x_random_size)) <= t)
    
    m.setObjective(quicksum( (x[i]-x_hat[i][0])*(x[i]-x_hat[i][0]) for i in range(x_random_size)), GRB.MINIMIZE)
    m.params.OutputFlag=0 

    m.optimize()
    
    x_val = m.getAttr('x', x)
    return x_val


# In[12]:


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


# In[13]:


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


# In[14]:


def g(x,ep):
    return np.dot(ep,x)
    
def dg(x,ep):
    m=np.dot(x,np.ones((1,ep.shape[0]))).transpose()
    return ep


# In[15]:


def Model_checking_sd(x_s):
    F= [0] * random_size
    for i in range(random_size):
        F[i] = -(max((bb[i]-sum( a[i,j]*x_s[j] for j in range(x_random_size))),0))/math.sqrt(sum(x_s[j]*x_s[j] for j in range(x_random_size))) # g(x) is \xi*x pointwise
     
        
    if delta_r == 0:
        V_s=[0] * random_size
        s_value=[]
        for i in range(random_size):
            s_value.append( V_s[i])
        sorted_array=np.sort(s_value)
        result = sorted_array[j_size-violation-1]

        
    else:
        if (epsilon%(1/random_size))<=1e-10:
            sorted_array=np.sort(F)
            truncated = sorted_array[random_size-violation:]
            truncated= np.array(truncated,dtype=np.float32)
            result = np.average(truncated) + delta_r / epsilon
        else:
            integer_part = epsilon//(1/random_size)
            fraction_part = epsilon%(1/random_size)
            sorted_array=np.sort(F)
            truncated = sorted_array[random_size-violation:]
            truncated= np.array(truncated,dtype=np.float32)
            result = delta_r / epsilon + (fraction_part * sorted_array[random_size-violation-1] + sum(truncated)*(1/random_size))/epsilon
    
    return result

  


# In[16]:


# alsox SD
start=time.time()

support=random_size
delta_1= 1e-2
t_1= -30
t_2= -40
t=(t_1+t_2)/float(2)
x_s=np.zeros((x_random_size,1))
x_s = abs(t)/np.sum(abs(c)) * np.ones((x_random_size,1))


while t_1-t_2>=delta_1:
#     t_last = t
    t=(t_1+t_2)/float(2)
#     x_s = x_s/t_last*t
    x_s = proj_d(x_s,x_random_size,c,t)

    # summation_1=0
    # for i in range(x_random_size):

    #     if summation_1+abs(c)[i]*1<=abs(t):
    #         x_s[i]=1
    #         summation_1=summation_1+abs(c)[i]
    #     else:
    #         x_s[i]=(abs(t)-summation_1)/float(abs(c)[i])
    #         summation_1=t
    #         break
    #     x_s[i] = x_value.values()[i]


    subgrad = np.zeros((random_size,x_random_size))
    obj=100
    iteration = 0 

    Difference=1e20

    maxiter=200


    while iteration<=maxiter:
#     while iteration<=maxiter and Difference>1e-5:
        
        iteration +=1
        g_val=np.dot(a,x_s)
        
        extended_g_val=np.dot(g_val,np.ones((1,a.shape[1])))
        subgrad=copy.deepcopy(a)

        subgrad[extended_g_val<=100+1e-5]=0

#         subgrad_coefficient = np.average(subgrad,axis=0)
        subgrad_coefficient = (x_s/math.sqrt(np.dot(x_s.transpose(),x_s))*delta_r).transpose()+np.average(subgrad,axis=0)

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
    

#     g_val=g(x_f,a)
#     support=np.count_nonzero(g_val>100+1e-5)

    checking = Model_checking_sd(x_s)
    
#     print('current checking is : ', checking)

    if checking <= 0:
        t_1=t
    else:
        t_2=t
    t=(t_1+t_2)/2.0   
#     print('current t is : ', t)
#     print('current a is: ', a)


modeltime = time.time() - start

alsox_value = t_1


# In[17]:


print('G-D')
print(f'The modeltime is {modeltime:.3f} s.')
print(f'The objective s is {alsox_value:.5f}.')


# In[ ]:





# In[18]:


CVAR = Model_cvar()


# In[19]:


start=time.time()
coefficient=delta_r/epsilon 
kkk=0.9

while alsox_value >= 1.001*CVAR[0]:
#     k1=(1+x/100)*coefficient
    kkk+=0.1
    rho=kkk*coefficient
    print (kkk, " of kkk")
    print (alsox_value, " solution")
    # Set initial values


    delta_1= 1e-2
    t_1= -30
    t_2= -40
    t=(t_1+t_2)/float(2)
    x_s=np.zeros((x_random_size,1))
    x_s = abs(t)/np.sum(abs(c)) * np.ones((x_random_size,1))


    while t_1-t_2>=delta_1:
    #     t_last = t
        t=(t_1+t_2)/float(2)
    #     x_s = x_s/t_last*t
        x_s = proj_d(x_s,x_random_size,c,t)

        # summation_1=0
        # for i in range(x_random_size):

        #     if summation_1+abs(c)[i]*1<=abs(t):
        #         x_s[i]=1
        #         summation_1=summation_1+abs(c)[i]
        #     else:
        #         x_s[i]=(abs(t)-summation_1)/float(abs(c)[i])
        #         summation_1=t
        #         break
        #     x_s[i] = x_value.values()[i]


        subgrad = np.zeros((random_size,x_random_size))
        obj=100
        iteration = 0 

        Difference=1e20

        maxiter=200


        while iteration<=maxiter:
    #     while iteration<=maxiter and Difference>1e-5:

            iteration +=1
            g_val=np.dot(a,x_s) + math.sqrt(np.dot(x_s.transpose(),x_s))*rho

            extended_g_val=np.dot(g_val,np.ones((1,a.shape[1])))
            
            subgrad=copy.deepcopy(a) + (x_s/math.sqrt(np.dot(x_s.transpose(),x_s))*rho).transpose()

            subgrad[extended_g_val<=100+1e-5]=0

    #         subgrad_coefficient = np.average(subgrad,axis=0)
            subgrad_coefficient = (x_s/math.sqrt(np.dot(x_s.transpose(),x_s))*delta_r).transpose()+np.average(subgrad,axis=0)

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


    #     g_val=g(x_f,a)
    #     support=np.count_nonzero(g_val>100+1e-5)

        checking = Model_checking_sd(x_s)

    #     print('current checking is : ', checking)

        if checking <= 0:
            t_1=t
        else:
            t_2=t
        t=(t_1+t_2)/2.0   
    #     print('current t is : ', t)
    #     print('current a is: ', a)
    
    
    alsox_value = t_2
    
    


modeltime_training = time.time() - start
    
    
    


# In[20]:


print('Training')
print(f'The modeltime is {modeltime_training:.3f} s.')
print(f'The  value of rho is {rho:.5f}.')


# In[21]:


# Testing
start=time.time()

delta_1= 1e-2
t_1= -30
t_2= -40
t=(t_1+t_2)/float(2)
x_s=np.zeros((x_random_size,1))
x_s = abs(t)/np.sum(abs(c)) * np.ones((x_random_size,1))


while t_1-t_2>=delta_1:
#     t_last = t
    t=(t_1+t_2)/float(2)
#     x_s = x_s/t_last*t
    x_s = proj_d(x_s,x_random_size,c,t)

    # summation_1=0
    # for i in range(x_random_size):

    #     if summation_1+abs(c)[i]*1<=abs(t):
    #         x_s[i]=1
    #         summation_1=summation_1+abs(c)[i]
    #     else:
    #         x_s[i]=(abs(t)-summation_1)/float(abs(c)[i])
    #         summation_1=t
    #         break
    #     x_s[i] = x_value.values()[i]


    subgrad = np.zeros((random_size,x_random_size))
    obj=100
    iteration = 0 

    Difference=1e20

    maxiter=200


    while iteration<=maxiter:
#     while iteration<=maxiter and Difference>1e-5:

        iteration +=1
        g_val=np.dot(a,x_s) + math.sqrt(np.dot(x_s.transpose(),x_s))*rho

        extended_g_val=np.dot(g_val,np.ones((1,a.shape[1])))

        subgrad=copy.deepcopy(a) + (x_s/math.sqrt(np.dot(x_s.transpose(),x_s))*rho).transpose()

        subgrad[extended_g_val<=100+1e-5]=0

#         subgrad_coefficient = np.average(subgrad,axis=0)
        subgrad_coefficient = (x_s/math.sqrt(np.dot(x_s.transpose(),x_s))*delta_r).transpose()+np.average(subgrad,axis=0)

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


#     g_val=g(x_f,a)
#     support=np.count_nonzero(g_val>100+1e-5)

    checking = Model_checking_sd(x_s)

#     print('current checking is : ', checking)

    if checking <= 0:
        t_1=t
    else:
        t_2=t
    t=(t_1+t_2)/2.0   
#     print('current t is : ', t)
#     print('current a is: ', a)


testing_time = time.time() - start
    
training_value = t_1    
    


# In[22]:


print('Testing')
print(f'The modeltime is {testing_time:.3f} s.')
print(f'The  value  is {training_value:.5f}.')


# In[ ]:




