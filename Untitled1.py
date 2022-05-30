#!/usr/bin/env python
# coding: utf-8

# In[10]:


from sympy.interactive import printing
printing.init_printing(use_latex=True)
from sympy import *
import sympy as sp
import math
import numpy as np
from scipy.integrate import odeint

r = sp.Symbol('r',real=True)
#r = float(input('Enter r variable :'))

r1 = 0.5

r2= 0.8

n = 3



f_r = (( r - r1 )/( r2 - r1 )) ** n

print('f(r) is : ',f_r)

print('----------------------------------------')

v = 0.3
E_0 = 210
E_a = 70
E_m = 30
E_avg = (E_0 + E_a + E_m) / 3
landa_M = (E_avg * math.pow(10,9)) / ((1 + v) * (1 - (2 * v)))
landa_0 = (E_avg* math.pow(10,9)) / ((1 + v) * (1 - (2 * v)))
landa_A = (E_avg* math.pow(10,9)) / ((1 + v) * (1 - (2 * v)))

print(f'Landa(0) is : {landa_0}')
print(f'Landa(M) is : {landa_M}')
print(f'Landa(A) is : {landa_A}')

print('----------------------------------------')

g_0 = (E_avg * math.pow(10,9)) / (2 * (1 + v))
g_A = (E_avg* math.pow(10,9)) / (2 * (1 + v))
g_M = (E_avg* math.pow(10,9)) / (2 * (1 + v))

print(f'G(0) is : { g_0}')
print(f'G(A) is : { g_A}')
print(f'G(M) is : { g_M}')

print('----------------------------------------')

zay = 0
alpha_0 = 0.000004
alpha_A = 0.000011
alpha_M = 0.0000066
epsilon_T = 0.069

# T_r_1 = float(input("Enter T(r 1) Varibale : "))

# T_r_2 = float(input("Enter T(r 2) Varibale : "))

delta_T_r = 0

print('Delta T(r) is :',delta_T_r)


U = sp.Function('U')

epsilon_r_z = 0

#-----epsilon_r_r-----
eq1 = U(r).diff(r)

#-----epsilon_teta_r
eq2 = U(r)/r

#-----sigma_r_r_0
sec1 = (landa_0 + 2*g_0)*eq1 + landa_0*( eq2 + epsilon_r_z) - (3*landa_0 + 2*g_0)*alpha_0*delta_T_r

#-----sigma_teta_r_0
sec2 = (landa_0 + 2*g_0)*eq2 + landa_0*( eq1 + epsilon_r_z) - (3*landa_0 + 2*g_0)*alpha_0*delta_T_r

#----sig_z_r_0
sec3 = (landa_0 + 2*g_0)*epsilon_r_z + landa_0*( eq2 + eq1) - (3*landa_0 + 2*g_0)*alpha_0*delta_T_r

#----sigma_r_r_1
sec4 = (1-zay)*((landa_A + 2*g_A)*eq1 + landa_A*(eq2 + epsilon_r_z) - (3*landa_A + 2*g_A)*alpha_A*delta_T_r) + zay*((landa_M + 2*g_M)*eq1 + landa_M*(eq2 + epsilon_r_z) - (3*landa_M + 2*g_M)*(alpha_M*delta_T_r + epsilon_T))

#----sigma_teta_r_1
sec5 = (1-zay)*((landa_A + 2*g_A)*eq2 + landa_A*(eq1 + epsilon_r_z) - (3*landa_A + 2*g_A)*alpha_A*delta_T_r) + zay*((landa_M + 2*g_M)*eq2 + landa_M*(eq1 + epsilon_r_z) - (3*landa_M + 2*g_M)*(alpha_M*delta_T_r + epsilon_T))

#----sigma_z_r_1
sec6 = (1-zay)*((landa_A + 2*g_A)* epsilon_r_z + landa_A*(eq2 + eq1) - (3*landa_A + 2*g_A)*alpha_A*delta_T_r) + zay*((landa_M + 2*g_M)* epsilon_r_z + landa_M*(eq2 + eq1) - (3*landa_M + 2*g_M)*(alpha_M*delta_T_r + epsilon_T))

derivative = sp.Derivative(f_r,r)

#F = f_r

#F =  sp.Function('F')(r)

#dif_f_r = sp.diff(r)

#second_diff_f_r = sp.diff(dif_f_r)
#-----sigma_r_r

S_r = (1-f_r)*sec1 + f_r*sec4

#-----sigma_teta_r
S_t = (1-f_r)*sec2 + f_r*sec5


#display(derivative)

equilibrium = sp.Derivative(S_r,r) + ((S_r - S_t)/r)
equilibrium_1 = Eq(equilibrium,0)
display(equilibrium_1)
# dsolve(equilibrium_1)


a0= sp.Symbol('a0', real=True)
b0 = sp.Symbol('b0', real=True)

# r = sp.Symbol('r',real=True)

# u = sp.Function('U')

# Eq1 = sp.Eq(sp.diff(u(r),r,2)+(1/r*(sp.diff(u(r),r,1))) - ((1/r**2) * u(r)))
Eq1 = equilibrium_1

display(sp.simplify(Eq1))

u_sl0 = sp.dsolve(Eq1,U(r)).rhs
display(sp.Eq(U(r),u_sl0))

cond0 = sp.Eq(u_sl0.subs(r,0),a0)
cond1 = sp.Eq(u_sl0.subs(r,0),b0)
C1 = sp.Symbol("C1")
C2 = sp.Symbol("C2")
C1C2_sl = sp.solve([cond0,cond1],(C1,C2))
u_sl1 = u_sl0.subs(C1C2_sl)
display(sp.Eq(U(r),u_sl1))


# In[ ]:




