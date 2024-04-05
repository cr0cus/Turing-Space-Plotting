import math
from functools import reduce

import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt

import sympy as sp
from sympy import symbols, Function, Eq, sympify, parse_expr
from sympy.parsing.sympy_parser import standard_transformations, implicit_multiplication_application


# Функция возвращает точку равновесия системы в виде кортежа (u0, v0)
def equilibrium_point(f, g, u, v):
    # метод solve возвращает множества решений в виде списка кортежей,
    # берём первый элемент (пока предполагается, что равновесие единственно)
    return sp.solve([f, g], [u, v])[0]


# Функция возвращает список элементов матрицы Якоби на равновесии (u0, v0)
def jacobian(f, g, u, v):
    fu = sp.diff(f, u)
    fv = sp.diff(f, v)
    gu = sp.diff(g, u)
    gv = sp.diff(g, v)

    diff_list = [fu, fv, gu, gv]
    u0, v0 = equilibrium_point(f, g, u, v)

    return [func.subs({u: u0, v: v0}) for func in diff_list]


# ! пока не используется, но позже будет прикручено к plot_turing_space для достаточных условий
# Возвращает многочлен h(mu)
def h(jac, d, mu):
    fu, fv, gu, gv = jac
    det_j = sp.simplify(fu*gv - fv*gu)

    return d*mu**2 + sp.simplify((d*fu + gv))*mu + det_j

def necessary_conds(jac, vars_list, d, a, b):
    fu, fv, gu, gv = jac
    det_j = sp.simplify(fu*gv - fv*gu)
    discr = (d*fu + gv) - 2*sp.sqrt(d)*sp.sqrt(det_j)

    detj, f_u = sp.symbols('DetJ, f_u')

    trace_cond = (fu + gv) < 0
    det_j_cond = det_j > 0
    discr_h_cond = discr >= 0

    if vars_list == ['a', 'b']:
        trace_cond = sp.lambdify([a, b, d], trace_cond, 'numpy') 
        det_j_cond =  sp.lambdify([a, b, d],  det_j_cond, 'numpy') 
        discr_h_cond = sp.lambdify([a, b, d], discr_h_cond, 'numpy')
    elif vars_list == ['DetJ', 'f_u']:    
        trace_cond = trace_cond.subs({det_j: detj, fu: f_u})
        det_j_cond = det_j_cond.subs({det_j : detj, fu: f_u})
        discr_h_cond = discr_h_cond.subs({det_j : detj, fu: f_u}) 

        trace_cond = sp.lambdify([detj, f_u, d], trace_cond, 'numpy') 
        det_j_cond =  sp.lambdify([detj, f_u, d],  det_j_cond, 'numpy') 
        discr_h_cond = sp.lambdify([detj, f_u, d], discr_h_cond, 'numpy')
    
    return [trace_cond, det_j_cond, discr_h_cond] 


# ! пока не используется, но позже будет прикручено к plot_turing_space для достаточных условий
# Возвращает список из n собственных значений на интервале (0, l)
def compute_mu_list(n, l):
    mu = [(math.pi * k / l)**2 for k in range(1, n+1)]
    return mu


