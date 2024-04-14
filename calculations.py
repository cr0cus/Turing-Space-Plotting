import math
from functools import reduce

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

class symbols():
    def __init__():
        pass
        

class TuringRegionVisualizer():
    u, v, a, b = sp.symbols('u, v, a, b')
    detj, f_u = sp.symbols('DetJ, f_u')
    d = sp.symbols('d')
    mu = sp.symbols('mu')
    def __init__(self, sys, x_lim, y_lim):
        #self.sys = sys
        self.f = None
        self.g = None
        if sys == "Брюсселятор":
            self.f = type(self).a - (type(self).b + 1) * type(self).u + type(self).u ** 2 * type(self).v
            self.g = type(self).b * type(self).u - type(self).u ** 2 * type(self).v
        elif sys == "Система Шнакенберга":
            self.f = type(self).u ** 2 * type(self).v - type(self).u + type(self).a
            self.g = -type(self).u ** 2 * type(self).v + type(self).b
        self.x_lim = x_lim
        self.y_lim = y_lim

    # Функция возвращает точку равновесия системы в виде кортежа (u0, v0)
    def equilibrium_point(self):
        # метод solve возвращает множества решений в виде списка кортежей,
        # берём первый элемент (пока предполагается, что равновесие единственно)
        return sp.solve([self.f, self.g], [type(self).u, type(self).v])[0]

    # Функция возвращает список элементов матрицы Якоби на равновесии (u0, v0)
    def jacobian(self):
        fu = sp.diff(self.f, type(self).u)
        fv = sp.diff(self.f, type(self).v)
        gu = sp.diff(self.g, type(self).u)
        gv = sp.diff(self.g, type(self).v)
    
        diff_list = [fu, fv, gu, gv]
        u0, v0 = self.equilibrium_point()
        
        return [func.subs({type(self).u: u0, type(self).v: v0}) for func in diff_list]

    def get_necessary_conds(self, vars_list):
        fu, fv, gu, gv = self.jacobian()
        det_j = sp.simplify(fu*gv - fv*gu)
        discr = (type(self).d*fu + gv) - 2*sp.sqrt(type(self).d)*sp.sqrt(det_j)
    
        trace_cond = (fu + gv) < 0
        det_j_cond = det_j > 0
        discr_h_cond = discr >= 0
    
        if vars_list == ['a', 'b']:
            trace_cond = sp.lambdify([type(self).a, type(self).b, type(self).d], trace_cond, 'numpy') 
            det_j_cond =  sp.lambdify([type(self).a, type(self).b, type(self).d],  det_j_cond, 'numpy') 
            discr_h_cond = sp.lambdify([type(self).a, type(self).b, type(self).d], discr_h_cond, 'numpy')
        elif vars_list == ['DetJ', 'f_u']:    
            trace_cond = trace_cond.subs({det_j: type(self).detj, fu: type(self).f_u})
            det_j_cond = det_j_cond.subs({det_j : type(self).detj, fu: type(self).f_u})
            discr_h_cond = discr_h_cond.subs({det_j : detj, fu: type(self).f_u}) 
    
            trace_cond = sp.lambdify([type(self).detj, type(self).f_u, type(self).d], trace_cond, 'numpy') 
            det_j_cond =  sp.lambdify([type(self).detj, type(self).f_u, type(self).d],  det_j_cond, 'numpy') 
            discr_h_cond = sp.lambdify([type(self).detj, type(self).f_u, type(self).d], discr_h_cond, 'numpy')
        
        return [trace_cond, det_j_cond, discr_h_cond] 
        
    def get_necessary_region(self, vars_list, ox_grid, oy_grid, d_var):
        necessary_conds = [cond(ox_grid, oy_grid, d_var) for cond in self.get_necessary_conds(vars_list)]
        region = np.ones_like(ox_grid)
        for ineq in necessary_conds:
            region[ineq == False] = 0
        return region, necessary_conds

    # ! пока не используется, но позже будет прикручено к plot_turing_space для достаточных условий
    # Вычисляет k-тое собственное значение оператора лапласа на интервале (0, l)
    def get_mu_k(self, k, l):                  
        return math.pi * k / l

    # Возвращает многочлен h(mu)
    def h(self):
        fu, fv, gu, gv = self.jacobian()
        det_j = sp.simplify(fu*gv - fv*gu)
        
        return d*mu**2 + sp.simplify((d*fu + gv))*mu + det_j
        
    # Возвращает многочлен h в точке mu_k
    def h_mu_k(self, mu_k):
        return self.h().subs({mu: mu_k})

    def suff_curves_intersection_point(self):

        pass

    def get_sufficient_region(self, vars_list, ox_grid, oy_grid, d_var, l_var):
        '''
         # И достаточные:
        suff_cond_lines = []
        for k in range(len(mu)):
            if k==0:
                a_min = 0
                a_max = math.sqrt(d*mu[k]*mu[k+1])
            elif k<len(mu)-1:
                a_min = math.sqrt(d*mu[k-1]*mu[k])
                a_max = math.sqrt(d*mu[k]*mu[k+1])
            else:
                a_min = math.sqrt(d*mu[k-1]*mu[k])
                a_max = max(np.max(A), a_min)
        suff_cond = (B >= mu[k] + A**2/(d*mu[k]) + A**2/d + 1) & (A>=a_min) & (A<a_max)
        suff_cond_lines.append(suff_cond)
    
        sufficient_region = reduce(lambda x, y: x | y, suff_cond_lines)
    
        sufficient_region[ineq1 == False] = False
        sufficient_region[ineq2 == False] = False
        '''
        pass
           
    


