import math
from functools import reduce

import numpy as np
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import sv_ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import plottingmath as pm

import sympy as sp
from sympy import symbols, Function, Eq, sympify, parse_expr
from sympy.parsing.sympy_parser import standard_transformations, implicit_multiplication_application


class MainWindow:
    def __init__(self):
        self.w = 720
        self.h = 800
        self.root = tk.Tk()
        self.root.resizable(False, False)
        self.root.title("TuringPlotsUI")
        self.root.protocol('WM_DELETE_WINDOW', self.main_window_was_closed)
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        self.root.geometry(str(self.w) + "x" + str(self.h) + "+"
                           + str(int((screen_width - self.w) / 2)) + "+" + str(int((screen_height - self.h) / 2)))

        self.systems = ["Брюсселятор", "Система Шнакенберга"]
        self.cb_systems = ttk.Combobox(values=self.systems, width=25, state="readonly")
        self.cb_systems.current(0)
        self.cb_systems.place(x=10, y=30)

        self.variables = ["(a, b)", "(DetJ, f_u)"]
        self.cb_variables = ttk.Combobox(values=self.variables, width=25, state="readonly")
        self.cb_variables.current(0)
        self.cb_variables.place(x=10, y=80)

        self.lbl_l = ttk.Label(text="Длина отрезка L")
        self.lbl_l.place(x=300, y=10)
        self.entry_l = ttk.Entry(width=25)
        self.entry_l.place(x=300, y=30)
        self.entry_l.insert(0, "0")

        self.lbl_coef_diff = ttk.Label(text="Коэффициент диффузии")
        self.lbl_coef_diff.place(x=300, y=60)
        self.entry_coef_diff = ttk.Entry(width=25)
        self.entry_coef_diff.place(x=300, y=80)
        self.entry_coef_diff.insert(0, "100")

        self.lbl_ax_thresholds = ttk.Label(text="Границы осей")
        self.lbl_ax_thresholds.place(x=300, y=110)

        self.lbl_x_min = ttk.Label(text="x min")
        self.lbl_x_min.place(x=300, y=130)
        self.lbl_x_max = ttk.Label(text="x max")
        self.lbl_x_max.place(x=350, y=130)

        self.lbl_y_min = ttk.Label(text="y min")
        self.lbl_y_min.place(x=400, y=130)
        self.lbl_y_max = ttk.Label(text="y max")
        self.lbl_y_max.place(x=450, y=130)

        self.entry_x_min = ttk.Entry(width=5)
        self.entry_x_min.place(x=300, y=150)
        self.entry_x_max = ttk.Entry(width=5)
        self.entry_x_max.place(x=350, y=150)
        self.entry_y_min = ttk.Entry(width=5)
        self.entry_y_min.place(x=400, y=150)
        self.entry_y_max = ttk.Entry(width=5)
        self.entry_y_max.place(x=450, y=150)

        self.entry_x_min.insert(0, "0")
        self.entry_x_max.insert(0, "5")
        self.entry_y_min.insert(0, "0")
        self.entry_y_max.insert(0, "5")

        self.btn_confirm = ttk.Button(text="Применить", command=self.confirm_plotting)
        self.btn_confirm.place(x=10, y=200)

        self.btn_save_plot = ttk.Button(text="Сохранить график", command=self.save_plot)
        self.btn_save_plot.place(x=300, y=200)

        self.fig, self.ax = plt.subplots(figsize=(7, 5))
        self.frame_plot = tk.Frame(self.root)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas.get_tk_widget().pack()
        # self.canvas.draw()

        self.frame_plot.place(x=10, y=230)

    def main_window_was_closed(self):
        plt.close(self.fig)
        self.root.destroy()

    def confirm_plotting(self):
        self.ax.cla()

        sys = self.cb_systems.get()
        x_lim = [float(self.entry_x_min.get()), float(self.entry_x_max.get())]
        y_lim = [float(self.entry_y_min.get()), float(self.entry_y_max.get())]
        d_var = float(self.entry_coef_diff.get())
        l_var = float(self.entry_l.get())

        ax_labels = self.cb_variables.get()
        vars_list = [ax_labels.split(', ')[0][1:], ax_labels.split(', ')[1][:-1]]

        u, v, a, b = sp.symbols('u, v, a, b')
        detj, f_u = sp.symbols('DetJ, f_u')
        d = sp.symbols('d')
        mu = sp.symbols('mu')

        f = None
        g = None
        if sys == "Брюсселятор":
            f = a - (b + 1) * u + u ** 2 * v
            g = b * u - u ** 2 * v
        elif sys == "Система Шнакенберга":
            f = u ** 2 * v - u + a
            g = -u ** 2 * v + b

        jac = pm.jacobian(f, g, u, v)

        ox_lnsp = np.linspace(x_lim[0], x_lim[1], 500)
        oy_lnsp = np.linspace(y_lim[0], y_lim[1], 500)
        ox_grid, oy_grid = np.meshgrid(ox_lnsp, oy_lnsp)

        necessary_inequals = [cond(ox_grid, oy_grid, d_var) for cond in pm.necessary_conds(jac, vars_list, d, a, b)]
        region = np.ones_like(ox_grid)
        for ineq in necessary_inequals:
            region[ineq == False] = 0
        necessary_region = region

        self.ax.contourf(ox_grid, oy_grid, necessary_region, hatches=['', '\\\\'], alpha=0, levels=[-0.5, 0.5, 1.5])
        styles_list = ['solid', 'dashed', 'dashdot']
        for i, ineq in enumerate(necessary_inequals):
            self.ax.contour(ox_grid, oy_grid, ineq, linestyles=styles_list[i], linewidths=1, colors='k', levels=[0])

        self.ax.set(xlabel=vars_list[0], ylabel=vars_list[1])
        self.ax.grid(False)
        self.canvas.draw()

    def save_plot(self):
        types = [("PNG", "*.png")]
        file = filedialog.asksaveasfile(title="Save File",
                                        initialdir=".",
                                        defaultextension=".png",
                                        filetypes=types)
        if file:
            fname = file.name
            file.close()
            self.fig.savefig(fname)


if __name__ == "__main__":
    main_window = MainWindow()
    sv_ttk.set_theme("light")
    main_window.root.mainloop()
