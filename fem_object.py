#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#                       Описание объекта расчета
###################################################################

import os
import sys
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
from math import fabs
from fem_mesh import TMesh
from fem_fem import TFEM
from fem_params import TFEMParams
from fem_static import TFEMStatic
from fem_dynamic import TFEMDynamic
from fem_defs import eps


# Вывод сообщения об ошибке
def error(err_msg):
    print('\033[1;31m%s\033[1;m' % err_msg)


class TObject:
    def __init__(self):
        self.__params__ = TFEMParams()  # Параметры расчета
        self.__mesh__ = TMesh()         # КЭ-модель
        self.__results__ = []           # Список результатов расчета для перемещений, деформаций, ...

    def set_mesh(self, name):
        self.__mesh__.load(name)
        print('Object: %s' % self.object_name())
        print('Points: %d' % len(self.__mesh__.x))
        print('FE: %d - %s' % (len(self.__mesh__.fe), self.__mesh__.fe_name()))
#        print('FE type: %s' % self.__mesh__.fe_type)

    # Название объекта
    def object_name(self):
        return os.path.splitext(os.path.basename(self.__mesh__.mesh_file))[0]

    def set_problem_type(self, problem_type):
        self.__params__.problem_type = problem_type

    def set_solve_method(self, solve_method):
        self.__params__.solve_method = solve_method

    def set_eps(self, e):
        self.__params__.eps = e

    def set_width(self, width):
        self.__params__.width = width

    def set_precision(self, precision):
        self.__params__.precision = precision

    def set_elasticity(self, e, m):
        self.__params__.e = e
        self.__params__.m = m

    def set_density(self, density):
        self.__params__.density = density

    def set_time(self, t0, t1, th):
        self.__params__.t0 = t0
        self.__params__.t1 = t1
        self.__params__.th = th

    def set_damping(self, damping):
        self.__params__.damping = damping

    def set_names(self, names):
        self.__params__.names = names

    def add_boundary_condition(self, e, p, d):
        self.__params__.add_boundary_condition(e, p, d)

    def add_initial_condition(self, e, d):
        self.__params__.add_initial_condition(e, '', d)

    def add_volume_load(self, e, p, d):
        self.__params__.add_volume_load(e, p, d)

    def add_surface_load(self, e, p, d):
        self.__params__.add_surface_load(e, p, d)

    def add_concentrated_load(self, e, p, d):
        self.__params__.add_concentrated_load(e, p, d)

    def add_variable(self, var, val):
        self.__params__.add_variable(var, val)

    def calc(self):
        fem = TFEM()
        if self.__params__.problem_type == 'static':
            fem = TFEMStatic()
        elif self.__params__.problem_type == 'dynamic':
            fem = TFEMDynamic()
        fem.set_mesh(self.__mesh__)
        fem.set_params(self.__params__)
        ret = fem.calc()
        if ret:
            self.__results__ = fem.get_result()
        return ret

    # Вывод результатов расчета
    def print_result(self, *argv):
        file = sys.stdout
        try:
            if len(argv) == 1:
                file = open(argv[0], 'w')
        except IOError:
            error('Error: unable to open file %s' % argv[0])
            return
        if self.__params__.problem_type == 'static':
            self.__print__(file)
        else:
            t = self.__params__.t0
            while t <= self.__params__.t1:
                file.write('t = %5.2f\n' % t)
                self.__print__(file, t)
                t += self.__params__.th
                if fabs(t - self.__params__.t1) < self.__params__.eps:
                    t = self.__params__.t1

        file.close()

    # Вывод результатов расчета для одного момента времени
    def __print__(self, file, t=0):
        # Определение ширины позиции
        len1 = len('%+*.*E' % (self.__params__.width, self.__params__.precision, 3.14159))
        len2 = len('%d' % len(self.__mesh__.x))
        # Вывод заголовка
        file.write('| %*s  (' % (len2, 'N'))
        for i in range(0, self.__mesh__.freedom):
            file.write(' %*s' % (len1, self.__params__.names[i]))
            if i < self.__mesh__.freedom - 1:
                file.write(',')
        file.write(') |')
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %*s |' % (len1, self.__results__[i].name))
        file.write('\n')
        for i in range(0, len(self.__mesh__.x)):
            file.write('| %*d  (' % (len2, i + 1))
            file.write(' %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.x[i]))
            if len(self.__mesh__.y):
                file.write(', %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.y[i]))
            if len(self.__mesh__.z):
                file.write(', %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.z[i]))
            file.write(') | ')
            for k in range(0, len(self.__results__)):
                if self.__results__[k].t == t:
                    file.write('%+*.*E' %
                               (self.__params__.width, self.__params__.precision, self.__results__[k].results[i]))
                    file.write(' | ')
            file.write('\n')
        file.write('\n')
        # Печать итогов
        file.write('|  %*s  ' % (len2, ' '))
        for i in range(0, self.__mesh__.freedom):
            file.write(' %*s' % (len1, ' '))
            if i < self.__mesh__.freedom - 1:
                file.write(' ')
        file.write('  |')
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %*s |' % (len1, self.__results__[i].name))
        file.write('\n')
        file.write('|   %*s  |' % (self.__mesh__.freedom*(len1 + 1) + self.__mesh__.freedom + len2, 'min:'))
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %+*.*E ' % (self.__params__.width, self.__params__.precision, self.__results__[i].min()))
                file.write('|')
        file.write('\n')
        file.write('|   %*s  |' % (self.__mesh__.freedom*(len1 + 1) + self.__mesh__.freedom + len2, 'max:'))
        for i in range(0, len(self.__results__)):
            if self.__results__[i].t == t:
                file.write(' %+*.*E ' % (self.__params__.width, self.__params__.precision, self.__results__[i].max()))
                file.write('|')
        file.write('\n\n\n')

    # Визуализация заданной функции
    def plot(self, fun_name, t=0):
        # Проверка корректности задания времени
        if self.__params__.problem_type == 'dynamic' and \
                ((t < self.__params__.t0 or t > self.__params__.t1) or t % self.__params__.th > eps):
            error('Error: incorrectly specified the time: %5.2f' % t)
            return
        # Визуализация результата
        if self.__mesh__.fe_type == 'fe_2d_3':
            self.__plot_2d_tri__(fun_name, t)
        elif self.__mesh__.fe_type == 'fe_3d_4':
            self.__plot_3d_tet__(fun_name, t)

    # Визуализация заданной функции в случае плоской треугольной сетки
    def __plot_2d_tri__(self, fun_name, t=0):
        # Поиск индекса функции в списке результатов
        index = -1
        for i in range(0, len(self.__results__)):
            if self.__results__[i].name == fun_name and self.__results__[i].t == t:
                index = i
                break
        if index == -1:
            error('Error: \'%s\' is not a recognized function name' % fun_name)
            return

        plt.figure()
        plt.gca().set_aspect('equal')

        min_u = self.__results__[index].min()
        max_u = self.__results__[index].max()
        triang = tri.Triangulation(self.__mesh__.x, self.__mesh__.y)
        refiner = tri.UniformTriRefiner(triang)
        tri_refi, z_test_refi = refiner.refine_field(self.__results__[index].results, subdiv=3)

        plt.triplot(triang, lw=0.5, color='white')
        levels = np.arange(min_u, max_u, (max_u - min_u)/16.0)
#        cmap = cm.get_cmap(name='terrain', lut=None)
        cmap = cm.get_cmap(name='spectral', lut=None)
        plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap)
#        plt.tricontour(tri_refi, z_test_refi, levels=levels,
#                       colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
#                       linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])

        plt.colorbar()
        if self.__params__.problem_type == 'dynamic':
            fun_name += ' (t = %5.2f)' % t
        plt.title(fun_name)
        plt.show()

    # Визуализация заданной функции в случае КЭ в форме тетраэдра
    def __plot_3d_tet__(self, fun_name, t=0):
        # Поиск индекса функции в списке результатов
        index = -1
        for i in range(0, len(self.__results__)):
            if self.__results__[i].name == fun_name and self.__results__[i].t == t:
                index = i
                break
        if index == -1:
            error('Error: \'%s\' is not a recognized function name' % fun_name)
            return

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        surf = ax.plot_trisurf(self.__mesh__.x, self.__mesh__.y, self.__mesh__.z, triangles=self.__mesh__.surface)
        ax.set_zlim(min(self.__mesh__.z), max(self.__mesh__.z))
        plt.colorbar(surf)

        if self.__params__.problem_type == 'dynamic':
            fun_name += ' (t = %5.2f)' % t
        plt.title(fun_name)
        plt.show()
