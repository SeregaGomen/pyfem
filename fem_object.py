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
import numpy as np
from math import fabs
from fem_mesh import TMesh
from fem_fem import TFEM
from fem_params import TFEMParams
from fem_static import TFEMStatic
from fem_dynamic import TFEMDynamic


# Вывод сообщения об ошибке
def error(err_msg):
    print('\033[1;31m%s\033[1;m' % err_msg)


class TObject:
    def __init__(self):
        self.params = TFEMParams()                      # Параметры расчета
        self.mesh = TMesh()                             # КЭ-модель
        self.result = []                                # Список результатов расчета для перемещений, деформаций, ...

    def set_mesh(self, name):
        self.mesh.load(name)
        print('Object: %s' % self.object_name())
        print('Points: %d' % len(self.mesh.x))
        print('FE: %d - %s' % (len(self.mesh.fe), self.mesh.fe_name()))
#        print('FE type: %s' % self.mesh.fe_type)

    # Название объекта
    def object_name(self):
        return os.path.splitext(os.path.basename(self.mesh.mesh_file))[0]

    def set_problem_type(self, problem_type):
        self.params.problem_type = problem_type

    def set_solve_method(self, solve_method):
        self.params.solve_method = solve_method

    def set_eps(self, e):
        self.params.eps = e

    def set_width(self, width):
        self.params.width = width

    def set_precision(self, precision):
        self.params.precision = precision

    def set_elasticity(self, e, m):
        self.params.e = e
        self.params.m = m

    def set_density(self, density):
        self.params.density = density

    def set_time(self, t0, t1, th):
        self.params.t0 = t0
        self.params.t1 = t1
        self.params.th = th

    def set_damping(self, damping):
        self.params.damping = damping

    def set_names(self, names):
        self.params.names = names

    def add_boundary_condition(self, e, p, d):
        self.params.add_boundary_condition(e, p, d)

    def add_initial_condition(self, e, d):
        self.params.add_initial_condition(e, '', d)

    def add_volume_load(self, e, p, d):
        self.params.add_volume_load(e, p, d)

    def add_surface_load(self, e, p, d):
        self.params.add_surface_load(e, p, d)

    def add_concentrated_load(self, e, p, d):
        self.params.add_concentrated_load(e, p, d)

    def add_variable(self, var, val):
        self.params.add_variable(var, val)

    def calc(self):
        fem = TFEM()
        if self.params.problem_type == 'static':
            fem = TFEMStatic()
        elif self.params.problem_type == 'dynamic':
            fem = TFEMDynamic()
        fem.set_mesh(self.mesh)
        fem.set_params(self.params)
        ret = fem.calc()
        if ret:
            self.result = fem.get_result()
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
        if self.params.problem_type == 'static':
            self.__print__(file)
        else:
            t = self.params.t0
            while t <= self.params.t1:
                file.write('t = %5.2f\n' % t)
                self.__print__(file, t)
                t += self.params.th
                if fabs(t - self.params.t1) < self.params.eps:
                    t = self.params.t1

        file.close()

    # Вывод результатов расчета для одного момента времени
    def __print__(self, file, t=0):
        # Определение ширины позиции
        len1 = len('%+*.*E' % (self.params.width, self.params.precision, 3.14159))
        len2 = len('%d' % len(self.mesh.x))
        # Вывод заголовка
        file.write('| %*s  (' % (len2, 'N'))
        for i in range(0, self.mesh.freedom):
            file.write(' %*s' % (len1, self.params.names[i]))
            if i < self.mesh.freedom - 1:
                file.write(',')
        file.write(') |')
        for i in range(0, len(self.result)):
            if self.result[i].t == t:
                file.write(' %*s |' % (len1, self.result[i].name))
        file.write('\n')
        for i in range(0, len(self.mesh.x)):
            file.write('| %*d  (' % (len2, i + 1))
            file.write(' %+*.*E' % (self.params.width, self.params.precision, self.mesh.x[i]))
            if len(self.mesh.y):
                file.write(', %+*.*E' % (self.params.width, self.params.precision, self.mesh.y[i]))
            if len(self.mesh.z):
                file.write(', %+*.*E' % (self.params.width, self.params.precision, self.mesh.z[i]))
            file.write(') | ')
            for k in range(0, len(self.result)):
                if self.result[k].t == t:
                    file.write('%+*.*E' % (self.params.width, self.params.precision, self.result[k].results[i]))
                    file.write(' | ')
            file.write('\n')
        file.write('\n')
        # Печать итогов
        file.write('|  %*s  ' % (len2, ' '))
        for i in range(0, self.mesh.freedom):
            file.write(' %*s' % (len1, ' '))
            if i < self.mesh.freedom - 1:
                file.write(' ')
        file.write('  |')
        for i in range(0, len(self.result)):
            if self.result[i].t == t:
                file.write(' %*s |' % (len1, self.result[i].name))
        file.write('\n')
        file.write('|   %*s  |' % (self.mesh.freedom*(len1 + 1) + self.mesh.freedom + len2, 'min:'))
        for i in range(0, len(self.result)):
            if self.result[i].t == t:
                file.write(' %+*.*E ' % (self.params.width, self.params.precision, self.result[i].min()))
                file.write('|')
        file.write('\n')
        file.write('|   %*s  |' % (self.mesh.freedom*(len1 + 1) + self.mesh.freedom + len2, 'max:'))
        for i in range(0, len(self.result)):
            if self.result[i].t == t:
                file.write(' %+*.*E ' % (self.params.width, self.params.precision, self.result[i].max()))
                file.write('|')
        file.write('\n\n\n')

    # Визуализация заданной функции
    def plot(self, fun_name, t=0):
        # Поиск индекса функции в списке результатов
        index = -1
        for i in range(0, len(self.result)):
            if self.result[i].name == fun_name and self.result[i].t == t:
                index = i
                break
        if index == -1:
            error('Error: \'%s\' is not a recognized function name' % fun_name)
            return

        plt.figure()
        plt.gca().set_aspect('equal')

        min_u = self.result[index].min()
        max_u = self.result[index].max()
        triang = tri.Triangulation(self.mesh.x, self.mesh.y)
        refiner = tri.UniformTriRefiner(triang)
        tri_refi, z_test_refi = refiner.refine_field(self.result[index].results, subdiv=3)

        plt.triplot(triang, lw=0.5, color='white')
        levels = np.arange(min_u, max_u, (max_u - min_u)/16.0)
        cmap = cm.get_cmap(name='terrain', lut=None)
        plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap)
        plt.tricontour(tri_refi, z_test_refi, levels=levels,
                       colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
                       linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])

        if self.params.problem_type == 'dynamic':
            fun_name += ' (t = %5.2f)' % t
        plt.title(fun_name)
        plt.show()
