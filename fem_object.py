#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#                       Описание объекта расчета
###################################################################

import os
import sys
from fem_mesh import TMesh
from fem_params import TFEMParams
from fem_error import TFEMException
from fem_static import TFEMStatic


class TObject:
    def __init__(self):
        self.params = TFEMParams()                      # Параметры расчета
        self.mesh = TMesh()                             # КЭ-модель
        self.result = []                                # Список результатов расчета для перемещений, деформаций, ...

    def set_mesh(self, name):
        self.mesh.load(name)
        print('Object: %s' % self.object_name())
        print('Points: %d' % len(self.mesh.x))
        print('FE: %d - %s' % (len(self.mesh.fe), self.mesh.fe_name(self.mesh.fe_type)))
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

    def add_initial_condition(self, e, p, d):
        self.params.add_initial_condition(e, p, d)

    def add_volume_load(self, e, p, d):
        self.params.add_volume_load(e, p, d)

    def add_surface_load(self, e, p, d):
        self.params.add_surface_load(e, p, d)

    def add_concentrated_load(self, e, p, d):
        self.params.add_concentrated_load(e, p, d)

    def add_variable(self, var, val):
        self.params.add_variable(var, val)

    def calc(self):
        if self.params.solve_method == '':
            raise TFEMException('solve_method_err')
        if self.params.problem_type == '':
            raise TFEMException('problem_type_err')
        if not len(self.params.e) or self.params.e[0] == 0 or self.params.m[0] == 0:
            raise TFEMException('elasticity_err')
        if self.params.problem_type == 'dynamic':
            if self.params.t0 == self.params.t1 or self.params.th <= 0:
                raise TFEMException('time_err')
        if self.params.problem_type == 'static':
            fem = TFEMStatic()
#        elif self.params.problem_type == 'dynamic':
#            pass
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
            raise TFEMException('read_file_err')
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
            file.write(' %*s |' % (len1, self.result[i].name))
        file.write('\n')
        file.write('|   %*s  |' % (self.mesh.freedom*(len1 + 1) + self.mesh.freedom + len2, 'min:'))
        for i in range(0, len(self.result)):
            file.write(' %+*.*E ' % (self.params.width, self.params.precision, self.result[i].min()))
            file.write('|')
        file.write('\n')
        file.write('|   %*s  |' % (self.mesh.freedom*(len1 + 1) + self.mesh.freedom + len2, 'max:'))
        for i in range(0, len(self.result)):
            file.write(' %+*.*E ' % (self.params.width, self.params.precision, self.result[i].max()))
            file.write('|')
        file.write('\n')
        file.close()
