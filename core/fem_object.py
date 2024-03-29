#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#                       Описание объекта расчета
###################################################################

import os
import sys
import json
from timeit import default_timer as timer
from datetime import *
from math import fabs
from core.fem_mesh import TMesh, TMeshTRPA
from core.fem_fem import TFEM
from core.fem_params import TFEMParams
from core.fem_static import TFEMStatic
from core.fem_dynamic import TFEMDynamic
from core.fem_error import TException, print_error


class TObject:
    def __init__(self):
        self.__fem = TFEM()             # Метод расчета
        self.__params = TFEMParams()    # Параметры расчета
        self.__mesh = TMesh()           # КЭ-модель

    def set_mesh(self, name):
        try:
            if os.path.splitext(name)[1].upper() == '.TRPA':
                self.__mesh = TMeshTRPA()
                self.__mesh.load(name)
                print('Object: %s' % self.object_name())
                print('Points: %d' % len(self.__mesh.x))
                print('FE: %d - %s' % (len(self.__mesh.fe), self.__mesh.fe_name()))
            else:
                print_error('Unknown file extension!')
                return False
        except TException as err:
            err.print_error()
            return False
        return True

    # Название объекта
    def object_name(self):
        return os.path.splitext(os.path.basename(self.__mesh.mesh_file))[0]

    def set_problem_type(self, problem_type):
        self.__params.problem_type = problem_type

    def set_solve_method(self, solve_method):
        self.__params.solve_method = solve_method

    def set_eps(self, e):
        self.__params.eps = e

    def set_output(self, width, precision):
        self.__params.width = width
        self.__params.precision = precision

    def set_names(self, names):
        self.__params.names = names

    def set_time(self, t0, t1, th):
        self.__params.t0 = t0
        self.__params.t1 = t1
        self.__params.th = th

    def add_density(self, e, p=None):
        self.__params.add_density(e, p)

    def add_damping(self, e, p=None):
        self.__params.add_damping(e, p)

    def add_temperature(self, e, p=None):
        self.__params.add_temperature(e, p)

    def add_alpha(self, e, p=None):
        self.__params.add_alpha(e, p)

    def add_young_modulus(self, e, p=None):
        self.__params.add_young_modulus(e, p)

    def add_shear_modulus(self, e, p=None):
        self.__params.add_shear_modulus(e, p)

    def add_poisson_ratio(self, e, p=None):
        self.__params.add_poisson_ratio(e, p)

    def add_thickness(self, e, p=None):
        self.__params.add_thickness(e, p)

    def add_boundary_condition(self, d, e, p=None):
        self.__params.add_boundary_condition(e, p, d)

    def add_initial_condition(self, d, e):
        self.__params.add_initial_condition(e, None, d)

    def add_volume_load(self, d, e, p=None):
        self.__params.add_volume_load(e, p, d)

    def add_concentrated_load(self, d, e, p=None):
        self.__params.add_concentrated_load(e, p, d)

    def add_surface_load(self, d, e, p=None):
        self.__params.add_surface_load(e, p, d)

    def add_pressure_load(self, e, p=None):
        self.__params.add_pressure_load(e, p)

    def add_variable(self, var, val):
        self.__params.add_variable(var, val)

    def calc(self):
        ret = False
        start = timer()
        if self.__params.problem_type == 'static':
            self.__fem = TFEMStatic()
        elif self.__params.problem_type == 'dynamic':
            self.__fem = TFEMDynamic()
        self.__fem.set_mesh(self.__mesh)
        self.__fem.set_params(self.__params)
        try:
            ret = self.__fem.calc()
            if ret:
                print('Lead time %f sec' % (timer() - start))
        except TException as err:
            err.print_error()
        return ret

    # Вывод результатов расчета
    def print_result(self, *argv):
        file = sys.stdout
        try:
            if len(argv) == 1:
                file = open(argv[0], 'w')
        except IOError:
            print_error('Error: unable to open file %s' % argv[0])
            return
        if self.__params.problem_type == 'static':
            self.__print(file)
        else:
            t = self.__params.t0
            while t <= self.__params.t1:
                file.write('t = %5.2f\n' % t)
                self.__print(file, t)
                t += self.__params.th
                if fabs(t - self.__params.t1) < self.__params.eps:
                    t = self.__params.t1
        if len(argv) == 1:
            file.close()

    # Запись результатов расчета в файл (json)
    def save_result(self, file_name):
        if len(file_name) < 5 or file_name[len(file_name) - 5:] != '.json':
            file_name += '.json'
        # Определяем текущее время и дату
        now = datetime.now()
        dt = str(now.day) + '.' + str(now.month) + '.' + str(now.year)
        tm = str(now.hour) + ':' + str(now.minute) + ':' + str(now.second)
        # Формирование структуры файла
        header = dict()
        header['name'] = self.__mesh.mesh_file[self.__mesh.mesh_file.rfind('/') + 1:
                                                self.__mesh.mesh_file.rfind('.')]
        header['type'] = self.__params.problem_type
        header['date_time'] = dt + " " + tm

        mesh = dict()
        mesh['fe_type'] = self.__mesh.fe_type
        mesh['vertex'] = self.__mesh.x
        mesh['fe'] = self.__mesh.fe
        mesh['be'] = self.__mesh.be

        results = list()
        for i in range(0, len(self.__fem.get_results())):
            res = dict()
            res['function'] = self.__fem.get_results(i).name
            res['t'] = self.__fem.get_results(i).t
            res['results'] = self.__fem.get_results(i).results
            results.append(res)

        data = dict()
        data['header'] = header
        data['mesh'] = mesh
        data['results'] = results

        json_data = json.dumps(data)
        with open(file_name, 'w') as file:
            json.dump(json_data, file)

    # Вывод результатов расчета для одного момента времени
    def __print(self, file, t=0):
        # Определение ширины позиции
        len1 = len('%+*.*E' % (self.__params.width, self.__params.precision, 3.14159))
        len2 = len('%d' % len(self.__mesh.x))
        # Вывод заголовка
        file.write('| %*s  |' % (len2, 'N'))
        for i in range(0, self.__mesh.dimension):
            file.write(' %*s' % (len1, self.__params.names[i]))
            if i < self.__mesh.dimension - 1:
                file.write('|')
        file.write(' |')
        for i in range(0, len(self.__fem.get_results())):
            if self.__fem.get_results(i).t == t:
                file.write(' %*s |' % (len1, self.__fem.get_results(i).name))
        file.write('\n')
        for i in range(0, len(self.__mesh.x)):
            file.write('| %*d  |' % (len2, i + 1))
            file.write(' %+*.*E' % (self.__params.width, self.__params.precision, self.__mesh.x[i][0]))
            if self.__mesh.dimension > 1:
                file.write('| %+*.*E' % (self.__params.width, self.__params.precision, self.__mesh.x[i][1]))
            if self.__mesh.dimension > 2:
                file.write('| %+*.*E' % (self.__params.width, self.__params.precision, self.__mesh.x[i][2]))
            file.write(' | ')
            for k in range(0, len(self.__fem.get_results())):
                if self.__fem.get_results(k).t == t:
                    file.write('%+*.*E' %
                               (self.__params.width, self.__params.precision, self.__fem.get_results(k).results[i]))
                    file.write(' | ')
            file.write('\n')
        file.write('\n')
        # Печать итогов
        file.write('|  %*s  ' % (len2, ' '))
        for i in range(0, self.__mesh.dimension):
            file.write(' %*s' % (len1, ' '))
            if i < self.__mesh.dimension - 1:
                file.write(' ')
        file.write('  |')
        for i in range(0, len(self.__fem.get_results())):
            if self.__fem.get_results(i).t == t:
                file.write(' %*s |' % (len1, self.__fem.get_results(i).name))
        file.write('\n')
        file.write('|   %*s  |' % (self.__mesh.dimension * (len1 + 1) + self.__mesh.dimension + len2, 'min:'))
        for i in range(0, len(self.__fem.get_results())):
            if self.__fem.get_results(i).t == t:
                file.write(' %+*.*E ' % (self.__params.width, self.__params.precision, self.__fem.get_results(i).min()))
                file.write('|')
        file.write('\n')
        file.write('|   %*s  |' % (self.__mesh.dimension * (len1 + 1) + self.__mesh.dimension + len2, 'max:'))
        for i in range(0, len(self.__fem.get_results())):
            if self.__fem.get_results(i).t == t:
                file.write(' %+*.*E ' % (self.__params.width, self.__params.precision, self.__fem.get_results(i).max()))
                file.write('|')
        file.write('\n\n\n')
