#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#                       Описание объекта расчета
###################################################################

import os
import sys
import simplejson as json
from datetime import *
from math import fabs
from core.fem_mesh import TMesh
from core.fem_fem import TFEM
from core.fem_params import TFEMParams
from core.fem_static import TFEMStatic
from core.fem_dynamic import TFEMDynamic
from core.fem_error import TFEMException


# Вывод сообщения об ошибке
def print_error(err_msg):
    print('\033[1;31m%s\033[1;m' % err_msg)


class TObject:
    def __init__(self):
        self.__params__ = TFEMParams()          # Параметры расчета
        self.__mesh__ = TMesh()                 # КЭ-модель
        self.__results__ = []                   # Список результатов расчета для перемещений, деформаций, ...

    def set_mesh(self, name):
        try:
            self.__mesh__.load(name)
            print('Object: %s' % self.object_name())
            print('Points: %d' % len(self.__mesh__.x))
            print('FE: %d - %s' % (len(self.__mesh__.fe), self.__mesh__.fe_name()))
        except TFEMException as err:
            err.print_error()
            return False
        return True

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
            print_error('Error: unable to open file %s' % argv[0])
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
        header['name'] = self.__mesh__.mesh_file[self.__mesh__.mesh_file.rfind('/') + 1:
                                                 self.__mesh__.mesh_file.rfind('.')]
        header['type'] = self.__params__.problem_type
        header['date_time'] = dt + " " + tm

        mesh = dict()
        mesh['fe_type'] = self.__mesh__.fe_type
        mesh['vertex'] = self.__mesh__.x
        mesh['fe'] = self.__mesh__.fe
        mesh['be'] = self.__mesh__.be

        results = list()
        for i in range(0, len(self.__results__)):
            res = dict()
            res['function'] = self.__results__[i].name
            res['t'] = self.__results__[i].t
            res['results'] = self.__results__[i].results
            results.append(res)

        data = dict()
        data['header'] = header
        data['mesh'] = mesh
        data['results'] = results

        json_data = json.dumps(data)
        with open(file_name, 'w') as file:
            json.dump(json_data, file)

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
            file.write(' %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.x[i][0]))
            if self.__mesh__.freedom > 1:
                file.write(', %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.x[i][1]))
            if self.__mesh__.freedom > 2:
                file.write(', %+*.*E' % (self.__params__.width, self.__params__.precision, self.__mesh__.x[i][2]))
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
