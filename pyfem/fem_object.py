#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#            Реализация расчета задачи с помощью МКЭ
###################################################################

import os
import math
from fem_defs import eps, DIR_X, DIR_Y, DIR_Z
from fem_mesh import TMesh
from fem_params import TFEMParams
from fem_error import TFEMException
from fem_parser import TParser
from fem_fe import TFE, TFE1D2


class TObject:
    def __init__(self):
        self.file_name = ''                 # Имя файла с данными о геометрической модели
        self.object_name = ''               # Название объекта
        self.params = TFEMParams()          # Параметры расчета
        self.mesh = TMesh()                 # КЭ-модель
        self.__volume_force__ = []          # Узловые объемная, поверхностная и сосредоточенная нагрузки
        self.__surface_force__ = []
        self.__concentrated_force__ = []
        self.__fe__ = TFE()                 # Конечный элемент
        self.__global_matrix__ = []         # Глобальная матрица жесткости (ГМЖ)
        self.__global_vector__ = []         # Глобальная правая часть

    def set_mesh(self, name):
        self.file_name = name
        self.object_name = os.path.basename(name)
        self.mesh.load(name)

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

    def add_volume_condition(self, e, p, d):
        self.params.add_volume_condition(e, p, d)

    def add_surface_condition(self, e, p, d):
        self.params.add_surface_condition(e, p, d)

    def add_concentrated_condition(self, e, p, d):
        self.params.add_concentrated_condition(e, p, d)

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
            self.calc_static_problem()
        elif self.params.problem_type == 'dynamic':
            self.calc_dynamic_problem()

    # Добавление локальной матрицы жесткости (ЛМЖ) к ГМЖ
    def __ansamble__(self, index):
        # Добавление матрицы
        for i in range(0, self.__fe__.size):
            k = self.mesh.fe[index][i//self.mesh.freedom]*self.mesh.freedom + i % self.mesh.freedom
            for j in range(i, self.__fe__.size):
                l = self.mesh.fe[index][j//self.mesh.freedom]*self.mesh.freedom + j % self.mesh.freedom
                self.__global_matrix__[k][l] += self.__fe__.K[i][j]
                if k != l:
                    self.__global_matrix__[l][k] += self.__fe__.K[i][j]
            self.__global_vector__[k] += self.__fe__.K[i][self.__fe__.size]

    # Создание нужного типа КЭ
    def create_fe(self):
        if self.mesh.fe_type == 'fe_1d_2':
            self.__fe__ = TFE1D2()

    # Расчет статической задачи методом конечных элементов
    def calc_static_problem(self):
        # Создание ГМЖ
        for i in range(0, len(self.mesh.x)*self.mesh.freedom):
            row = [0]*len(self.mesh.x)*self.mesh.freedom
            self.__global_matrix__.append(row)
        self.__global_vector__ = [0]*len(self.mesh.x)*self.mesh.freedom

        self.create_fe()
        self.__fe__.set_elasticity(self.params.e, self.params.m)
        # Вычисление компонент нагрузки
        self.prepare_force()
        # Формирование глобальной матрицы жесткости
        for i in range(0, len(self.mesh.fe)):
            x = [0]*len(self.mesh.fe[i])
            y = [0]*len(self.mesh.fe[i])
            z = [0]*len(self.mesh.fe[i])
            vx = [0]*len(self.mesh.fe[i])*self.mesh.freedom
            vy = [0]*len(self.mesh.fe[i])*self.mesh.freedom
            vz = [0]*len(self.mesh.fe[i])*self.mesh.freedom
            # Настройка КЭ
            for j in range(len(self.mesh.fe[i])):
                x[j] = self.mesh.x[self.mesh.fe[i][j]]
                y[j] = self.mesh.y[self.mesh.fe[i][j]] if (len(self.mesh.y)) else 0
                z[j] = self.mesh.z[self.mesh.fe[i][j]] if (len(self.mesh.z)) else 0
                vx[j] = self.__volume_force__[self.mesh.fe[i][j]*self.mesh.freedom + 0]
                vy[j] = self.__volume_force__[self.mesh.fe[i][j]*self.mesh.freedom + 1] if (len(self.mesh.y)) else 0
                vz[j] = self.__volume_force__[self.mesh.fe[i][j]*self.mesh.freedom + 2] if (len(self.mesh.z)) else 0
            self.__fe__.set_coord(x, y, z)
            self.__fe__.set_volume_force(vx, vy, vz)
            self.__fe__.generate()
            # Ансамблирование ЛМЖ к ГМЖ
            self.__ansamble__(i)
        # Учет сосредоточенной и поверхностной нагрузок
        self.use_force_condition()
        # Учет краевых условий
        self.use_boundary_condition()
        # Решение СЛАУ
        if not self.solve():
            print('The system of equations is not solved!')
        else:
            print(self.__global_vector__)

    # Предварительное вычисление нагрузок
    def prepare_force(self):
        parser = TParser()
        self.__volume_force__ = [0]*len(self.mesh.x)
        self.__surface_force__ = [0]*len(self.mesh.x)
        self.__concentrated_force__ = [0]*len(self.mesh.x)
        for i in range(0, len(self.params.names)):
            parser.add_variable(self.params.names[i])
        for i in range(0, len(self.params.bc_list)):
            if not (self.params.bc_list[i].type == 'volume' or self.params.bc_list[i].type == 'surface' or
               self.params.bc_list[i].type == 'concentrated'):
                continue
            for j in range(0, len(self.mesh.x)):
                x = self.mesh.x[j]
                y = self.mesh.y[j] if (len(self.mesh.y)) else 0
                z = self.mesh.z[j] if (len(self.mesh.z)) else 0
                parser.set_variable(self.params.names[0], x)
                parser.set_variable(self.params.names[1], y)
                parser.set_variable(self.params.names[2], z)
                if len(self.params.bc_list[i].predicate):
                    parser.set_code(self.params.bc_list[i].predicate)
                    if parser.error != '':
                        return
                    if parser.run() == 0:
                        continue
                parser.set_code(self.params.bc_list[i].expression)
                val = parser.run()
                index = -1
                if self.params.bc_list[i].direct & DIR_X:
                    index = j*self.mesh.freedom + 0
                    self.add_force(i, index, val)
                if self.params.bc_list[i].direct & DIR_Y:
                    index = j*self.mesh.freedom + 1
                    self.add_force(i, index, val)
                if self.params.bc_list[i].direct & DIR_Z:
                    index = j*self.mesh.freedom + 2
                    self.add_force(i, index, val)

    # Учет соответствующей нагрузки
    def add_force(self, i, j, value):
        if self.params.bc_list[i].type == 'volume':
            self.__volume_force__[j] += value
        elif self.params.bc_list[i].type == 'surface':
            self.__surface_force__[j] += value
        elif self.params.bc_list[i].type == 'concentrated':
            self.__concentrated_force__[j] += value

    # Учет сосредоточенной и поверхностной нагрузок
    def use_force_condition(self):
        # Учет сосредоточенной нагрузки
        for i in range(0, len(self.__concentrated_force__)):
            if self.__concentrated_force__[i]:
                self.__global_vector__[i] += self.__concentrated_force__[i]
        # Учет поверхностной нагрузки
        if self.mesh.freedom == 1:
            return
        x = [0]*len(self.mesh.surface[0])
        y = [0]*len(self.mesh.surface[0])
        z = [0]*len(self.mesh.surface[0])
        for i in range(0, len(self.mesh.surface)):
            for j in range(0, len(self.mesh.surface[0])):
                x[j] = self.mesh.x[self.mesh.surface[i][j]]
                y[j] = self.mesh.y[self.mesh.surface[i][j]]
                z[j] = self.mesh.z[self.mesh.surface[i][j]] if (len(self.mesh.z)) else 0
            se = self.square(x, y, z)
            for j in range(0, len(self.mesh.surface[0])):
                for k in range(0, self.mesh.freedom):
                    l = self.mesh.surface[i][j]*self.mesh.freedom + k
                    if self.__surface_force__[l]:
                        if self.mesh.fe_type == 'fe_3d_10':
                            self.__global_vector__[l] += self.__surface_force__[k]*se/float(len(self.mesh.surface[0]))
                        elif j > 2:
                            self.__global_vector__[l] += self.__surface_force__[k]*se/3.0

    # Вычисление длины (площади) граничного элемента
    @staticmethod
    def square(x, y, z):
        s = 0
        if len(x) == 2:     # Граничный элемент - отрезок
            s = math.sqrt((x[0] - x[1])*(x[0] - x[1]) + (y[0] - y[1])*(y[0] - y[1]))
        elif len(x) == 3:   # Граничный элемент - треугольник
            a = math.sqrt((x[0] - x[1])*(x[0] - x[1]) + (y[0] - y[1])*(y[0] - y[1]) + (z[0] - z[1])*(z[0] - z[1]))
            b = math.sqrt((x[0] - x[2])*(x[0] - x[2]) + (y[0] - y[2])*(y[0] - y[2]) + (z[0] - z[2])*(z[0] - z[2]))
            c = math.sqrt((x[2] - x[1])*(x[2] - x[1]) + (y[2] - y[1])*(y[2] - y[1]) + (z[2] - z[1])*(z[2] - z[1]))
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif len(x) == 4:   # Граничный элемент - четырехугольник
            a = math.sqrt((x[0] - x[1])*(x[0] - x[1]) + (y[0] - y[1])*(y[0] - y[1]) + (z[0] - z[1])*(z[0] - z[1]))
            b = math.sqrt((x[0] - x[2])*(x[0] - x[2]) + (y[0] - y[2])*(y[0] - y[2]) + (z[0] - z[2])*(z[0] - z[2]))
            c = math.sqrt((x[2] - x[1])*(x[2] - x[1]) + (y[2] - y[1])*(y[2] - y[1]) + (z[2] - z[1])*(z[2] - z[1]))
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))

            a = math.sqrt((x[0] - x[3])*(x[0] - x[3]) + (y[0] - y[3])*(y[0] - y[3]) + (z[0] - z[3])*(z[0] - z[3]))
            b = math.sqrt((x[0] - x[2])*(x[0] - x[2]) + (y[0] - y[2])*(y[0] - y[2]) + (z[0] - z[2])*(z[0] - z[2]))
            c = math.sqrt((x[2] - x[3])*(x[2] - x[3]) + (y[2] - y[3])*(y[2] - y[3]) + (z[2] - z[3])*(z[2] - z[3]))
            p = 0.5*(a + b + c)
            s += math.sqrt(p*(p - a)*(p - b)*(p - c))
        return s

    # Учет граничных условий
    def use_boundary_condition(self):
        parser = TParser()
        for i in range(0, len(self.params.names)):
            parser.add_variable(self.params.names[i])
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'boundary':
                for j in range(0, len(self.mesh.x)):
                    x = self.mesh.x[j]
                    y = self.mesh.y[j] if (len(self.mesh.y)) else 0
                    z = self.mesh.z[j] if (len(self.mesh.z)) else 0
                    parser.set_variable(self.params.names[0], x)
                    parser.set_variable(self.params.names[1], y)
                    parser.set_variable(self.params.names[2], z)
                    if len(self.params.bc_list[i].predicate):
                        parser.set_code(self.params.bc_list[i].predicate)
                        if parser.error != '' or parser.run() == 0:
                            continue
                    parser.set_code(self.params.bc_list[i].expression)
                    val = parser.run()
                    direct = self.params.bc_list[i].direct
                    if direct & DIR_X:
                        self.set_boundary_condition(j, 0, val)
                    if direct & DIR_Y:
                        self.set_boundary_condition(j, 1, val)
                    if direct & DIR_Z:
                        self.set_boundary_condition(j, 2, val)

    # Задание граничных условий
    def set_boundary_condition(self, i, j, val):
        for k in range(0, len(self.mesh.x)*self.mesh.freedom):
            l = i*self.mesh.freedom + j
            if l != k:
                self.__global_matrix__[l][k] = self.__global_matrix__[k][l] = 0
            self.__global_vector__[l] = val*self.__global_matrix__[l][l]

    # Расчет динамической задачи методом конечных элементов
    def calc_dynamic_problem(self):
        pass

    # Решение СЛАУ
    def solve(self):
        ret = False
        if self.params.solve_method == 'direct':
            ret = self.solve_direct()
        elif self.params.solve_method == 'iterative':
            ret = self.solve_iterative()
        return ret

    # Решение СЛАУ методом Гаусса
    def solve_direct(self):
        size = len(self.__global_matrix__)
        result = [0]*size
        # Прямой ход метода Гаусса
        for i in range(0, size - 1):
            if math.fabs(self.__global_matrix__[i][i]) < eps:
                for l in range(i + 1, size):
                    if math.fabs(self.__global_matrix__[l][i]) < eps:
                        continue
                    for j in range(0, size):
                        self.__global_matrix__[l][j], self.__global_matrix__[i][j] = \
                            self.__global_matrix__[i][j], self.__global_matrix__[l][j]
                    self.__global_vector__[l], self.__global_vector__[i] = \
                        self.__global_vector__[i], self.__global_vector__[l]
            val1 = self.__global_matrix__[i][i]
            for j in range(i + 1, size):
                val2 = self.__global_matrix__[j][i]
                if math.fabs(val2) < eps:
                    continue
                for k in range(i, size):
                    self.__global_matrix__[j][k] -= val2*self.__global_matrix__[i][k]/val1
                self.__global_vector__[j] -= val2*self.__global_vector__[i]/val1
            if math.fabs(self.__global_matrix__[size - 1][size - 1]) < eps:
                return False
        # Обратный ход метода Гаусса
        result[size - 1] = self.__global_vector__[size - 1]/self.__global_matrix__[size - 1][size - 1]
        for k in range(0, size - 1):
            i = size - k - 2
            s = self.__global_vector__[i]
            for j in range(i + 1, size):
                s -= result[j]*self.__global_matrix__[i][j]
                if math.fabs(self.__global_matrix__[i][i]) < eps:
                    return False
                result[i] = s/self.__global_matrix__[i][i]
        self.__global_vector__ = result
        return True

    # Решение СЛАУ методом простой итерации
    def solve_iterative(self):
        return True

