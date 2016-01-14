#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
# Реализация вычислений задач заданного типа методо конечных элементов
#######################################################################

import math
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve, bicgstab, ArpackError
from abc import abstractmethod
from fem_mesh import TMesh
from fem_params import TFEMParams
from fem_progress import TProgress
from fem_fe import TFE, TFE1D2, TFE2D3, TFE2D4, TFE3D4, TFE3D8
from fem_parser import TParser
from fem_defs import DIR_X, DIR_Y, DIR_Z
from fem_error import TFEMException


# Абстрактный базовый класс, реализующий МКЭ
class TFEM:
    def __init__(self):
        self.__mesh__ = TMesh()                         # Дискретная модель объекта
        self.__params__ = TFEMParams()                  # Параметры расчета
        self.__global_matrix__ = lil_matrix((0, 0))     # Глобальная матрица жесткости (ГМЖ)
        self.__global_vector__ = []                     # Глобальный вектор нагрузок (правая часть)
        self.__progress__ = TProgress()                 # Индикатор прогресса расчета
        self.__result__ = []                            # Список результатов расчета для перемещений, деформаций, ...

    # Запуск процедуры расчета
    def calc(self):
        try:
            ret = self.__calc_problem__()
        except TFEMException as err:
            ret = False
            err.print_error()
        return ret

    # Добавление локальной матрицы жесткости (масс, демпфирования) к глобальной
    @abstractmethod
    def __assembly__(self, fe, index):
        raise NotImplementedError('Method TFEM.__assembly__ is pure virtual')

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...) 
    @abstractmethod
    def calc_results(self):
        raise NotImplementedError('Method TFEM.calc_results is pure virtual')

    # Решение СЛАУ
    def solve(self):
        ret = False
        if self.__params__.solve_method == 'direct':
            ret = self.solve_direct()
        elif self.__params__.solve_method == 'iterative':
            ret = self.solve_iterative()
        return ret

    # Прямое решение СЛАУ
    def solve_direct(self):
        self.__progress__.set_process('Solving of equation system...', 1, 1)
        self.__global_matrix__ = self.__global_matrix__.tocsr()
        try:
            self.__global_vector__ = spsolve(self.__global_matrix__, self.__global_vector__)
        except ArpackError:
            return False
        self.__progress__.set_progress(1)
        return True

    # Приближенное решение СЛАУ
    def solve_iterative(self):
        self.__progress__.set_process('Solving of equation system...', 1, 1)
        self.__global_matrix__ = self.__global_matrix__.tocsr()
        self.__global_vector__, info = bicgstab(self.__global_matrix__, self.__global_vector__, self.__global_vector__,
                                                self.__params__.eps)
        self.__progress__.set_progress(1)
        return True if not info else False

    # Создание нужного типа КЭ
    def create_fe(self):
        fe = TFE()
        if self.__mesh__.fe_type == 'fe_1d_2':
            fe = TFE1D2()
        elif self.__mesh__.fe_type == 'fe_2d_3':
            fe = TFE2D3()
        elif self.__mesh__.fe_type == 'fe_2d_4':
            fe = TFE2D4()
        elif self.__mesh__.fe_type == 'fe_3d_4':
            fe = TFE3D4()
        elif self.__mesh__.fe_type == 'fe_3d_8':
            fe = TFE3D8()
        return fe

    # Вычисление длины (площади) заданного граничного элемента
    def square(self, index):
        x = [0]*len(self.__mesh__.surface[index])
        y = [0]*len(self.__mesh__.surface[index])
        z = [0]*len(self.__mesh__.surface[index])
        for i in range(0, len(self.__mesh__.surface[index])):
            x[i], y[i], z[i] = self.__mesh__.get_coord(self.__mesh__.surface[index][i])
        s = 0
        if len(x) == 2:     # Граничный элемент - отрезок
            s = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2)
        elif len(x) == 3:   # Граничный элемент - треугольник
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif len(x) == 4:   # Граничный элемент - четырехугольник
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            s = math.sqrt(p*(p - a)*(p - b)*(p - c))

            a = math.sqrt((x[0] - x[3])**2 + (y[0] - y[3])**2 + (z[0] - z[3])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[3])**2 + (y[2] - y[3])**2 + (z[2] - z[3])**2)
            p = 0.5*(a + b + c)
            s += math.sqrt(p*(p - a)*(p - b)*(p - c))
        return s

    # Вычисление объема (длины, площади) заданного конечного элемента
    def volume(self, index):
        x = [0]*len(self.__mesh__.fe[index])
        y = [0]*len(self.__mesh__.fe[index])
        z = [0]*len(self.__mesh__.fe[index])
        for i in range(0, len(self.__mesh__.fe[index])):
            x[i], y[i], z[i] = self.__mesh__.get_coord(self.__mesh__.fe[index][i])
        v = 0
        if self.__mesh__.fe_type == 'fe_1d_2':
            v = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2)
        elif self.__mesh__.fe_type == 'fe_2d_3':
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            v = math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.__mesh__.fe_type == 'fe_2d_4':
            a = math.sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2 + (z[2] - z[1])**2)
            p = 0.5*(a + b + c)
            v = math.sqrt(p*(p - a)*(p - b)*(p - c))
            a = math.sqrt((x[0] - x[3])**2 + (y[0] - y[3])**2 + (z[0] - z[3])**2)
            b = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2 + (z[0] - z[2])**2)
            c = math.sqrt((x[2] - x[3])**2 + (y[2] - y[3])**2 + (z[2] - z[3])**2)
            p = 0.5*(a + b + c)
            v += math.sqrt(p*(p - a)*(p - b)*(p - c))
        elif self.__mesh__.fe_type == 'fe_3d_4':
            a = (x[1] - x[0])*(y[2] - y[0])*(z[3] - z[0]) + (x[3] - x[0])*(y[1] - y[0])*(z[2] - z[0]) + \
                (x[2] - x[0])*(y[3] - y[0])*(z[1] - z[0])
            b = (x[3] - x[0])*(y[2] - y[0])*(z[1] - z[0]) + (x[2] - x[0])*(y[1] - y[0])*(z[3] - z[0]) + \
                (x[1] - x[0])*(y[3] - y[0])*(z[2] - z[0])
            v = math.fabs(a - b)/6.0
        elif self.__mesh__.fe_type == 'fe_3d_8':
            ref = [[0, 1, 4, 7], [4, 1, 5, 7], [1, 2, 6, 7], [1, 5, 6, 7], [1, 2, 3, 7], [0, 3, 1, 7]]
            for i in range(0, 6):
                a = (x[ref[i][1]] - x[ref[i][0]])*(y[ref[i][2]] - y[ref[i][0]])*(z[ref[i][3]] - z[ref[i][0]]) + \
                    (x[ref[i][3]] - x[ref[i][0]])*(y[ref[i][1]] - y[ref[i][0]])*(z[ref[i][2]] - z[ref[i][0]]) + \
                    (x[ref[i][2]] - x[ref[i][0]])*(y[ref[i][3]] - y[ref[i][0]])*(z[ref[i][1]] - z[ref[i][0]])
                b = (x[ref[i][3]] - x[ref[i][0]])*(y[ref[i][2]] - y[ref[i][0]])*(z[ref[i][1]] - z[ref[i][0]]) + \
                    (x[ref[i][2]] - x[ref[i][0]])*(y[ref[i][1]] - y[ref[i][0]])*(z[ref[i][3]] - z[ref[i][0]]) + \
                    (x[ref[i][1]] - x[ref[i][0]])*(y[ref[i][3]] - y[ref[i][0]])*(z[ref[i][2]] - z[ref[i][0]])
                v += math.fabs(a - b)/6.0
        return v

    # Учет граничных условий
    def use_boundary_condition(self):
        parser = self.create_parser()
        counter = 0
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'boundary':
                counter += 1
        self.__progress__.set_process('Use of boundary conditions...', 1, counter*len(self.__mesh__.x))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'boundary':
                for j in range(0, len(self.__mesh__.x)):
                    self.__progress__.set_progress(counter)
                    counter += 1
                    x, y, z = self.__mesh__.get_coord(j)
                    parser.set_variable(self.__params__.names[0], x)
                    parser.set_variable(self.__params__.names[1], y)
                    parser.set_variable(self.__params__.names[2], z)
                    if len(self.__params__.bc_list[i].predicate):
                        parser.set_code(self.__params__.bc_list[i].predicate)
                        if parser.run() == 0:
                            continue
                    parser.set_code(self.__params__.bc_list[i].expression)
                    val = parser.run()
                    direct = self.__params__.bc_list[i].direct
                    if direct & DIR_X:
                        self.set_boundary_condition(j, 0, val)
                    if direct & DIR_Y:
                        self.set_boundary_condition(j, 1, val)
                    if direct & DIR_Z:
                        self.set_boundary_condition(j, 2, val)

    # Проверка соответствия граничного элемента предикату отбора (всех его вершин)
    def check_boundary_elements(self, i, predicate):
        if not len(predicate):
            return True
        parser = self.create_parser()
        for k in range(0, len(self.__mesh__.surface[0])):
            x, y, z = self.__mesh__.get_coord(self.__mesh__.surface[i][k])
            parser.set_variable(self.__params__.names[0], x)
            parser.set_variable(self.__params__.names[1], y)
            parser.set_variable(self.__params__.names[2], z)
            parser.set_code(predicate)
            if parser.run() == 0:
                return False
        return True

    # Задание граничных условий
    def set_boundary_condition(self, i, j, val):
        l = i*self.__mesh__.freedom + j
        for k in self.__global_matrix__[l].nonzero()[1]:
            if l != k:
                self.__global_matrix__[l, k] = self.__global_matrix__[k, l] = 0
        self.__global_vector__[l] = val*self.__global_matrix__[l, l]

    # Определение кол-ва результатов в зависимости от типа и размерности задачи
    def num_result(self):
        res = 0
        if self.__params__.problem_type == 'static':
            if self.__mesh__.freedom == 1:
                # u, Exx, Sxx
                res = 3
            elif self.__mesh__.freedom == 2:
                # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy
                res = 8
            elif self.__mesh__.freedom == 3:
                # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz
                res = 15
        elif self.__params__.problem_type == 'dynamic':
            if self.__mesh__.freedom == 1:
                # u, Exx, Sxx, ut, utt
                res = 5
            elif self.__mesh__.freedom == 2:
                # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
                res = 12
            elif self.__mesh__.freedom == 3:
                # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
                res = 21
        return res

    # Индекс функции в зависимости от типа и размерности задачи
    def index_result(self, i):
        ret = 0
        # u, Exx, Sxx
        index1 = [4, 7, 13]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy
        index2 = [4, 5, 7, 8, 10, 13, 14, 16]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz
        index3 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        # u, Exx, Sxx, ut, utt
        index4 = [4, 7, 13, 19, 22]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
        index5 = [4, 5, 7, 8, 10, 13, 14, 16, 19, 20, 22, 23]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
        index6 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        if self.__params__.problem_type == 'static':
            if self.__mesh__.freedom == 1:
                ret = index1[i]
            elif self.__mesh__.freedom == 2:
                ret = index2[i]
            elif self.__mesh__.freedom == 3:
                ret = index3[i]
        elif self.__params__.problem_type == 'dynamic':
            if self.__mesh__.freedom == 1:
                ret = index4[i]
            elif self.__mesh__.freedom == 2:
                ret = index5[i]
            elif self.__mesh__.freedom == 3:
                ret = index6[i]
        return ret

    # Задание сетки
    def set_mesh(self, mesh):
        self.__mesh__ = mesh

    # Задание параметров расчета
    def set_params(self, params):
        self.__params__ = params

    # Возврат результатов расчета
    def get_result(self):
        return self.__result__

    # Настройка парсера
    def create_parser(self):
        parser = TParser()
        for i in range(0, len(self.__params__.names)):
            parser.add_variable(self.__params__.names[i])
        for key, value in self.__params__.var_list.items():
            parser.add_variable(key, value)
        return parser
