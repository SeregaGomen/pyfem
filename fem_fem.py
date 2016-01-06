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
from fem_result import TResult


# Абстрактный базовый класс, реализующий МКЭ
class TFEM:
    def __init__(self):
        self.__mesh__ = TMesh()                         # Дискретная модель объекта
        self.__params__ = TFEMParams()                  # Параметры расчета
        self.__volume_load__ = []                       # Узловые объемная, поверхностная и сосредоточенная нагрузки
        self.__surface_load__ = []
        self.__concentrated_load__ = []
        self.__fe__ = TFE()                             # Конечный элемент
        self.__global_matrix__ = lil_matrix((0, 0))     # Глобальная матрица жесткости (ГМЖ)
        self.__global_vector__ = []                     # Глобальная правая часть
        self.__progress__ = TProgress()                 # Индикатор прогресса расчета
        self.__result__ = []                            # Список результатов расчета для перемещений, деформаций, ...
        self.__surface_attribute__ = []                 # Признак того, что ГЭ нужно учитывать при вычислении нагрузки

    # Запуск процедуры расчета
    @abstractmethod
    def calc(self):
        pass

    # Добавление локальной матрицы жесткости (масс, демпфирования) к глобальной
    @abstractmethod
    def __assembly__(self, index):
        pass

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...) 
    def calc_results(self):
        # Выделяем память для хранения результатов
        res = []
        for i in range(0, self.num_result()):
            r = [0]*len(self.__mesh__.x)
            res.append(r)
        uvw = [0]*len(self.__mesh__.fe[0])*self.__mesh__.freedom
        counter = [0]*len(self.__mesh__.x)  # Счетчик кол-ва вхождения узлов для осреднения результатов
        # Копируем полученные перемещения
        for i in range(0, len(self.__mesh__.x)):
            for j in range(0, self.__mesh__.freedom):
                res[j][i] = self.__global_vector__[i*self.__mesh__.freedom + j]
        # Вычисляем стандартные результаты по всем КЭ
        self.__progress__.set_process('Calculation results...', 1, len(self.__mesh__.fe))
        for i in range(0, len(self.__mesh__.fe)):
            self.__progress__.set_progress(i + 1)
            x = [0]*len(self.__mesh__.fe[i])
            y = [0]*len(self.__mesh__.fe[i])
            z = [0]*len(self.__mesh__.fe[i])
            for j in range(len(self.__mesh__.fe[i])):
                x[j], y[j], z[j] = self.__mesh__.get_coord(self.__mesh__.fe[i][j])
            self.__fe__.set_coord(x, y, z)
            for j in range(0, len(self.__mesh__.fe[i])):
                for k in range(0, self.__mesh__.freedom):
                    uvw[j*self.__mesh__.freedom + k] = \
                        self.__global_vector__[self.__mesh__.freedom*self.__mesh__.fe[i][j] + k]
            r = self.__fe__.calc(uvw)
            for m in range(0, len(r)):
                for j in range(0, len(r[0])):
                    res[self.__mesh__.freedom + m][self.__mesh__.fe[i][j]] += r[m][j]
                    if not m:
                        counter[self.__mesh__.fe[i][j]] += 1
        # Осредняем результаты
        for i in range(self.__mesh__.freedom, self.num_result()):
            for j in range(0, len(self.__mesh__.x)):
                res[i][j] /= counter[j]
        # Сохраняем полученные результаты в списке
        for i in range(0, self.num_result()):
            r = TResult()
            r.name = self.__params__.names[self.index_result(i)]
            r.results = res[i]
            self.__result__.append(r)
        
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
        if self.__mesh__.fe_type == 'fe_1d_2':
            self.__fe__ = TFE1D2()
        elif self.__mesh__.fe_type == 'fe_2d_3':
            self.__fe__ = TFE2D3()
        elif self.__mesh__.fe_type == 'fe_2d_4':
            self.__fe__ = TFE2D4()
        elif self.__mesh__.fe_type == 'fe_3d_4':
            self.__fe__ = TFE3D4()
        elif self.__mesh__.fe_type == 'fe_3d_8':
            self.__fe__ = TFE3D8()

    # Учет нагрузки
    def add_load(self, i, j, value):
        if self.__params__.bc_list[i].type == 'volume':
            self.__volume_load__[j] += value
        elif self.__params__.bc_list[i].type == 'surface':
            self.__surface_load__[j] += value
        elif self.__params__.bc_list[i].type == 'concentrated':
            self.__concentrated_load__[j] += value

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
        for i in range(0, len(self.__params__.names)):
            parser.add_variable(self.__params__.names[i])
        for key, value in self.__params__.var_list.items():
            parser.add_variable(key, value)
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
                        if parser.error != '':
                            return
                        if parser.run() == 0:
                            continue
                    parser.set_code(self.__params__.bc_list[i].expression)
                    if parser.error != '':
                        return
                    val = parser.run()
                    direct = self.__params__.bc_list[i].direct
                    if direct & DIR_X:
                        self.set_boundary_condition(j, 0, val)
                    if direct & DIR_Y:
                        self.set_boundary_condition(j, 1, val)
                    if direct & DIR_Z:
                        self.set_boundary_condition(j, 2, val)

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
