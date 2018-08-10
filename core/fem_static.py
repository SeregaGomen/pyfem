#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve, bicgstab, ArpackError
from core.fem_fem import TFEM
from core.fem_defs import DIR_1, DIR_2, DIR_3
from core.fem_result import TResult


class TFEMStatic(TFEM):
    def __init__(self):
        super().__init__()
        self.__global_matrix_stiffness__ = lil_matrix((0, 0))   # Глобальная матрица жесткости (ГМЖ)
        self.__global_load__ = []                               # Глобальный вектор нагрузок (правая часть)

    # Расчет статической задачи методом конечных элементов
    def __calc_problem__(self):
        # Создание ГМЖ
        size = len(self.__mesh__.x)*self.__mesh__.freedom
        self.__global_matrix_stiffness__ = lil_matrix((size, size))
        self.__global_load__ = [0]*size

        fe = self.__create_fe__()
        fe.set_elasticity(self.__params__.e, self.__params__.m)
        fe.set_h(self.__params__.h)
        # Вычисление компонент нагрузки
        self.__prepare_concentrated_load__()
        self.__prepare_surface_load__()
        self.__prepare_volume_load__()
        # Формирование глобальной матрицы жесткости
        self.__progress__.set_process('Assembling global stiffness matrix...', 1, len(self.__mesh__.fe))
        for i in range(0, len(self.__mesh__.fe)):
            self.__progress__.set_progress(i + 1)
            x = self.__mesh__.get_fe_coord(i)
            fe.set_coord(x)
            fe.generate()
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly__(fe, i)
        # Учет краевых условий
        self.__use_boundary_condition__()
        # Решение СЛАУ
        if not self.__solve__():
            print('The system of equations is not solved!')
            return False
        self.__calc_results__()
        print('**************** Success! ****************')
        return True

    # Добавление локальной матрицы жесткости (ЛМЖ) к ГМЖ
    def __assembly__(self, fe, index):
        # Добавление матрицы
        for i in range(0, len(fe.K)):
            k = self.__mesh__.fe[index][i//self.__mesh__.freedom]*self.__mesh__.freedom + i % self.__mesh__.freedom
            for j in range(i, len(fe.K)):
                r = self.__mesh__.fe[index][j//self.__mesh__.freedom]*self.__mesh__.freedom + j % self.__mesh__.freedom
                self.__global_matrix_stiffness__[k, r] += fe.K[i][j]
                if k != r:
                    self.__global_matrix_stiffness__[r, k] += fe.K[i][j]

    # Вычисление сосредоточенных нагрузок
    def __prepare_concentrated_load__(self, t=0):
        parser = self.__create_parser__()
        counter = 0
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'concentrated':
                counter += 1
        if not counter:
            return
        self.__progress__.set_process('Computation of concentrated load...', 1, counter*len(self.__mesh__.x))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if not self.__params__.bc_list[i].type == 'concentrated':
                continue
            for j in range(0, len(self.__mesh__.x)):
                self.__progress__.set_progress(counter)
                counter += 1
                x = self.__mesh__.get_coord(j)
                for k in range(0, len(x)):
                    parser.set_variable(self.__params__.names[k], x[k])
                parser.set_variable(self.__params__.names[3], t)
                if len(self.__params__.bc_list[i].predicate):
                    parser.set_code(self.__params__.bc_list[i].predicate)
                    if parser.run() == 0:
                        continue
                parser.set_code(self.__params__.bc_list[i].expression)
                val = parser.run()
                if self.__params__.bc_list[i].direct & DIR_1:
                    self.__global_load__[j*self.__mesh__.freedom + 0] += val
                if self.__params__.bc_list[i].direct & DIR_2:
                    self.__global_load__[j*self.__mesh__.freedom + 1] += val
                if self.__params__.bc_list[i].direct & DIR_3:
                    self.__global_load__[j*self.__mesh__.freedom + 2] += val

    # Вычисление поверхностных нагрузок
    def __prepare_surface_load__(self, t=0):
        if not len(self.__mesh__.be):
            return
        parser = self.__create_parser__()
        counter = 0
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'surface':
                counter += 1
        if not counter:
            return
        self.__progress__.set_process('Computation of surface load...', 1, counter*len(self.__mesh__.be))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type != 'surface':
                continue
            for j in range(0, len(self.__mesh__.be)):
                self.__progress__.set_progress(counter)
                counter += 1
                if not self.__check_boundary_elements__(j, self.__params__.bc_list[i].predicate):
                    continue
                rel_se = self.__mesh__.square(j)/float(len(self.__mesh__.be[j]))
                for k in range(0, len(self.__mesh__.be[j])):
                    x = self.__mesh__.get_coord(self.__mesh__.be[j][k])
                    for l in range(0, len(x)):
                        parser.set_variable(self.__params__.names[l], x[l])
                    parser.set_variable(self.__params__.names[3], t)
                    parser.set_code(self.__params__.bc_list[i].expression)
                    val = parser.run()
                    if self.__params__.bc_list[i].direct & DIR_1:
                        self.__global_load__[self.__mesh__.be[j][k]*self.__mesh__.freedom + 0] += val*rel_se
                    if self.__params__.bc_list[i].direct & DIR_2:
                        self.__global_load__[self.__mesh__.be[j][k]*self.__mesh__.freedom + 1] += val*rel_se
                    if self.__params__.bc_list[i].direct & DIR_3:
                        self.__global_load__[self.__mesh__.be[j][k]*self.__mesh__.freedom + 2] += val*rel_se

    # Вычисление объемных нагрузок
    def __prepare_volume_load__(self, t=0):
        parser = self.__create_parser__()
        counter = 0
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'volume':
                counter += 1
        if not counter:
            return
        self.__progress__.set_process('Computation of volume load...', 1, counter*len(self.__mesh__.fe))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type != 'volume':
                continue
            for j in range(0, len(self.__mesh__.fe)):
                self.__progress__.set_progress(counter)
                counter += 1
                rel_ve = self.__mesh__.volume(j)/float(len(self.__mesh__.fe[j]))
                for k in range(0, len(self.__mesh__.fe[j])):
                    x = self.__mesh__.get_coord(self.__mesh__.fe[j][k])
                    for l in range(0, len(x)):
                        parser.set_variable(self.__params__.names[l], x[l])
                    parser.set_variable(self.__params__.names[3], t)
                    parser.set_code(self.__params__.bc_list[i].expression)
                    val = parser.run()
                    if self.__params__.bc_list[i].direct & DIR_1:
                        self.__global_load__[self.__mesh__.fe[j][k]*self.__mesh__.freedom + 0] += val*rel_ve
                    if self.__params__.bc_list[i].direct & DIR_2:
                        self.__global_load__[self.__mesh__.fe[j][k]*self.__mesh__.freedom + 1] += val*rel_ve
                    if self.__params__.bc_list[i].direct & DIR_3:
                        self.__global_load__[self.__mesh__.fe[j][k]*self.__mesh__.freedom + 2] += val*rel_ve

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...)
    def __calc_results__(self, t=0):
        # Выделяем память для хранения результатов
        res = []
        for i in range(0, self.__num_result__()):
            r = [0]*len(self.__mesh__.x)
            res.append(r)
        uvw = [0]*len(self.__mesh__.fe[0])*self.__mesh__.freedom
        counter = [0]*len(self.__mesh__.x)  # Счетчик кол-ва вхождения узлов для осреднения результатов
        # Копируем полученные перемещения
        for i in range(0, len(self.__mesh__.x)):
            for j in range(0, self.__mesh__.freedom):
                if self.__mesh__.is_plate():  # Для пластин перемещения меняются местами и вычисляются u и v
                    if j < 2:
                        res[j][i] = -self.__global_load__[i * self.__mesh__.freedom + j + 1] * \
                                    self.__global_load__[i * self.__mesh__.freedom + 0]
                    else:
                        res[j][i] = self.__global_load__[i * self.__mesh__.freedom + 0]
                else:
                    res[j][i] = self.__global_load__[i * self.__mesh__.freedom + j]
        # Вычисляем стандартные результаты по всем КЭ
        fe = self.__create_fe__()
        fe.set_elasticity(self.__params__.e, self.__params__.m)
        self.__progress__.set_process('Calculation results...', 1, len(self.__mesh__.fe))
        for i in range(0, len(self.__mesh__.fe)):
            self.__progress__.set_progress(i + 1)
            x = self.__mesh__.get_fe_coord(i)
            fe.set_coord(x)
            for j in range(0, len(self.__mesh__.fe[i])):
                for k in range(0, self.__mesh__.freedom):
                    uvw[j * self.__mesh__.freedom + k] = \
                        self.__global_load__[self.__mesh__.freedom*self.__mesh__.fe[i][j] + k]
            r = fe.calc(uvw)
            for m in range(0, len(r)):
                for j in range(0, len(r[0])):
                    res[self.__mesh__.freedom + m][self.__mesh__.fe[i][j]] += r[m][j]
                    if not m:
                        counter[self.__mesh__.fe[i][j]] += 1
        # Осредняем результаты
        for i in range(self.__mesh__.freedom, self.__num_result__()):
            for j in range(0, len(self.__mesh__.x)):
                res[i][j] /= counter[j]
        # Сохраняем полученные результаты в списке
        for i in range(0, self.__num_result__()):
            r = TResult()
            r.name = self.__params__.names[self.__index_result__(i)]
            r.results = res[i]
            r.t = t
            self.__result__.append(r)

    # Задание граничных условий
    def __set_boundary_condition__(self, i, j, val):
        r = i*self.__mesh__.freedom + j
        for k in self.__global_matrix_stiffness__[r].nonzero()[1]:
            if r != k:
                self.__global_matrix_stiffness__[r, k] = self.__global_matrix_stiffness__[k, r] = 0
        self.__global_load__[r] = val*self.__global_matrix_stiffness__[r, r]

    # Учет граничных условий
    def __use_boundary_condition__(self):
        parser = self.__create_parser__()
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
                    x = self.__mesh__.get_coord(j)
                    for k in range(0, len(x)):
                        parser.set_variable(self.__params__.names[k], x[k])
                    if len(self.__params__.bc_list[i].predicate):
                        parser.set_code(self.__params__.bc_list[i].predicate)
                        if parser.run() == 0:
                            continue
                    parser.set_code(self.__params__.bc_list[i].expression)
                    val = parser.run()
                    direct = self.__params__.bc_list[i].direct
                    if direct & DIR_1:
                        self.__set_boundary_condition__(j, 0, val)
                    if direct & DIR_2:
                        self.__set_boundary_condition__(j, 1, val)
                    if direct & DIR_3:
                        self.__set_boundary_condition__(j, 2, val)

    # Прямое решение СЛАУ
    def __solve_direct__(self):
        self.__progress__.set_process('Solving of equation system...', 1, 1)
        self.__global_matrix_stiffness__ = self.__global_matrix_stiffness__.tocsr()
        try:
            self.__global_load__ = spsolve(self.__global_matrix_stiffness__, self.__global_load__)
        except ArpackError:
            return False
        self.__progress__.set_progress(1)
        return True

    # Приближенное решение СЛАУ
    def __solve_iterative__(self):
        self.__progress__.set_process('Solving of equation system...', 1, 1)
        self.__global_matrix_stiffness__ = self.__global_matrix_stiffness__.tocsr()
        self.__global_load__, info = bicgstab(self.__global_matrix_stiffness__, self.__global_load__,
                                              self.__global_load__, self.__params__.eps)
        self.__progress__.set_progress(1)
        return True if not info else False

    # Проверка соответствия граничного элемента предикату отбора (всех его вершин)
    def __check_boundary_elements__(self, i, predicate):
        if not len(predicate):
            return True
        parser = self.__create_parser__()
        for k in range(0, len(self.__mesh__.be[0])):
            x = self.__mesh__.get_coord(self.__mesh__.be[i][k])
            for l in range(0, len(x)):
                parser.set_variable(self.__params__.names[l], x[l])
            parser.set_code(predicate)
            if parser.run() == 0:
                return False
        return True

    # Определение кол-ва результатов в зависимости от размерности задачи
    def __num_result__(self):
        res = 0
        if self.__mesh__.freedom == 1:
            # u, Exx, Sxx
            res = 3
        elif self.__mesh__.freedom == 2:
            # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy
            res = 8
        elif self.__mesh__.freedom == 3:
            if self.__mesh__.is_plate():
                # u, v, w, Exx, Eyy, Exy, Sxx, Syy, Sxy
                res = 9
            else:
                # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz
                res = 15
        elif self.__mesh__.freedom == 6:
            # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz
            res = 15
        return res

    # Индекс функции в зависимости от размерности задачи
    def __index_result__(self, i):
        ret = 0
        # u, Exx, Sxx
        index1 = [4, 7, 13]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy
        index2 = [4, 5, 7, 8, 10, 13, 14, 16]
        # u, v, w, Exx, Eyy, Exy, Sxx, Syy, Sxy
        index3 = [4, 5, 6, 7, 8, 10, 13, 14, 16]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz
        index4 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        if self.__mesh__.freedom == 1:
            ret = index1[i]
        elif self.__mesh__.freedom == 2:
            ret = index2[i]
        elif self.__mesh__.freedom == 3:
            if self.__mesh__.is_plate():
                ret = index3[i]
            else:
                ret = index4[i]
        elif self.__mesh__.freedom == 6:
            ret = index4[i]
        return ret
