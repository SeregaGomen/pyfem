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
from numpy import array


class TFEMStatic(TFEM):
    def __init__(self):
        super().__init__()
        self._global_matrix_stiffness = lil_matrix((0, 0))   # Глобальная матрица жесткости (ГМЖ)
        self._global_load = []                               # Глобальный вектор нагрузок (правая часть)

    # Расчет статической задачи методом конечных элементов
    def _calc_problem(self):
        # Создание ГМЖ
        size = len(self.mesh.x) * self.mesh.freedom
        self._global_matrix_stiffness = lil_matrix((size, size))
        self._global_load = [0] * size

        fe = self.create_fe()
        fe.set_params(self.params)
        # Вычисление компонент нагрузки
        self._prepare_concentrated_load()
        self._prepare_surface_load()
        self._prepare_volume_load()
        # Формирование глобальной матрицы жесткости
        self._progress.set_process('Assembling global stiffness matrix...', 1, len(self.mesh.fe))
        for i in range(0, len(self.mesh.fe)):
            self._progress.set_progress(i + 1)
            x = self.mesh.get_fe_coord(i)
            fe.set_coord(x)
            fe.generate()
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly(fe, i)
        # Учет краевых условий
        self._use_boundary_condition()
        # Решение СЛАУ
        if not self._solve():
            print('The system of equations is not solved!')
            return False
        self._calc_results()
        print('**************** Success! ****************')
        return True

    # Добавление локальной матрицы жесткости (ЛМЖ) к ГМЖ
    def __assembly(self, fe, index):
        # Добавление матрицы
        for i in range(0, len(fe.K)):
            k = self.mesh.fe[index][i // self.mesh.freedom] * self.mesh.freedom + i % self.mesh.freedom
            for j in range(i, len(fe.K)):
                r = self.mesh.fe[index][j // self.mesh.freedom] * self.mesh.freedom + j % self.mesh.freedom
                self._global_matrix_stiffness[k, r] += fe.K[i][j]
                if k != r:
                    self._global_matrix_stiffness[r, k] += fe.K[i][j]

    # Вычисление сосредоточенных нагрузок
    def _prepare_concentrated_load(self, t=0):
        counter = 0
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'concentrated':
                counter += 1
        if not counter:
            return
        self._progress.set_process('Computation of concentrated load...', 1, counter * len(self.mesh.x))
        counter = 1
        for i in range(0, len(self.params.bc_list)):
            if not self.params.bc_list[i].type == 'concentrated':
                continue
            for j in range(0, len(self.mesh.x)):
                self._progress.set_progress(counter)
                counter += 1
                parser = self.create_parser(self.mesh.get_coord(j), t)
                if len(self.params.bc_list[i].predicate):
                    parser.set_code(self.params.bc_list[i].predicate)
                    if parser.run() == 0:
                        continue
                parser.set_code(self.params.bc_list[i].expression)
                load = parser.run()
                if self.params.bc_list[i].direct & DIR_1:
                    self._global_load[j * self.mesh.freedom + 0] += load
                if self.params.bc_list[i].direct & DIR_2:
                    self._global_load[j * self.mesh.freedom + 1] += load
                if self.params.bc_list[i].direct & DIR_3:
                    self._global_load[j * self.mesh.freedom + 2] += load

    # Вычисление поверхностных нагрузок
    def _prepare_surface_load(self, t=0):
        if not len(self.mesh.be):
            return
        counter = 0
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'surface':
                counter += 1
        if not counter:
            return
        self._progress.set_process('Computation of surface load...', 1, counter * len(self.mesh.be))
        counter = 1
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type != 'surface':
                continue
            for j in range(0, len(self.mesh.be)):
                self._progress.set_progress(counter)
                counter += 1
                if not self.__check_be(j, self.params.bc_list[i].predicate):
                    continue
                share = self.__surface_load_share(j)
                for k in range(0, len(self.mesh.be[j])):
                    parser = self.create_parser(self.mesh.get_coord(self.mesh.be[j][k]), t)
                    parser.set_code(self.params.bc_list[i].expression)
                    load = parser.run()
                    if self.params.bc_list[i].direct & DIR_1:
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 0] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_2:
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 1] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_3:
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 2] += load * share[k]

    # Вычисление объемных нагрузок
    def _prepare_volume_load(self, t=0):
        counter = 0
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'volume':
                counter += 1
        if not counter:
            return
        self._progress.set_process('Computation of volume load...', 1, counter * len(self.mesh.fe))
        counter = 1
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type != 'volume':
                continue
            for j in range(0, len(self.mesh.fe)):
                self._progress.set_progress(counter)
                counter += 1
                if not self.__check_fe(j, self.params.bc_list[i].predicate):
                    continue
                share = self.__volume_load_share(j)
                for k in range(0, len(self.mesh.fe[j])):
                    parser = self.create_parser(self.mesh.get_coord(self.mesh.fe[j][k]), t)
                    if len(self.params.bc_list[i].predicate):
                        parser.set_code(self.params.bc_list[i].predicate)
                        if parser.run() == 0:
                            continue
                    parser.set_code(self.params.bc_list[i].expression)
                    load = parser.run()
                    if self.params.bc_list[i].direct & DIR_1:
                        self._global_load[self.mesh.fe[j][k] * self.mesh.freedom + 0] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_2:
                        self._global_load[self.mesh.fe[j][k] * self.mesh.freedom + 1] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_3:
                        self._global_load[self.mesh.fe[j][k] * self.mesh.freedom + 2] += load * share[k]

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...)
    def _calc_results(self, t=0):
        # Выделяем память для хранения результатов
        res = []
        for i in range(0, self.__num_result()):
            r = [0]*len(self.mesh.x)
            res.append(r)
        uvw = [0] * len(self.mesh.fe[0]) * self.mesh.freedom
        counter = [0]*len(self.mesh.x)  # Счетчик кол-ва вхождения узлов для осреднения результатов
        # Копируем полученные перемещения
        for i in range(0, len(self.mesh.x)):
            if self.mesh.is_plate():  # Для пластин перемещения меняются местами и вычисляются u и v
                res[0][i] = self._global_load[i * self.mesh.freedom + 1]   # u (Tau_x)
                res[1][i] = self._global_load[i * self.mesh.freedom + 2]   # v (Tau_y)
                res[2][i] = self._global_load[i * self.mesh.freedom + 0]   # w
            else:
                for j in range(0, self.mesh.freedom):
                    res[j][i] = self._global_load[i * self.mesh.freedom + j]
        # Вычисляем стандартные результаты по всем КЭ
        fe = self.create_fe()
        fe.set_params(self.params)
        self._progress.set_process('Calculation results...', 1, len(self.mesh.fe))
        for i in range(0, len(self.mesh.fe)):
            self._progress.set_progress(i + 1)
            x = self.mesh.get_fe_coord(i)
            fe.set_coord(x)
            for j in range(0, len(self.mesh.fe[i])):
                for k in range(0, self.mesh.freedom):
                    uvw[j * self.mesh.freedom + k] = \
                        self._global_load[self.mesh.freedom * self.mesh.fe[i][j] + k]
            r = fe.calc(uvw)
            offset = 3 if self.mesh.is_shell() else self.mesh.freedom
            for m in range(0, len(r)):
                for j in range(0, len(r[0])):
                    res[offset + m][self.mesh.fe[i][j]] += r[m][j]
                    if not m:
                        counter[self.mesh.fe[i][j]] += 1
        # Осредняем результаты
        for i in range(self.mesh.freedom, self.__num_result()):
            for j in range(0, len(self.mesh.x)):
                res[i][j] /= counter[j]
        # Сохраняем полученные результаты в списке
        for i in range(0, self.__num_result()):
            r = TResult()
            r.name = self.params.names[self.__index_result(i)]
            r.results = res[i]
            r.t = t
            self._result.append(r)

    # Задание граничных условий
    def _set_boundary_condition(self, i, j, val):
        r = i * self.mesh.freedom + j
        for k in self._global_matrix_stiffness[r].nonzero()[1]:
            if r != k:
                self._global_matrix_stiffness[r, k] = self._global_matrix_stiffness[k, r] = 0
        self._global_load[r] = val * self._global_matrix_stiffness[r, r]

    # Учет граничных условий
    def _use_boundary_condition(self):
        counter = 0
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'boundary':
                counter += 1
        self._progress.set_process('Use of boundary conditions...', 1, counter * len(self.mesh.x))
        counter = 1
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'boundary':
                for j in range(0, len(self.mesh.x)):
                    self._progress.set_progress(counter)
                    counter += 1
                    parser = self.create_parser(self.mesh.get_coord(j))
                    if len(self.params.bc_list[i].predicate):
                        parser.set_code(self.params.bc_list[i].predicate)
                        if parser.run() == 0:
                            continue
                    parser.set_code(self.params.bc_list[i].expression)
                    val = parser.run()
                    direct = self.params.bc_list[i].direct
                    if direct & DIR_1:
                        self._set_boundary_condition(j, 0, val)
                    if direct & DIR_2:
                        self._set_boundary_condition(j, 1, val)
                    if direct & DIR_3:
                        self._set_boundary_condition(j, 2, val)

    # Прямое решение СЛАУ
    def _solve_direct(self):
        self._progress.set_process('Solving of equation system...', 1, 1)
        self._global_matrix_stiffness = self._global_matrix_stiffness.tocsr()
        try:
            self._global_load = spsolve(self._global_matrix_stiffness, self._global_load)
        except ArpackError:
            return False
        self._progress.set_progress(1)
        return True

    # Приближенное решение СЛАУ
    def _solve_iterative(self):
        self._progress.set_process('Solving of equation system...', 1, 1)
        self._global_matrix_stiffness = self._global_matrix_stiffness.tocsr()
        self._global_load, info = bicgstab(self._global_matrix_stiffness, self._global_load,
                                           self._global_load, self.params.eps)
        self._progress.set_progress(1)
        return True if not info else False

    # Проверка соответствия граничного элемента предикату отбора (всех его вершин)
    def __check_be(self, i, predicate):
        if not len(predicate):
            return True
        for k in range(0, len(self.mesh.be[0])):
            parser = self.create_parser(self.mesh.get_coord(self.mesh.be[i][k]))
            parser.set_code(predicate)
            if parser.run() == 0:
                return False
        return True

    # Проверка соответствия элемента предикату отбора (всех его вершин)
    def __check_fe(self, i, predicate):
        if not len(predicate):
            return True
        for k in range(0, len(self.mesh.fe[0])):
            parser = self.create_parser(self.mesh.get_coord(self.mesh.fe[i][k]))
            parser.set_code(predicate)
            if parser.run() == 0:
                return False
        return True

    # Определение кол-ва результатов в зависимости от размерности задачи
    def __num_result(self):
        res = 0
        if self.mesh.freedom == 1:
            # u, Exx, Sxx
            res = 3
        elif self.mesh.freedom == 2:
            # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy
            res = 8
        elif self.mesh.freedom == 3 or self.mesh.freedom == 6:
            # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz
            res = 15
        return res

    # Индекс функции в зависимости от размерности задачи
    def __index_result(self, i):
        ret = 0
        # u, Exx, Sxx
        index1 = [4, 7, 13]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy
        index2 = [4, 5, 7, 8, 10, 13, 14, 16]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz
        index3 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        if self.mesh.freedom == 1:
            ret = index1[i]
        elif self.mesh.freedom == 2:
            ret = index2[i]
        elif self.mesh.freedom == 3 or self.mesh.freedom == 6:
            ret = index3[i]
        return ret

    # Определение компонент поверхностной нагрузки в зависимости от типа КЭ
    def __surface_load_share(self, index):
        share = array([])
        if self.mesh.fe_type == 'fe_1d_2':
            share = array([1 * self.params.thickness])
        elif self.mesh.fe_type == 'fe_2d_3' or self.mesh.fe_type == 'fe_2d_4':
            share = array([1 / 2, 1 / 2]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_2d_6':
            share = array([1 / 6, 1 / 6, 2 / 3]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_3d_4' or self.mesh.fe_type == 'fe_2d_3_p' or self.mesh.fe_type == 'fe_2d_3_s':
            share = array([1 / 3, 1 / 3, 1 / 3]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_3d_8' or self.mesh.fe_type == 'fe_2d_4_p' or self.mesh.fe_type == 'fe_2d_4_s':
            share = array([1, 1, 1, 1]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_3d_10':
            share = array([0, 0, 0, 1 / 6, 1 / 6, 1/ 6]) * self.mesh.square(index)
        return share

    # Определение компонент объемной нагрузки в зависимости от типа КЭ
    def __volume_load_share(self, index):
        share = array([])
        if self.mesh.fe_type == 'fe_1d_2':
            share = array([1 / 2]) * self.mesh.volume(index) * self.params.thickness
        elif self.mesh.fe_type == 'fe_2d_3' or self.mesh.fe_type == 'fe_2d_3_p' or self.mesh.fe_type == 'fe_2d_3_s':
            share = array([1 / 6, 1 / 6, 1 / 6]) * self.mesh.volume(index) * self.params.thickness
        elif self.mesh.fe_type == 'fe_2d_4' or self.mesh.fe_type == 'fe_2d_4_p' or self.mesh.fe_type == 'fe_2d_4_s':
            share = array([1, 1, 1, 1]) * self.mesh.volume(index) * self.params.thickness
        elif self.mesh.fe_type == 'fe_2d_6':
            share = array([0, 0, 0, 1 / 6, 1 / 6, 1 / 6]) * self.mesh.volume(index) * self.params.thickness
        elif self.mesh.fe_type == 'fe_3d_4':
            share = array([1 / 24, 1 / 24, 1 / 24, 1 / 24]) * self.mesh.volume(index)
        elif self.mesh.fe_type == 'fe_3d_8':
            share = array([1, 1, 1, 1, 1, 1, 1, 1]) * self.mesh.volume(index)
        elif self.mesh.fe_type == 'fe_3d_10':
            share = array([-1 / 120, -1 / 120, -1 / 120, -1 / 120, 1 / 30, 1 / 30, 1 / 30, 1 / 30, 1 / 30, 1 / 30]) * \
                    self.mesh.volume(index)
        return share
