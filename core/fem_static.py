#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve, bicgstab
from numpy import array
from core.fem_fem import TFEM
from core.fem_defs import DIR_X, DIR_Y, DIR_Z
from core.fem_result import TResult
from core.fem_progress import TThreadProgress


class TFEMStatic(TFEM):
    def __init__(self):
        super().__init__()
        self._global_matrix_stiffness = lil_matrix((0, 0))  # Глобальная матрица жесткости (ГМЖ)
        self._global_load = []                              # Глобальный вектор нагрузок (правая часть)
        self._fe_thickness = []                             # Толщина каждого КЭ

    # Расчет статической задачи методом конечных элементов
    def _calc_problem(self):
        # Создание ГМЖ
        size = len(self.mesh.x) * self.mesh.freedom
        self._global_matrix_stiffness = lil_matrix((size, size))
        self._global_load = [0] * size

        # Формирование глобальной матрицы жесткости
        self._progress.set_process('Assembling global stiffness matrix...', 1, len(self.mesh.fe))
        fe = self.create_fe()
        for i in range(len(self.mesh.fe)):
            self._progress.set_progress(i + 1)
            # Настройка КЭ
            self._set_fe(fe, i)
            fe.generate()
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly(fe, i)
        # Учет состредоточенных, поверхностных и объемных нагрузок
        self._use_surface_load()
        self._use_volume_load()
        self._use_concentrated_load()
        self._use_pressure_load()
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
            # Учет нагрузки
            self._global_load[k] += fe.load[i]

    # Вычисление сосредоточенных нагрузок
    def _use_concentrated_load(self, t=0):
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
                x = self.mesh.get_coord(j)
                if t != 0:
                    x.append(t)
                if self.params.bc_list[i].predicate is not None:
                    if not self.params.bc_list[i].predicate(x):
                        continue
                load = self.params.bc_list[i].expression(x)
                if self.params.bc_list[i].direct & DIR_X:
                    self._global_load[j * self.mesh.freedom + 0] += load
                if self.params.bc_list[i].direct & DIR_Y:
                    self._global_load[j * self.mesh.freedom + 1] += load
                if self.params.bc_list[i].direct & DIR_Z:
                    self._global_load[j * self.mesh.freedom + 2] += load

    # Вычисление поверхностных нагрузок
    def _use_surface_load(self, t=0):
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
                x = self.mesh.get_be_center(j)
                if t != 0:
                    x.append(t)
                load = self.params.bc_list[i].expression(x)
                for k in range(0, len(self.mesh.be[j])):
                    if self.params.bc_list[i].direct & DIR_X:
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 0] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_Y:
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 1] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_Z:
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 2] += load * share[k]

    # Вычисление объемных нагрузок
    def _use_volume_load(self, t=0):
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
                x = self.mesh.get_fe_center(j)
                if t != 0:
                    x.append(t)
                share = self.__volume_load_share(j)
                load = self.params.bc_list[i].expression(x)
                for k in range(0, len(self.mesh.fe[j])):
                    if self.params.bc_list[i].direct & DIR_X:
                        self._global_load[self.mesh.fe[j][k] * self.mesh.freedom + 0] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_Y:
                        self._global_load[self.mesh.fe[j][k] * self.mesh.freedom + 1] += load * share[k]
                    if self.params.bc_list[i].direct & DIR_Z:
                        self._global_load[self.mesh.fe[j][k] * self.mesh.freedom + 2] += load * share[k]

    # Вычисление нагрузок поверхностного давления
    def _use_pressure_load(self, t=0):
        counter = 0
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'pressure':
                counter += 1
        if not counter:
            return
        self._progress.set_process('Computation of pressure load...', 1, counter * len(self.mesh.be))
        counter = 1
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type != 'pressure':
                continue
            for j in range(len(self.mesh.be)):
                self._progress.set_progress(counter)
                counter += 1
                if not self.__check_be(j, self.params.bc_list[i].predicate):
                    continue

                x = self.mesh.get_be_center(j)
                share = self.__surface_load_share(j)
                load = self.params.bc_list[i].expression(x)
                # Вычисление нормали к ГЭ
                v = self.mesh.be_normal(j)
                for k in range(0, len(self.mesh.be[j])):
                    if self.mesh.is_plate():
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 0] += load * share[k]
                    else:
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 0] += load * share[k] * v[0]
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 1] += load * share[k] * v[1]
                        self._global_load[self.mesh.be[j][k] * self.mesh.freedom + 2] += load * share[k] * v[2]

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
        self._progress.set_process('Calculation results...', 1, len(self.mesh.fe))
        for i in range(len(self.mesh.fe)):
            self._progress.set_progress(i + 1)
            self._set_fe(fe, i)
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
            self.results.append(r)

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
                    x = self.mesh.get_coord(j)
                    if self.params.bc_list[i].predicate is not None:
                        val = self.params.bc_list[i].predicate(x)
                        if not val:
                            continue

                    val = self.params.bc_list[i].expression(x)
                    direct = self.params.bc_list[i].direct
                    if direct & DIR_X:
                        self._set_boundary_condition(j, 0, val)
                    if direct & DIR_Y:
                        self._set_boundary_condition(j, 1, val)
                    if direct & DIR_Z:
                        self._set_boundary_condition(j, 2, val)

    # Прямое решение СЛАУ
    def _solve_direct(self):
        # self._progress.set_process('Solving of equation system...', 1, 1)
        progress = TThreadProgress('Solving of equation system...')
        progress.start()
        self._global_matrix_stiffness = self._global_matrix_stiffness.tocsr()
        try:
            self._global_load = spsolve(self._global_matrix_stiffness, self._global_load)
            progress.stop()
        except ValueError:
            progress.stop()
            return False
        except RuntimeError:
            progress.stop()
            return False
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
        if predicate is None:
            return True
        for k in range(self.mesh.base_be_size()):
            x = self.mesh.get_coord(self.mesh.be[i][k])
            if not predicate(x):
                return False
        return True

    # Проверка соответствия элемента предикату отбора (всех его вершин)
    def __check_fe(self, i, predicate):
        if predicate is None:
            return True
        for k in range(self.mesh.base_fe_size()):
            x = self.mesh.get_coord(self.mesh.fe[i][k])
            if not predicate(x):
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
            share = array([1]) * 1 if not len(self._fe_thickness) else self._fe_thickness[index]
        elif self.mesh.fe_type == 'fe_2d_3' or self.mesh.fe_type == 'fe_2d_4':
            share = array([1 / 2, 1 / 2]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_2d_6':
            share = array([1 / 6, 1 / 6, 2 / 3]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_3d_4' or self.mesh.fe_type == 'fe_2d_3_p' or self.mesh.fe_type == 'fe_3d_3_s':
            share = array([1 / 3, 1 / 3, 1 / 3]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_3d_8' or self.mesh.fe_type == 'fe_2d_4_p' or self.mesh.fe_type == 'fe_3d_4_s':
            share = array([1 / 4, 1 / 4, 1 / 4, 1 / 4]) * self.mesh.square(index)
        elif self.mesh.fe_type == 'fe_3d_10' or self.mesh.fe_type == 'fe_2d_6_p' or self.mesh.fe_type == 'fe_3d_6_s':
            share = array([0, 0, 0, 1 / 3, 1 / 3, 1 / 3]) * self.mesh.square(index)
        return share

    # Определение компонент объемной нагрузки в зависимости от типа КЭ
    def __volume_load_share(self, index):
        thickness = 1 if not len(self._fe_thickness) else self._fe_thickness[index]
        share = array([])
        if self.mesh.fe_type == 'fe_1d_2':
            share = array([1 / 2, 1 / 2]) * self.mesh.volume(index) * thickness
        elif self.mesh.fe_type == 'fe_2d_3' or self.mesh.fe_type == 'fe_2d_3_p' or self.mesh.fe_type == 'fe_3d_3_s':
            share = array([1 / 3, 1 / 3, 1 / 3]) * self.mesh.volume(index) * thickness
        elif self.mesh.fe_type == 'fe_2d_4' or self.mesh.fe_type == 'fe_2d_4_p' or self.mesh.fe_type == 'fe_3d_4_s':
            share = array([1 / 4, 1 / 4, 1 / 4, 1 / 4]) * self.mesh.volume(index) * thickness
        elif self.mesh.fe_type == 'fe_2d_6' or self.mesh.fe_type == 'fe_2d_6_p' or self.mesh.fe_type == 'fe_3d_6_s':
            share = array([0, 0, 0, 1 / 3, 1 / 3, 1 / 3]) * self.mesh.volume(index) * thickness
        elif self.mesh.fe_type == 'fe_3d_4':
            share = array([1 / 4, 1 / 4, 1 / 4, 1 / 4]) * self.mesh.volume(index)
        elif self.mesh.fe_type == 'fe_3d_8':
            share = array([1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8]) * self.mesh.volume(index)
        elif self.mesh.fe_type == 'fe_3d_10':
            share = array([-1 / 20, -1 / 20, -1 / 20, -1 / 20, 1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5, 1 / 5]) * \
                    self.mesh.volume(index)
        return share

    # Определение координат центра КЭ
    def _fe_center(self, x):
        cx = array([0., 0., 0.])
        n = self.mesh.base_fe_size()
        for i in range(n):
            for j in range(len(x[i])):
                cx[j] += x[i][j]
        return cx / n

    # Настройка параметров КЭ
    def _set_fe(self, fe, fe_index):
        coord = self.mesh.get_fe_coord(fe_index)
        fe.set_coord(coord)
        x = self._fe_center(coord)
        # Определение переменных параметров КЭ
        for i in range(len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'thickness' or self.params.bc_list[i].type == 'young_modulus' or \
                    self.params.bc_list[i].type == 'poisson_ratio' or self.params.bc_list[i].type == 'temperature' or \
                    self.params.bc_list[i].type == 'alpha' or self.params.bc_list[i].type == 'density' or \
                    self.params.bc_list[i].type == 'damping' or self.params.bc_list[i].type == 'shear_modulus':
                param = self.params.bc_list[i].type
            else:
                continue

            if self.params.bc_list[i].predicate is not None:
                val = self.params.bc_list[i].predicate(x)
                if not val:
                    continue

            val = self.params.bc_list[i].expression(x)
            if param == 'thickness':
                fe.set_thickness(val)
                self._fe_thickness.append(val)
            elif param == 'young_modulus':
                fe.set_young_modulus(val)
            elif param == 'shear_modulus':
                fe.set_shear_modulus(val)
            elif param == 'poisson_ratio':
                fe.set_poisson_ratio(val)
            elif param == 'temperature':
                fe.set_temperature(val)
            elif param == 'alpha':
                fe.set_alpha(val)
            elif param == 'density':
                fe.set_density(val)
            elif param == 'damping':
                fe.set_damping(val)
