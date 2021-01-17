#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

import math
from scipy.sparse import lil_matrix, coo_matrix
from numpy import zeros, savez, load
from core.fem_defs import INIT_U, INIT_V, INIT_W, INIT_UT, INIT_VT, INIT_WT, INIT_UTT, INIT_VTT, INIT_WTT
from core.fem_static import TFEMStatic
from core.fem_parser import TParser


# Сохранение разреженной матрицы в файл
def save_matrix(file_name, matrix):
    matrix_coo = matrix.tocoo()
    row = matrix_coo.row
    col = matrix_coo.col
    data = matrix_coo.data
    shape = matrix_coo.shape
    savez(file_name, row=row, col=col, data=data, shape=shape)


# Загрузка разреженной матрицы из файла
def load_matrix(file_name):
    matrix = load(file_name)
    return coo_matrix((matrix['data'], (matrix['row'], matrix['col'])), shape=matrix['shape']).tolil()


class TFEMDynamic(TFEMStatic):
    def __init__(self):
        super().__init__()
        self._global_matrix_mass = lil_matrix((0, 0))        # Глобальная матрица масс (ГММ)
        self._global_matrix_damping = lil_matrix((0, 0))     # Глобальная матрица демпфирования (ГМД)

    # Расчет динамической задачи методом конечных элементов
    def _calc_problem(self):
        size = len(self.mesh.x) * self.mesh.freedom
        self._global_matrix_stiffness = lil_matrix((size, size))
        self._global_matrix_mass = lil_matrix((size, size))
        self._global_matrix_damping = lil_matrix((size, size))
        self._global_load = [0] * size
        fe = self.create_fe()
        # Создание глобальных матриц жесткости, масс и демпфирования
        self._progress.set_process('Assembling global stiffness, mass and damping matrix...', 1, len(self.mesh.fe))
        for i in range(0, len(self.mesh.fe)):
            self._progress.set_progress(i + 1)
            # Настройка КЭ
            self._set_fe(fe, i)
            fe.generate(False)
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly(fe, i)
        # Формирование левой части СЛАУ
        self.__create_dynamic_matrix()
        # Учет начальных условий
        u0, ut0, utt0 = self.__prepare_initial_condition()
        # Итерационный процесс по времени
        t = self.params.t0
        while t <= self.params.t1:
            print('t = %5.2f' % t)
            # Формирование правой части СЛАУ
            self.__create_dynamic_vector(u0, ut0, utt0, t)
            # Учет краевых условий
            self._use_boundary_condition()
            # Решение СЛАУ
            if not self._solve():
                print('The system of equations is not solved!')
                return False
            u0, ut0, utt0 = self.__calc_dynamic_results(u0, ut0, utt0, t)
            t += self.params.th
            if math.fabs(t - self.params.t1) < self.params.eps:
                t = self.params.t1
        print('**************** Success! ****************')
        return True

    # Извлечение начальных условий
    def __prepare_initial_condition(self):
        u0 = zeros(len(self.mesh.x) * self.mesh.freedom)
        ut0 = zeros(len(self.mesh.x) * self.mesh.freedom)
        utt0 = zeros(len(self.mesh.x) * self.mesh.freedom)
        parser = TParser()
        counter = 0
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'initial':
                counter += 1
        self._progress.set_process('Using initial conditions...', 1, counter * len(self.mesh.x))
        counter = 1
        for i in range(0, len(self.params.bc_list)):
            if self.params.bc_list[i].type == 'initial':
                parser.set_code(self.params.bc_list[i].expression)
                value = parser.run()
                direct = self.params.bc_list[i].direct
                for j in range(0, len(self.mesh.x)):
                    self._progress.set_progress(counter)
                    counter += 1
                    if direct & INIT_U:
                        u0[j * self.mesh.freedom + 0] = value
                    if direct & INIT_V:
                        u0[j * self.mesh.freedom + 1] = value
                    if direct & INIT_W:
                        u0[j * self.mesh.freedom + 2] = value
                    if direct & INIT_UT:
                        ut0[j * self.mesh.freedom + 0] = value
                    if direct & INIT_VT:
                        ut0[j * self.mesh.freedom + 1] = value
                    if direct & INIT_WT:
                        ut0[j * self.mesh.freedom + 2] = value
                    if direct & INIT_UTT:
                        utt0[j * self.mesh.freedom + 0] = value
                    if direct & INIT_VTT:
                        utt0[j * self.mesh.freedom + 1] = value
                    if direct & INIT_WTT:
                        utt0[j * self.mesh.freedom + 2] = value
        return u0, ut0, utt0

    # Вычисление напряжений, деформаций, скоростей и ускорений
    def __calc_dynamic_results(self, u0, ut0, utt0, t):
        # Вычисление деформаций и напряжений
        super()._calc_results(t)
        # Вычисление скоростей и ускорений (конечными разностями)
        th = self.params.th
        for i in range(0, len(self.mesh.x)):
            for j in range(0, self.mesh.freedom):
                u_0 = u0[i * self.mesh.freedom + j]
                u_1 = self._global_load[i * self.mesh.freedom + j]
                u_t_0 = ut0[i * self.mesh.freedom + j]
                u_t_1 = (u_1 - u_0) / th
                u_t_t_1 = (u_t_1 - u_t_0) / th
                u0[i * self.mesh.freedom + j] = u_1
                ut0[i * self.mesh.freedom + j] = u_t_1
                utt0[i * self.mesh.freedom + j] = u_t_t_1
                self.results[len(self.results) - 2 * self.mesh.freedom + j].results[i] = u_t_1
                self.results[len(self.results) - self.mesh.freedom + j].results[i] = u_t_t_1
        return u0, ut0, utt0

    # Добавление ЛМЖ, ЛММ и ЛМД к ГМЖ
    def __assembly(self, fe, index):
        # Добавление матрицы
        for i in range(0, len(fe.K)):
            k = self.mesh.fe[index][i // self.mesh.freedom] * self.mesh.freedom + i % self.mesh.freedom
            for j in range(i, len(fe.K)):
                r = self.mesh.fe[index][j // self.mesh.freedom] * self.mesh.freedom + j % self.mesh.freedom
                self._global_matrix_stiffness[k, r] += fe.K[i][j]
                self._global_matrix_mass[k, r] += fe.M[i][j]
                self._global_matrix_damping[k, r] += fe.C[i][j]
                if k != r:
                    self._global_matrix_stiffness[r, k] += fe.K[i][j]
                    self._global_matrix_mass[r, k] += fe.M[i][j]
                    self._global_matrix_damping[r, k] += fe.C[i][j]

    # Формирование левой части (матрицы) уравнения квазистатического равновесия
    def __create_dynamic_matrix(self):
        theta = 1.37
        self._progress.set_process('Creating static part of global matrix...', 1, 2)
        self._global_matrix_stiffness += \
            self._global_matrix_mass.dot(6.0 / (theta ** 2 ** self.params.th ** 2)) + \
            self._global_matrix_damping.dot(6.0 / (3.0 / (theta * self.params.th)))
        self._progress.set_progress(2)

    # Формирование правой части (вектора) уравнения квазистатического равновесия
    def __create_dynamic_vector(self, u0, ut0, utt0, t):
        theta = 1.37
        k1 = 6.0 / (theta ** 2 ** self.params.th ** 2)
        k2 = 3.0 / (theta * self.params.th)
        k3 = 0.5 * (theta * self.params.th)
        self._global_load = zeros(len(self.mesh.x) * self.mesh.freedom)
        # Вычисление компонент нагрузки для текущего момента времени
        self._use_concentrated_load(t)
        self._use_surface_load(t)
        self._use_volume_load(t)
        self._global_load += \
            self._global_matrix_mass.dot(u0.dot(k1) + ut0.dot(2.0 * k2) + utt0.dot(2.0)) + \
            self._global_matrix_damping.dot(u0.dot(k2) + ut0.dot(2.0) + utt0.dot(k3))

    # Определение кол-ва результатов в зависимости от размерности задачи
    def __num_result(self):
        res = 0
        if self.mesh.freedom == 1:
            # u, Exx, Sxx, ut, utt
            res = 5
        elif self.mesh.freedom == 2:
            # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
            res = 12
        elif self.mesh.freedom == 3:
            # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
            res = 21
        return res

    # Индекс функции в зависимости от размерности задачи
    def __index_result(self, i):
        ret = 0
        # u, Exx, Sxx, ut, utt
        index4 = [4, 7, 13, 19, 22]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
        index5 = [4, 5, 7, 8, 10, 13, 14, 16, 19, 20, 22, 23]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
        index6 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        if self.mesh.freedom == 1:
            ret = index4[i]
        elif self.mesh.freedom == 2:
            ret = index5[i]
        elif self.mesh.freedom == 3:
            ret = index6[i]
        return ret
