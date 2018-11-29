#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

import math
from scipy.sparse import lil_matrix, coo_matrix
from numpy import zeros, savez, load
from core.fem_defs import INIT_U, INIT_V, INIT_W, INIT_U_T, INIT_V_T, INIT_W_T, INIT_U_T_T, INIT_V_T_T, INIT_W_T_T
from core.fem_static import TFEMStatic


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
        self.__global_matrix_mass__ = lil_matrix((0, 0))        # Глобальная матрица масс (ГММ)
        self.__global_matrix_damping__ = lil_matrix((0, 0))     # Глобальная матрица демпфирования (ГМД)

    # Расчет динамической задачи методом конечных элементов
    def _calc_problem_(self):
        size = len(self._mesh_.x) * self._mesh_.freedom
        self.__global_matrix_stiffness__ = lil_matrix((size, size))
        self.__global_matrix_mass__ = lil_matrix((size, size))
        self.__global_matrix_damping__ = lil_matrix((size, size))
        self.__global_load__ = [0]*size
        fe = self._create_fe_()
        fe.set_elasticity(self._params_.e, self._params_.m)
        fe.set_damping(self._params_.damping)
        fe.set_density(self._params_.density)

        # Создание глобальных матриц жесткости, масс и демпфирования
        self._progress_.set_process('Assembling global stiffness, mass and damping matrix...', 1,
                                    len(self._mesh_.fe))
        for i in range(0, len(self._mesh_.fe)):
            self._progress_.set_progress(i + 1)
            # Настройка КЭ
            x = self._mesh_.get_fe_coord(i)
            fe.set_coord(x)
            fe.generate(False)
            # Ансамблирование ЛМЖ к ГМЖ
            self._assembly_(fe, i)
        # Формирование левой части СЛАУ
        self.__create_dynamic_matrix__()
        # Учет начальных условий
        u0, ut0, utt0 = self.__prepare_initial_condition__()
        # Итерационный процесс по времени
        t = self._params_.t0
        while t <= self._params_.t1:
            print('t = %5.2f' % t)
            # Формирование правой части СЛАУ
            self.__create_dynamic_vector__(u0, ut0, utt0, t)
            # Учет краевых условий
            self._use_boundary_condition_()
            # Решение СЛАУ
            if not self._solve_():
                print('The system of equations is not solved!')
                return False
            u0, ut0, utt0 = self.__calc_dynamic_results__(u0, ut0, utt0, t)
            t += self._params_.th
            if math.fabs(t - self._params_.t1) < self._params_.eps:
                t = self._params_.t1
        print('**************** Success! ****************')
        return True

    # Извлечение начальных условий
    def __prepare_initial_condition__(self):
        u0 = zeros(len(self._mesh_.x) * self._mesh_.freedom)
        ut0 = zeros(len(self._mesh_.x) * self._mesh_.freedom)
        utt0 = zeros(len(self._mesh_.x) * self._mesh_.freedom)
        parser = self._create_parser_()
        counter = 0
        for i in range(0, len(self._params_.bc_list)):
            if self._params_.bc_list[i].type == 'initial':
                counter += 1
        self._progress_.set_process('Using initial conditions...', 1, counter * len(self._mesh_.x))
        counter = 1
        for i in range(0, len(self._params_.bc_list)):
            if self._params_.bc_list[i].type == 'initial':
                parser.set_code(self._params_.bc_list[i].expression)
                value = parser.run()
                direct = self._params_.bc_list[i].direct
                for j in range(0, len(self._mesh_.x)):
                    self._progress_.set_progress(counter)
                    counter += 1
                    if direct & INIT_U:
                        u0[j * self._mesh_.freedom + 0] = value
                    if direct & INIT_V:
                        u0[j * self._mesh_.freedom + 1] = value
                    if direct & INIT_W:
                        u0[j * self._mesh_.freedom + 2] = value
                    if direct & INIT_U_T:
                        ut0[j * self._mesh_.freedom + 0] = value
                    if direct & INIT_V_T:
                        ut0[j * self._mesh_.freedom + 1] = value
                    if direct & INIT_W_T:
                        ut0[j * self._mesh_.freedom + 2] = value
                    if direct & INIT_U_T_T:
                        utt0[j * self._mesh_.freedom + 0] = value
                    if direct & INIT_V_T_T:
                        utt0[j * self._mesh_.freedom + 1] = value
                    if direct & INIT_W_T_T:
                        utt0[j * self._mesh_.freedom + 2] = value
        return u0, ut0, utt0

    # Вычисление напряжений, деформаций, скоростей и ускорений
    def __calc_dynamic_results__(self, u0, ut0, utt0, t):
        # Вычисление деформаций и напряжений
        super()._calc_results_(t)
        # Вычисление скоростей и ускорений (конечными разностями)
        th = self._params_.th
        for i in range(0, len(self._mesh_.x)):
            for j in range(0, self._mesh_.freedom):
                u_0 = u0[i * self._mesh_.freedom + j]
                u_1 = self.__global_load__[i * self._mesh_.freedom + j]
                u_t_0 = ut0[i * self._mesh_.freedom + j]
                u_t_1 = (u_1 - u_0)/th
                u_t_t_1 = (u_t_1 - u_t_0)/th
                u0[i * self._mesh_.freedom + j] = u_1
                ut0[i * self._mesh_.freedom + j] = u_t_1
                utt0[i * self._mesh_.freedom + j] = u_t_t_1
                self._result_[len(self._result_) - 2 * self._mesh_.freedom + j].results[i] = u_t_1
                self._result_[len(self._result_) - self._mesh_.freedom + j].results[i] = u_t_t_1
        return u0, ut0, utt0

    # Добавление ЛМЖ, ЛММ и ЛМД к ГМЖ
    def _assembly_(self, fe, index):
        # Добавление матрицы
        for i in range(0, len(fe.K)):
            k = self._mesh_.fe[index][i // self._mesh_.freedom] * self._mesh_.freedom + i % self._mesh_.freedom
            for j in range(i, len(fe.K)):
                r = self._mesh_.fe[index][j // self._mesh_.freedom] * self._mesh_.freedom + j % self._mesh_.freedom
                self.__global_matrix_stiffness__[k, r] += fe.K[i][j]
                self.__global_matrix_mass__[k, r] += fe.M[i][j]
                self.__global_matrix_damping__[k, r] += fe.C[i][j]
                if k != r:
                    self.__global_matrix_stiffness__[r, k] += fe.K[i][j]
                    self.__global_matrix_mass__[r, k] += fe.M[i][j]
                    self.__global_matrix_damping__[r, k] += fe.C[i][j]

    # Формирование левой части (матрицы) уравнения квазистатического равновесия
    def __create_dynamic_matrix__(self):
        theta = 1.37
        self._progress_.set_process('Creating static part of global matrix...', 1, 2)
        self.__global_matrix_stiffness__ += \
            self.__global_matrix_mass__.dot(6.0 / (theta ** 2 ** self._params_.th ** 2)) + \
            self.__global_matrix_damping__.dot(6.0 / (3.0 / (theta * self._params_.th)))
        self._progress_.set_progress(2)

    # Формирование правой части (вектора) уравнения квазистатического равновесия
    def __create_dynamic_vector__(self, u0, ut0, utt0, t):
        theta = 1.37
        k1 = 6.0/(theta ** 2 ** self._params_.th ** 2)
        k2 = 3.0/(theta * self._params_.th)
        k3 = 0.5*(theta * self._params_.th)
        self.__global_load__ = zeros(len(self._mesh_.x) * self._mesh_.freedom)
        # Вычисление компонент нагрузки для текущего момента времени
        self._prepare_concentrated_load_(t)
        self._prepare_surface_load_(t)
        self._prepare_volume_load_(t)
        self.__global_load__ += \
            self.__global_matrix_mass__.dot(u0.dot(k1) + ut0.dot(2.0*k2) + utt0.dot(2.0)) + \
            self.__global_matrix_damping__.dot(u0.dot(k2) + ut0.dot(2.0) + utt0.dot(k3))

    # Определение кол-ва результатов в зависимости от размерности задачи
    def _num_result_(self):
        res = 0
        if self._mesh_.freedom == 1:
            # u, Exx, Sxx, ut, utt
            res = 5
        elif self._mesh_.freedom == 2:
            # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
            res = 12
        elif self._mesh_.freedom == 3:
            # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
            res = 21
        return res

    # Индекс функции в зависимости от размерности задачи
    def _index_result_(self, i):
        ret = 0
        # u, Exx, Sxx, ut, utt
        index4 = [4, 7, 13, 19, 22]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
        index5 = [4, 5, 7, 8, 10, 13, 14, 16, 19, 20, 22, 23]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
        index6 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        if self._mesh_.freedom == 1:
            ret = index4[i]
        elif self._mesh_.freedom == 2:
            ret = index5[i]
        elif self._mesh_.freedom == 3:
            ret = index6[i]
        return ret
