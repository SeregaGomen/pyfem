#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

import math
from os import remove
from scipy.sparse import lil_matrix, coo_matrix
from numpy import zeros, savez, load
from fem_defs import INIT_U, INIT_V, INIT_W, INIT_U_T, INIT_V_T, INIT_W_T, INIT_U_T_T, INIT_V_T_T, INIT_W_T_T
from fem_static import TFEMStatic


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
        self.__initial_condition__ = [[]]                       # Начальные условия

    # Расчет динамической задачи методом конечных элементов
    def __calc_problem__(self):
        size = len(self.__mesh__.x)*self.__mesh__.freedom
        self.__global_matrix_stiffness__ = lil_matrix((size, size))
        self.__global_matrix_mass__ = lil_matrix((size, size))
        self.__global_matrix_damping__ = lil_matrix((size, size))
        self.__global_load__ = [0]*size
        fe = self.__create_fe__()
        fe.set_elasticity(self.__params__.e, self.__params__.m)
        fe.set_damping(self.__params__.damping)
        fe.set_density(self.__params__.density)

        # Создание глобальных матриц жесткости, масс и демпфирования
        self.__progress__.set_process('Assembling global stiffness, mass and damping matrix...', 1,
                                      len(self.__mesh__.fe))
        for i in range(0, len(self.__mesh__.fe)):
            self.__progress__.set_progress(i + 1)
            x = [0]*len(self.__mesh__.fe[i])
            y = [0]*len(self.__mesh__.fe[i])
            z = [0]*len(self.__mesh__.fe[i])
            # Настройка КЭ
            for j in range(len(self.__mesh__.fe[i])):
                x[j], y[j], z[j] = self.__mesh__.get_coord(self.__mesh__.fe[i][j])
            fe.set_coord(x, y, z)
            fe.generate(False)
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly__(fe, i)
        # Формирование статической (левой) части СЛАУ
        self.__create_static_matrix__()
        # Сохранение матрицы для последующего использования
        save_matrix('tmp_matrix', self.__global_matrix_stiffness__)
        # Учет начальных условий
        self.__use_initial_condition__()
        # Итерационный процесс по времени
        t = self.__params__.t0
        while t <= self.__params__.t1:
            print('t = %5.2f' % t)
            # Вычисление компонент нагрузки для текущего момента времени
            self.__prepare_concentrated_load__(t)
            self.__prepare_surface_load__(t)
            self.__prepare_volume_load__(t)
            # Учет краевых условий
            self.__use_boundary_condition__()
            # Решение СЛАУ
            if not self.__solve__():
                print('The system of equations is not solved!')
                return False
            self.__calc_results__()
            t += self.__params__.th
            if math.fabs(t - self.__params__.t1) < self.__params__.eps:
                t = self.__params__.t1
            self.__global_matrix_stiffness__ = load_matrix('tmp_matrix.npz')
        # Удаляем временный файл с матрицей
        remove('tmp_matrix.npz')
        print('**************** Success! ****************')
        return True

    # Учет начальных условий
    def __use_initial_condition__(self):
        parser = self.__create_parser__()
        counter = 0
        self.__initial_condition__ = zeros((3, len(self.__mesh__.x)*self.__mesh__.freedom))
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'initial':
                counter += 1
        self.__progress__.set_process('Using initial conditions...', 1, counter*len(self.__mesh__.x))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'initial':
                parser.set_code(self.__params__.bc_list[i].expression)
                value = parser.run()
                direct = self.__params__.bc_list[i].direct
                for j in range(0, len(self.__mesh__.x)):
                    self.__progress__.set_progress(counter)
                    counter += 1
                    if direct & INIT_U:
                        self.__initial_condition__[0][j*self.__mesh__.freedom + 0] = value
                    if direct & INIT_V:
                        self.__initial_condition__[0][j*self.__mesh__.freedom + 1] = value
                    if direct & INIT_W:
                        self.__initial_condition__[0][j*self.__mesh__.freedom + 2] = value
                    if direct & INIT_U_T:
                        self.__initial_condition__[1][j*self.__mesh__.freedom + 0] = value
                    if direct & INIT_V_T:
                        self.__initial_condition__[1][j*self.__mesh__.freedom + 1] = value
                    if direct & INIT_W_T:
                        self.__initial_condition__[1][j*self.__mesh__.freedom + 2] = value
                    if direct & INIT_U_T_T:
                        self.__initial_condition__[2][j*self.__mesh__.freedom + 0] = value
                    if direct & INIT_V_T_T:
                        self.__initial_condition__[2][j*self.__mesh__.freedom + 1] = value
                    if direct & INIT_W_T_T:
                        self.__initial_condition__[2][j*self.__mesh__.freedom + 2] = value

    # Добавление ЛМЖ, ЛММ и ЛМД к ГМЖ
    def __assembly__(self, fe, index):
        # Добавление матрицы
        for i in range(0, len(fe.K)):
            k = self.__mesh__.fe[index][i//self.__mesh__.freedom]*self.__mesh__.freedom + i % self.__mesh__.freedom
            for j in range(i, len(fe.K)):
                l = self.__mesh__.fe[index][j//self.__mesh__.freedom]*self.__mesh__.freedom + j % self.__mesh__.freedom
                self.__global_matrix_stiffness__[k, l] += fe.K[i][j]
                self.__global_matrix_mass__[k, l] += fe.M[i][j]
                self.__global_matrix_damping__[k, l] += fe.D[i][j]
                if k != l:
                    self.__global_matrix_stiffness__[l, k] += fe.K[i][j]
                    self.__global_matrix_mass__[l, k] += fe.M[i][j]
                    self.__global_matrix_damping__[l, k] += fe.D[i][j]
            self.__global_load__[k] += fe.K[i][len(fe.K)]

    # Добавление к матрице жесткости матриц масс и демпфирования
    def __create_static_matrix__(self):
        theta = 1.37
        size1 = len(self.__global_matrix_mass__.nonzero()[0])
        size2 = len(self.__global_matrix_damping__.nonzero()[0])
        counter = 1
        k1 = 3.0/(theta*self.__params__.th)
        k2 = 6.0/(theta**2**self.__params__.th**2)
        self.__progress__.set_process('Creating static part of global matrix...', 1, size1 + size2)
        for m in range(0, size1):
            self.__progress__.set_progress(counter)
            i = self.__global_matrix_mass__.nonzero()[0][m]
            j = self.__global_matrix_mass__.nonzero()[1][m]
            self.__global_matrix_stiffness__[i, j] += k1*self.__global_matrix_mass__[i, j]
            counter += 1
        for m in range(0, size2):
            self.__progress__.set_progress(counter)
            i = self.__global_matrix_damping__.nonzero()[0][m]
            j = self.__global_matrix_damping__.nonzero()[1][m]
            self.__global_matrix_stiffness__[i, j] += k2*self.__global_matrix_damping__[i, j]
            counter += 1

    # Определение кол-ва результатов в зависимости от размерности задачи
    def __num_result__(self):
        res = 0
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

    # Индекс функции в зависимости от размерности задачи
    def __index_result__(self, i):
        ret = 0
        # u, Exx, Sxx, ut, utt
        index4 = [4, 7, 13, 19, 22]
        # u, v, Exx, Eyy, Exy, Sxx, Syy, Sxy, ut, vt, utt, vtt
        index5 = [4, 5, 7, 8, 10, 13, 14, 16, 19, 20, 22, 23]
        # u, v, w, Exx, Eyy, Ezz, Exy, Exz, Eyz, Sxx, Syy, Szz, Sxy, Sxz, Syz, ut, utt, vt, vtt, wt, wtt
        index6 = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        if self.__mesh__.freedom == 1:
            ret = index4[i]
        elif self.__mesh__.freedom == 2:
            ret = index5[i]
        elif self.__mesh__.freedom == 3:
            ret = index6[i]
        return ret
