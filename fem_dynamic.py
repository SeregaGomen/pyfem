#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

from fem_fem import TFEM
from scipy.sparse import lil_matrix


class TFEMDynamic(TFEM):
    def __init__(self):
        super().__init__()

    # Расчет динамической задачи методом конечных элементов
    def __calc_problem__(self):
        # Создание глобальных матриц жесткости, масс и демпфирования
        size = len(self.__mesh__.x)*self.__mesh__.freedom
        self.__global_matrix__ = lil_matrix((size, size))
        self.__global_vector__ = [0]*size

        self.__progress__.set_process('Assembling global stiffness matrix...', 1, len(self.__mesh__.fe))
        for i in range(0, len(self.__mesh__.fe)):
            self.__progress__.set_progress(i + 1)
            x = [0]*len(self.__mesh__.fe[i])
            y = [0]*len(self.__mesh__.fe[i])
            z = [0]*len(self.__mesh__.fe[i])
            vx = [0]*len(self.__mesh__.fe[i])
            vy = [0]*len(self.__mesh__.fe[i])
            vz = [0]*len(self.__mesh__.fe[i])
            # Настройка КЭ
            for j in range(len(self.__mesh__.fe[i])):
                x[j], y[j], z[j] = self.__mesh__.get_coord(self.__mesh__.fe[i][j])
            fe.set_coord(x, y, z)
            fe.generate()
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly__(fe, i)

