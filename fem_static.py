#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

from scipy.sparse import lil_matrix
from fem_fem import TFEM
from fem_defs import DIR_X, DIR_Y, DIR_Z
from fem_error import TFEMException


class TFEMStatic(TFEM):
    def __init__(self):
        super().__init__()

    # Запуск процедуры расчета
    def calc(self):
        try:
            ret = self.__calc_problem__()
        except TFEMException as err:
            ret = False
            err.print_error()
        return ret

    # Расчет статической задачи методом конечных элементов
    def __calc_problem__(self):
        # Создание ГМЖ
        size = len(self.__mesh__.x)*self.__mesh__.freedom
        self.__global_matrix__ = lil_matrix((size, size))
        self.__global_vector__ = [0]*size

        fe = self.create_fe()
        fe.set_elasticity(self.__params__.e, self.__params__.m)
        # Предварительное вычисление компонент нагрузки
        surface_attribute = self.prepare_load()
        # Формирование глобальной матрицы жесткости
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
                vx[j] = self.__volume_load__[self.__mesh__.fe[i][j]*self.__mesh__.freedom + 0]
                vy[j] = self.__volume_load__[self.__mesh__.fe[i][j]*self.__mesh__.freedom + 1] \
                    if (len(self.__mesh__.y)) else 0
                vz[j] = self.__volume_load__[self.__mesh__.fe[i][j]*self.__mesh__.freedom + 2] \
                    if (len(self.__mesh__.z)) else 0
            fe.set_coord(x, y, z)
            fe.set_volume_load(vx, vy, vz)
            fe.generate()
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly__(fe, i)
        # Учет сосредоточенной и поверхностной нагрузок
        self.use_load_condition(surface_attribute)
        # Учет краевых условий
        self.use_boundary_condition()
        # Решение СЛАУ
        if not self.solve():
            print('The system of equations is not solved!')
            return False
        self.calc_results()
        print('**************** Success! ****************')
        return True

    # Добавление локальной матрицы жесткости (ЛМЖ) к ГМЖ
    def __assembly__(self, fe, index):
        # Добавление матрицы
        for i in range(0, len(fe.K)):
            k = self.__mesh__.fe[index][i//self.__mesh__.freedom]*self.__mesh__.freedom + i % self.__mesh__.freedom
            for j in range(i, len(fe.K)):
                l = self.__mesh__.fe[index][j//self.__mesh__.freedom]*self.__mesh__.freedom + j % self.__mesh__.freedom
                self.__global_matrix__[k, l] += fe.K[i][j]
                if k != l:
                    self.__global_matrix__[l, k] += fe.K[i][j]
            self.__global_vector__[k] += fe.K[i][len(fe.K)]

    # Предварительное вычисление нагрузок
    def prepare_load(self):
        surface_attribute = []
        parser = self.create_parser()
        self.__volume_load__ = [0]*len(self.__mesh__.x)*self.__mesh__.freedom
        self.__surface_load__ = [0]*len(self.__mesh__.x)*self.__mesh__.freedom
        self.__concentrated_load__ = [0]*len(self.__mesh__.x)*self.__mesh__.freedom

        counter = 0
        for i in range(0, len(self.__params__.bc_list)):
            if not (self.__params__.bc_list[i].type == 'volume' or self.__params__.bc_list[i].type == 'surface' or
               self.__params__.bc_list[i].type == 'concentrated'):
                continue
            counter += 1

        self.__progress__.set_process('Computation of load...', 1, counter*len(self.__mesh__.x))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if not (self.__params__.bc_list[i].type == 'volume' or self.__params__.bc_list[i].type == 'surface' or
               self.__params__.bc_list[i].type == 'concentrated'):
                continue
            if self.__params__.bc_list[i].type == 'surface':
                # Проверяем, все ли узлы ГЭ удовлетворяют предикату отбора
                surface_attribute = self.check_boundary_elements(self.__params__.bc_list[i].predicate)

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
                        return surface_attribute
                    if parser.run() == 0:
                        continue
                parser.set_code(self.__params__.bc_list[i].expression)
                if parser.error != '':
                    return surface_attribute
                val = parser.run()
                if self.__params__.bc_list[i].direct & DIR_X:
                    index = j*self.__mesh__.freedom + 0
                    self.add_load(i, index, val)
                if self.__params__.bc_list[i].direct & DIR_Y:
                    index = j*self.__mesh__.freedom + 1
                    self.add_load(i, index, val)
                if self.__params__.bc_list[i].direct & DIR_Z:
                    index = j*self.__mesh__.freedom + 2
                    self.add_load(i, index, val)
        return surface_attribute

    # Учет сосредоточенной и поверхностной нагрузок
    def use_load_condition(self, surface_attribute):
        self.__progress__.set_process('Building the load vector-column...', 1,
                                      len(self.__concentrated_load__) + len(self.__mesh__.surface))
        counter = 1
        # Учет сосредоточенной нагрузки
        for i in range(0, len(self.__concentrated_load__)):
            self.__progress__.set_progress(counter)
            counter += 1
            self.__global_vector__[i] += self.__concentrated_load__[i]
        # Учет поверхностной нагрузки
        if self.__mesh__.freedom == 1:
            return
        x = [0]*len(self.__mesh__.surface[0])
        y = [0]*len(self.__mesh__.surface[0])
        z = [0]*len(self.__mesh__.surface[0])
        for i in range(0, len(self.__mesh__.surface)):
            self.__progress__.set_progress(counter)
            counter += 1
            if not len(surface_attribute) or not surface_attribute[i]:
                continue
            for j in range(0, len(self.__mesh__.surface[0])):
                x[j], y[j], z[j] = self.__mesh__.get_coord(self.__mesh__.surface[i][j])
            rel_se = self.square(x, y, z)/float(len(self.__mesh__.surface[0]))
            for j in range(0, len(self.__mesh__.surface[0])):
                for k in range(0, self.__mesh__.freedom):
                    l = self.__mesh__.surface[i][j]*self.__mesh__.freedom + k
                    self.__global_vector__[l] += self.__surface_load__[l]*rel_se