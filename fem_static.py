#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
#           Класс, реализующий расчет статической задачи
#######################################################################

from scipy.sparse import lil_matrix
from fem_fem import TFEM
from fem_defs import DIR_X, DIR_Y, DIR_Z
from fem_result import TResult


class TFEMStatic(TFEM):
    def __init__(self):
        super().__init__()

    # Расчет статической задачи методом конечных элементов
    def __calc_problem__(self):
        # Создание ГМЖ
        size = len(self.__mesh__.x)*self.__mesh__.freedom
        self.__global_matrix__ = lil_matrix((size, size))
        self.__global_vector__ = [0]*size

        fe = self.create_fe()
        fe.set_elasticity(self.__params__.e, self.__params__.m)
        # Вычисление компонент нагрузки
        self.prepare_concentrated_load()
        self.prepare_surface_load()
        self.prepare_volume_load()
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
            fe.set_coord(x, y, z)
            fe.generate()
            # Ансамблирование ЛМЖ к ГМЖ
            self.__assembly__(fe, i)
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

    # Вычисление сосредоточенных нагрузок
    def prepare_concentrated_load(self):
        parser = self.create_parser()
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
                if self.__params__.bc_list[i].direct & DIR_X:
                    self.__global_vector__[j*self.__mesh__.freedom + 0] += val
                if self.__params__.bc_list[i].direct & DIR_Y:
                    self.__global_vector__[j*self.__mesh__.freedom + 1] += val
                if self.__params__.bc_list[i].direct & DIR_Z:
                    self.__global_vector__[j*self.__mesh__.freedom + 2] += val

    # Вычисление поверхностных нагрузок
    def prepare_surface_load(self):
        x = [0]*len(self.__mesh__.surface[0])
        y = [0]*len(self.__mesh__.surface[0])
        z = [0]*len(self.__mesh__.surface[0])
        val = [0]*len(self.__mesh__.surface[0])
        parser = self.create_parser()
        counter = 0
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type == 'surface':
                counter += 1
        if not counter:
            return
        self.__progress__.set_process('Computation of surface load...', 1, counter*len(self.__mesh__.surface))
        counter = 1
        for i in range(0, len(self.__params__.bc_list)):
            if self.__params__.bc_list[i].type != 'surface':
                continue
            for j in range(0, len(self.__mesh__.surface)):
                self.__progress__.set_progress(counter)
                counter += 1
                if not self.check_boundary_elements(j, self.__params__.bc_list[i].predicate):
                    continue
                rel_se = self.square(j)/float(len(self.__mesh__.surface[j]))
                for k in range(0, len(self.__mesh__.surface[j])):
                    x[k], y[k], z[k] = self.__mesh__.get_coord(self.__mesh__.surface[j][k])
                    parser.set_variable(self.__params__.names[0], x[k])
                    parser.set_variable(self.__params__.names[1], y[k])
                    parser.set_variable(self.__params__.names[2], z[k])
                    parser.set_code(self.__params__.bc_list[i].expression)
                    val[k] = parser.run()
                    if self.__params__.bc_list[i].direct & DIR_X:
                        self.__global_vector__[self.__mesh__.surface[j][k]*self.__mesh__.freedom + 0] += val[k]*rel_se
                    if self.__params__.bc_list[i].direct & DIR_Y:
                        self.__global_vector__[self.__mesh__.surface[j][k]*self.__mesh__.freedom + 1] += val[k]*rel_se
                    if self.__params__.bc_list[i].direct & DIR_Z:
                        self.__global_vector__[self.__mesh__.surface[j][k]*self.__mesh__.freedom + 2] += val[k]*rel_se

    # Вычисление объемных нагрузок
    def prepare_volume_load(self):
        x = [0]*len(self.__mesh__.fe[0])
        y = [0]*len(self.__mesh__.fe[0])
        z = [0]*len(self.__mesh__.fe[0])
        val = [0]*len(self.__mesh__.fe[0])
        parser = self.create_parser()
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
                rel_ve = self.volume(j)/float(len(self.__mesh__.fe[j]))
                for k in range(0, len(self.__mesh__.fe[j])):
                    x[k], y[k], z[k] = self.__mesh__.get_coord(self.__mesh__.fe[j][k])
                    parser.set_variable(self.__params__.names[0], x[k])
                    parser.set_variable(self.__params__.names[1], y[k])
                    parser.set_variable(self.__params__.names[2], z[k])
                    parser.set_code(self.__params__.bc_list[i].expression)
                    val[k] = parser.run()
                    if self.__params__.bc_list[i].direct & DIR_X:
                        self.__global_vector__[self.__mesh__.fe[j][k]*self.__mesh__.freedom + 0] += val[k]*rel_ve
                    if self.__params__.bc_list[i].direct & DIR_Y:
                        self.__global_vector__[self.__mesh__.fe[j][k]*self.__mesh__.freedom + 1] += val[k]*rel_ve
                    if self.__params__.bc_list[i].direct & DIR_Z:
                        self.__global_vector__[self.__mesh__.fe[j][k]*self.__mesh__.freedom + 2] += val[k]*rel_ve

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
        fe = self.create_fe()
        fe.set_elasticity(self.__params__.e, self.__params__.m)
        self.__progress__.set_process('Calculation results...', 1, len(self.__mesh__.fe))
        for i in range(0, len(self.__mesh__.fe)):
            self.__progress__.set_progress(i + 1)
            x = [0]*len(self.__mesh__.fe[i])
            y = [0]*len(self.__mesh__.fe[i])
            z = [0]*len(self.__mesh__.fe[i])
            for j in range(len(self.__mesh__.fe[i])):
                x[j], y[j], z[j] = self.__mesh__.get_coord(self.__mesh__.fe[i][j])
            fe.set_coord(x, y, z)
            for j in range(0, len(self.__mesh__.fe[i])):
                for k in range(0, self.__mesh__.freedom):
                    uvw[j*self.__mesh__.freedom + k] = \
                        self.__global_vector__[self.__mesh__.freedom*self.__mesh__.fe[i][j] + k]
            r = fe.calc(uvw)
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
