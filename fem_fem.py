#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
# Реализация вычислений задач заданного типа методо конечных элементов
#######################################################################

from abc import abstractmethod
from fem_mesh import TMesh
from fem_params import TFEMParams
from fem_progress import TProgress
from fem_fe import TFE, TFE1D2, TFE2D3, TFE2D4, TFE3D4, TFE3D8
from fem_parser import TParser
from fem_error import TFEMException


# Абстрактный базовый класс, реализующий МКЭ
class TFEM:
    def __init__(self):
        self.__mesh__ = TMesh()                                 # Дискретная модель объекта
        self.__params__ = TFEMParams()                          # Параметры расчета
        self.__progress__ = TProgress()                         # Индикатор прогресса расчета
        self.__result__ = []                                    # Список результатов расчета

    @abstractmethod
    def __calc_problem__(self):
        raise NotImplementedError('Method TFEM.__calc_problem__ is pure virtual')

    # Добавление локальной матрицы жесткости (масс, демпфирования) к глобальной
    @abstractmethod
    def __assembly__(self, fe, index):
        raise NotImplementedError('Method TFEM.__assembly__ is pure virtual')

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...) 
    @abstractmethod
    def __calc_results__(self):
        raise NotImplementedError('Method TFEM.calc_results is pure virtual')

    # Решение СЛАУ
    def __solve__(self):
        ret = False
        if self.__params__.solve_method == 'direct':
            ret = self.__solve_direct__()
        elif self.__params__.solve_method == 'iterative':
            ret = self.__solve_iterative__()
        return ret

    # Прямое решение СЛАУ
    @abstractmethod
    def __solve_direct__(self):
        raise NotImplementedError('Method TFEM.__solve_direct__ is pure virtual')

    # Приближенное решение СЛАУ
    def __solve_iterative__(self):
        raise NotImplementedError('Method TFEM.__solve_iterative__ is pure virtual')

    # Создание нужного типа КЭ
    def __create_fe__(self):
        fe = TFE()
        if self.__mesh__.fe_type == 'fe_1d_2':
            fe = TFE1D2()
        elif self.__mesh__.fe_type == 'fe_2d_3':
            fe = TFE2D3()
        elif self.__mesh__.fe_type == 'fe_2d_4':
            fe = TFE2D4()
        elif self.__mesh__.fe_type == 'fe_3d_4':
            fe = TFE3D4()
        elif self.__mesh__.fe_type == 'fe_3d_8':
            fe = TFE3D8()
        return fe

    # Определение кол-ва результатов в зависимости от типа и размерности задачи
    def __num_result__(self):
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
    def __index_result__(self, i):
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

    # Настройка парсера
    def __create_parser__(self):
        parser = TParser()
        for i in range(0, len(self.__params__.names)):
            parser.add_variable(self.__params__.names[i])
        for key, value in self.__params__.var_list.items():
            parser.add_variable(key, value)
        return parser

    # Запуск процедуры расчета
    def calc(self):
        try:
            ret = self.__calc_problem__()
        except TFEMException as err:
            ret = False
            err.print_error()
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
