#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
# Реализация вычислений задач заданного типа методо конечных элементов
#######################################################################

from abc import abstractmethod
from core.fem_mesh import TMesh
from core.fem_params import TFEMParams
from core.fem_progress import TProgress
from core.fem_fe import TFE, TFE1D2, TFE2D3, TFE2D3P, TFE2D3S, TFE2D4, TFE2D6, TFE2D6P, TFE2D6S, TFE2D4P, TFE2D4S, \
    TFE3D4, TFE3D8, TFE3D10
from core.fem_parser import TParser
from core.fem_error import TException


# Абстрактный базовый класс, реализующий МКЭ
class TFEM:
    def __init__(self):
        self.mesh = TMesh()             # Дискретная модель объекта
        self.params = TFEMParams()      # Параметры расчета
        self.results = []               # Список результатов расчета
        self._progress = TProgress()    # Индикатор прогресса расчета

    @abstractmethod
    def _calc_problem(self):
        raise NotImplementedError('Method TFEM._calc_problem is pure virtual')

    # Добавление локальной матрицы жесткости (масс, демпфирования) к глобальной
    @abstractmethod
    def __assembly(self, fe, index):
        raise NotImplementedError('Method TFEM.__assembly is pure virtual')

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...) 
    @abstractmethod
    def _calc_results(self):
        raise NotImplementedError('Method TFEM._calc_results is pure virtual')

    # Прямое решение СЛАУ
    @abstractmethod
    def _solve_direct(self):
        raise NotImplementedError('Method TFEM._solve_direct is pure virtual')

    # Приближенное решение СЛАУ
    @abstractmethod
    def _solve_iterative(self):
        raise NotImplementedError('Method TFEM._solve_iterative is pure virtual')

    # Решение СЛАУ
    def _solve(self):
        ret = False
        if self.params.solve_method == 'direct':
            ret = self._solve_direct()
        elif self.params.solve_method == 'iterative':
            ret = self._solve_iterative()
        return ret

    # Создание нужного типа КЭ
    def create_fe(self):
        fe = TFE()
        if self.mesh.fe_type == 'fe_1d_2':
            fe = TFE1D2()
        elif self.mesh.fe_type == 'fe_2d_3':
            fe = TFE2D3()
        elif self.mesh.fe_type == 'fe_2d_4':
            fe = TFE2D4()
        elif self.mesh.fe_type == 'fe_2d_6':
            fe = TFE2D6()
        elif self.mesh.fe_type == 'fe_3d_4':
            fe = TFE3D4()
        elif self.mesh.fe_type == 'fe_3d_8':
            fe = TFE3D8()
        elif self.mesh.fe_type == 'fe_3d_10':
            fe = TFE3D10()
        elif self.mesh.fe_type == 'fe_2d_3_p':
            fe = TFE2D3P()
        elif self.mesh.fe_type == 'fe_2d_6_p':
            fe = TFE2D6P()
        elif self.mesh.fe_type == 'fe_2d_4_p':
            fe = TFE2D4P()
        elif self.mesh.fe_type == 'fe_3d_3_s':
            fe = TFE2D3S()
        elif self.mesh.fe_type == 'fe_3d_6_s':
            fe = TFE2D6S()
        elif self.mesh.fe_type == 'fe_3d_4_s':
            fe = TFE2D4S()
        return fe

    # Создание и настройка парсера
    def create_parser(self, x, t=0):
        parser = TParser()
        for i in range(0, len(x)):
            parser.add_variable(self.params.names[i], x[i])
        parser.add_variable(self.params.names[3], t)
        for key, value in self.params.var_list.items():
            parser.add_variable(key, value)
        return parser

    # Запуск процедуры расчета
    def calc(self):
        try:
            # Проверка наличия и соответствия необходимых параметров расчета
            self.params.check_params()
            ret = self._calc_problem()
        except TException as err:
            ret = False
            err.print_error()
        return ret

    # Задание сетки
    def set_mesh(self, mesh):
        self.mesh = mesh

    # Задание параметров расчета
    def set_params(self, params):
        self.params = params

    # Возврат результатов расчета
    def get_results(self, *params):
        if not len(params):
            return self.results
        return self.results[params[0]]
