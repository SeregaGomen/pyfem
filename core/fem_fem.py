#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################################
# Реализация вычислений задач заданного типа методо конечных элементов
#######################################################################

from abc import abstractmethod
from core.fem_mesh import TMesh
from core.fem_params import TFEMParams
from core.fem_progress import TProgress
from core.fem_fe import TFE, TFE1D2, TFE2D3, TFE2D3P, TFE2D3S, TFE2D4, TFE2D4P, TFE2D4S, TFE3D4, TFE3D8
from core.fem_parser import TParser
from core.fem_error import TFEMException


# Абстрактный базовый класс, реализующий МКЭ
class TFEM:
    def __init__(self):
        self._mesh_ = TMesh()                                 # Дискретная модель объекта
        self._params_ = TFEMParams()                          # Параметры расчета
        self._progress_ = TProgress()                         # Индикатор прогресса расчета
        self._result_ = []                                    # Список результатов расчета

    @abstractmethod
    def _calc_problem_(self):
        raise NotImplementedError('Method TFEM._calc_problem_ is pure virtual')

    # Добавление локальной матрицы жесткости (масс, демпфирования) к глобальной
    @abstractmethod
    def _assembly_(self, fe, index):
        raise NotImplementedError('Method TFEM._assembly_ is pure virtual')

    # Вычисление вспомогательных результатов (деформаций, напряжений, ...) 
    @abstractmethod
    def _calc_results_(self):
        raise NotImplementedError('Method TFEM.calc_results is pure virtual')

    # Прямое решение СЛАУ
    @abstractmethod
    def _solve_direct_(self):
        raise NotImplementedError('Method TFEM._solve_direct_ is pure virtual')

    # Приближенное решение СЛАУ
    @abstractmethod
    def _solve_iterative_(self):
        raise NotImplementedError('Method TFEM._solve_iterative_ is pure virtual')

    # Решение СЛАУ
    def _solve_(self):
        ret = False
        if self._params_.solve_method == 'direct':
            ret = self._solve_direct_()
        elif self._params_.solve_method == 'iterative':
            ret = self._solve_iterative_()
        return ret

    # Создание нужного типа КЭ
    def _create_fe_(self):
        fe = TFE()
        if self._mesh_.fe_type == 'fe_1d_2':
            fe = TFE1D2()
        elif self._mesh_.fe_type == 'fe_2d_3':
            fe = TFE2D3()
        elif self._mesh_.fe_type == 'fe_2d_4':
            fe = TFE2D4()
        elif self._mesh_.fe_type == 'fe_3d_4':
            fe = TFE3D4()
        elif self._mesh_.fe_type == 'fe_3d_8':
            fe = TFE3D8()
        elif self._mesh_.fe_type == 'fe_2d_3_p':
            fe = TFE2D3P()
        elif self._mesh_.fe_type == 'fe_2d_4_p':
            fe = TFE2D4P()
        elif self._mesh_.fe_type == 'fe_2d_3_s':
            fe = TFE2D3S()
        elif self._mesh_.fe_type == 'fe_2d_4_s':
            fe = TFE2D4S()
        return fe

    # Создание и настройка парсера
    def _create_parser_(self):
        parser = TParser()
        for i in range(0, len(self._params_.names)):
            parser.add_variable(self._params_.names[i])
        for key, value in self._params_.var_list.items():
            parser.add_variable(key, value)
        return parser

    # Вычисление значения нагрузки
    def _setup_parser_(self, parser, x, t, p_code, f_code):
        for i in range(0, len(x)):
            parser.set_variable(self._params_.names[i], x[i])
        parser.set_variable(self._params_.names[3], t)
        if len(p_code):
            parser.set_code(p_code)
            if parser.run() == 0:
                return False, 0
        parser.set_code(f_code)
        return True, parser.run()

    # Запуск процедуры расчета
    def calc(self):
        try:
            # Проверка наличия и соответствия необходимых параметров расчета
            self._params_.check_params()
            ret = self._calc_problem_()
        except TFEMException as err:
            ret = False
            err.print_error()
        return ret

    # Задание сетки
    def set_mesh(self, mesh):
        self._mesh_ = mesh

    # Задание параметров расчета
    def set_params(self, params):
        self._params_ = params

    # Возврат результатов расчета
    def get_result(self):
        return self._result_
