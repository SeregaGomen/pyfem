#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#              Параметры решения задачи с помощью МКЭ
###################################################################

from core.fem_error import TFEMException


# Тип задачи: статика или динамика
problem_type = [
    'static',
    'dynamic'
]

# Метод решения СЛАУ (точный, приближенный
solve_method = [
    'direct',
    'iterative'
]

# Стандартные имена функций (перемещения, деформации и напряжения) и их агрументов
std_name = [
    'x',    # 0  - идентификатор первого аргумента (x)
    'y',    # 1  - ... второго аргумента
    'z',    # 2  - ... третьего аргумента
    't',    # 3  - ... времени
    'U',    # 4  - компонента вектора перемещений по первому направлению (x)
    'V',    # 5  - ... по y
    'W',    # 6  - ... по z
    'Exx',  # 7  - компонента тензора нормальных деформаций по xx
    'Eyy',  # 8  - ... по yy
    'Ezz',  # 9  - ... по zz
    'Exy',  # 10 - компонента тензора тангенциальных деформаций по xу
    'Exz',  # 11 - ... по xz
    'Eyz',  # 12 - ... по yz
    'Sxx',  # 13 - компонента тензора нормальных напряжений по xx
    'Syy',  # 14 - ... по yy
    'Szz',  # 15 - ... по zz
    'Sxy',  # 16 - компонента тензора тангенциальных напряжений по xу
    'Sxz',  # 17 - ... по xz
    'Syz',  # 18 - ... по yz
    'Ut',   # 19 - скорость по x
    'Vt',   # 20 - ... по y
    'Wt',   # 21 - ... по z
    'Utt',  # 22 - ускорение по x
    'Vtt',  # 23 - ... по y
    'Wtt'   # 24 - ... по z
]

# Типы условий (краевых, упругих, ...)
type_condition = [
    'initial',
    'boundary',
    'volume',
    'surface',
    'concentrated',
    'thickness',
    'young_modulus',
    'poisson_ratio',
    'temperature',
    'alpha'
]


# Описание условия
class TCondition:
    def __init__(self):
        self.direct = 0         # Направление (номер функции, для которой задается условие): 0 - по x и т.д.
        self.expression = ''    # Функциональное выражение, определяющее значение условия (например, 10^5)
        self.predicate = ''     # Предикат отбора узлов
        self.type = ''          # Тип краевого условия


# Базовые параметры расчета задачи теории упругости с помощью МКЭ
class TFEMParams:
    def __init__(self):
        self.problem_type = ''  # Тип задачи
        self.solve_method = ''  # Метод решения СЛАУ
        self.width = 12         # Формат вывода результатов
        self.precision = 5
        self.eps = 1.0E-6       # Точность вычислений
        self.density = 0        # Плотность материала
        self.damping = 0        # Параметр демпфирования
        self.t0 = 0             # Начальный момент времени расчета
        self.t1 = 0             # Конечный момент времени расчета
        self.th = 0             # Шаг по времени
        self.alpha = 0          # Коэффициент теплового расширения
        self.dT = 0             # Разность температур
        self.names = std_name   # Список имен функций и их аргументов
        self.bc_list = []       # Список краевых условий
        self.var_list = {}      # Список вспомогательных переменных и их значений

    def __add_condition(self, t, e, p, d=0):
        c = TCondition()
        c.type = t
        c.direct = d
        c.expression = e
        c.predicate = p
        self.bc_list.append(c)

    def add_boundary_condition(self, e, p, d):
        self.__add_condition('boundary', e, p, d)

    def add_initial_condition(self, e, p, d):
        self.__add_condition('initial', e, p, d)

    def add_volume_load(self, e, p, d):
        self.__add_condition('volume', e, p, d)

    def add_surface_load(self, e, p, d):
        self.__add_condition('surface', e, p, d)

    def add_concentrated_load(self, e, p, d):
        self.__add_condition('concentrated', e, p, d)

    def add_thickness(self, e, p):
        self.__add_condition('thickness', e, p)

    def add_young_modulus(self, e, p):
        self.__add_condition('young_modulus', e, p)

    def add_poisson_ratio(self, e, p):
        self.__add_condition('poisson_ratio', e, p)

    def add_temperature(self, e, p):
        self.__add_condition('temperature', e, p)

    def add_alpha(self, e, p):
        self.__add_condition('alpha', e, p)

    def add_variable(self, var, val):
        self.var_list.setdefault(var, val)

    def check_params(self):
        if self.solve_method == '':
            raise TFEMException('solve_method_err')
        if self.problem_type == '':
            raise TFEMException('problem_type_err')
        is_young_modulus = False
        is_poisson_ratio = False
        for i in range(len(self.bc_list)):
            if self.bc_list[i].type == 'young_modulus':
                is_young_modulus = True
            if self.bc_list[i].type == 'poisson_ratio':
                is_poisson_ratio = True
        if not is_young_modulus or not is_poisson_ratio:
            raise TFEMException('elasticity_err')
        if self.problem_type == 'dynamic' and (self.t0 == self.t1 or self.th <= 0):
            raise TFEMException('time_err')
