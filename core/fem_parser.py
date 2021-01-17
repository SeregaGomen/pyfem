#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#        Обработка арифметических и логических выражений
###################################################################

from core.fem_error import TException


# Класс, реализующий разбор и выполнение арифметических и логических выражений
class TParser:
    def __init__(self):
        self.error = self.code = ''
        self.variables = {}

    def add_variable(self, var, val=0.0):
        if var in self.variables:
            self.error = 'redefinition_err'
        self.variables.setdefault(var, val)

    def set_code(self, c):
        self.code = c
#        try:
#            compile(self.code, 'fem_parser.py', 'single')
#        except Exception:
#            self.error = 'syntax_err'
#            raise TFEMException('syntax_err')

    def run(self):
        exec('from math import sin, cos, tan, exp, asin, acos, atan, atan2, sinh, cosh')
        try:
            for var in self.variables:
                exec(var + ' = ' + str(self.variables[var]))
            res = eval(self.code)
        except Exception:
            self.error = 'calc_err'
            raise TException('calc_err')
        return res
