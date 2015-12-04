#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#                   Реализация обработки ошибок
###################################################################

import sys


class TFEMException(Exception):
    def __init__(self, e):
        self.error = e

    def print_error(self):
        sys.stderr.write('Error: ')
        if self.error == 'brackets_err':
            sys.stderr.write('unbalanced or unexpected bracket\n')
        elif self.error == 'syntax_err':
            sys.stderr.write('syntax error\n')
        elif self.error == 'undef_err':
            sys.stderr.write('undefined variable or function\n')
        elif self.error == 'redefinition_err':
            sys.stderr.write('redefinition variable\n')
        elif self.error == 'incorrect_fe_err':
            sys.stderr.write('incorrect finite element\n')
        elif self.error == 'read_file_err':
            sys.stderr.write('read file error\n')
        elif self.error == 'unknown_fe_err':
            sys.stderr.write('unknown finite element type\n')
        elif self.error == 'solve_method_err':
            sys.stderr.write('not specified method for solving linear systems\n')
        elif self.error == 'problem_type_err':
            sys.stderr.write('unknown problem type (static or dynamic)\n')
        elif self.error == 'elasticity_err':
            sys.stderr.write('not specified the parameters of elasticity\n')
        elif self.error == 'time_err':
            sys.stderr.write('incorrectly specified time parameters\n')
        else:
            sys.stderr.write(self.error + '\n')
