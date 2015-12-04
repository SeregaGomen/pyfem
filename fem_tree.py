#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
# Реализация дерева разбора арифметических и логических выражений
###################################################################

import math
from abc import abstractmethod



# Абстрактный базовый класс значения выражения
class TNode:
    @abstractmethod
    def value(self):
        pass


# Класс, реализующий дерево разбора арифметических выражений
class TTree:
    def __init__(self, *args):
        if len(args) == 0:
            self.node = TRealNode(0)
        elif len(args) == 1:
            self.node = TRealNode(args[0])
        elif len(args) == 2:
            self.node = TUnaryNode(args[0], args[1])
        elif len(args) == 3:
            self.node = TBinaryNode(args[0], args[1], args[2])

    def value(self):
        return self.node.value()


# Вещественная переменная
class TRealNode(TNode):
    def __init__(self, val):
        self.val = val

    def value(self):
        return self.val


# Унарная операция
class TUnaryNode(TNode):
    def __init__(self, op, val):
        self.op = op
        self.val = val

    def value(self):
        if self.op == '-':
            return -self.val.value()
        elif self.op == '+':
            return +self.val.value()
        elif self.op == 'abs':
            return math.fabs(self.val.value())
        elif self.op == 'sin':
            return math.sin(self.val.value())
        elif self.op == 'cos':
            return math.cos(self.val.value())
        elif self.op == 'tan':
            return math.tan(self.val.value())
        elif self.op == 'exp':
            return math.exp(self.val.value())
        elif self.op == 'asin':
            return math.asin(self.val.value())
        elif self.op == 'acos':
            return math.acos(self.val.value())
        elif self.op == 'atan':
            return math.atan(self.val.value())
        elif self.op == 'sinh':
            return math.sinh(self.val.value())
        elif self.op == 'cosh':
            return math.cosh(self.val.value())
        elif self.op == 'not':
            return 1 if self.val.value() == 0 else 0


# Бинарная операция
class TBinaryNode(TNode):
    def __init__(self, left, op, right):
        self.left = left
        self.op = op
        self.right = right

    def value(self):
        if self.op == '+':
            return self.left.value() + self.right.value()
        elif self.op == '-':
            return self.left.value() - self.right.value()
        elif self.op == '*':
            return self.left.value()*self.right.value()
        elif self.op == '/':
            return self.left.value()/self.right.value()
        elif self.op == '^':
            return math.pow(self.left.value(), self.right.value())
        elif self.op == '=':
            return 1 if self.left.value() == self.right.value() else 0
        elif self.op == '<>':
            return 0 if self.left.value() == self.right.value() else 1
        elif self.op == '<':
            return 1 if self.left.value() < self.right.value() else 0
        elif self.op == '<=':
            return 1 if self.left.value() <= self.right.value() else 0
        elif self.op == '>':
            return 1 if self.left.value() > self.right.value() else 0
        elif self.op == '>=':
            return 1 if self.left.value() >= self.right.value() else 0
        elif self.op == 'or':
            return 0 if (self.left.value() or self.right.value()) == 0 else 1
        elif self.op == 'and':
            return 0 if (self.left.value() and self.right.value()) == 0 else 1
        elif self.op == 'atan2':
            return math.atan2(self.left.value(), self.right.value())
