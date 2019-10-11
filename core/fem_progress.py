#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################
#         Реализация отображения прогресса при вычислениях
###################################################################

import time
import sys
import itertools
from threading import Thread


class TProgress:
    def __init__(self):
        self.__process_id__ = ''
        self.__process_start__ = 0
        self.__process_stop__ = 0
        self.__process_current__ = 0
        self.__process_old__ = 0

    def set_process(self, pid, start, stop):
        self.__process_id__ = pid
        self.__process_start__ = start
        self.__process_stop__ = stop
        self.__process_old__ = 0
        if start <= stop:
            sys.stdout.write('\r' + self.__process_id__ + ' 0%')
            sys.stdout.flush()

    def set_progress(self, current):
        self.__process_current__ = current
        pos = int((100.0*float(self.__process_current__))/float(self.__process_stop__ - self.__process_start__ + 1))
        if pos == self.__process_old__ or not pos:
            return
        if pos % 2 == 0:
            sys.stdout.write('\r' + self.__process_id__ + ' ' + str(pos) + '%')
            if pos == 100:
                sys.stdout.write('\n')
        sys.stdout.flush()
        self.__process_old__ = pos


class TThreadProgress(Thread):
    def __init__(self, name):
        Thread.__init__(self)
        self.is_stopped = False
        self.progress_name = name

    def run(self):
        for c in itertools.cycle(['|', '/', '-', '\\']):
            if self.is_stopped:
                break
            sys.stdout.write('\r' + self.progress_name + ' ' + c)
            sys.stdout.flush()
            time.sleep(0.1)
        sys.stdout.write('\r' + self.progress_name + ' 100%\n')

    def stop(self):
        self.is_stopped = True
