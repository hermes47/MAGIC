#!/usr/bin/env python3

import os

branches = ('MCW5','FRFY','FRFYFRFY','MCW5FRFY','FRFYMCW5','','MCW5MCW5')

for steam in branches:
    os.system('./magic.py -s "P8EO(B16H%s)"' %steam)
