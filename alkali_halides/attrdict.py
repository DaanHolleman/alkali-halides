# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:35:06 2024

@author: dholl
"""

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
