# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 10:36:42 2024

@author: dholl
"""

from numpy import log10

def leading_zeros(x:list):
    """
    Calculates how many leading zeros should be present for lower numbers given a maximum number.
    """
    N = x if type(x) in [int,bool] else len(x)
    if N == 0: N = 1
    return int(log10(N))+1

def select(lst:list, chunksize:int = 20):
    """
    Select entries from a list via user feedback.
    """
    # Pick first item if length is 1
    if len(lst) == 1:
        return lst[0]
    
    # Raise error for an empty list since nothing can be selected
    if len(lst) == 0:
        raise Exception('Empty list!')
    
    # Print options to pick from if the list is longer than 1s
    chunks = [[]]
    lst_iter = enumerate(lst)
    for index, item in lst_iter:
        for c in range(chunksize):
            chunks[-1] += [(index, item)]
            try:
                index, item = next(lst_iter)
            except Exception as e:
                if type(e) is StopIteration:
                    break
                raise Exception(e)
        else:
            chunks += [[]]
    tot_chunks = len(chunks)
    l10 = leading_zeros(lst)
    
    # Loop through pages of possible choices
    print('Please select an item.')
    for jj, chunk in enumerate(chunks):
        print(f'PAGE {jj+1}/{tot_chunks}')
        for ii, item in chunk:
            print(f'[{ii:0{l10}d}] \t{item}')
        user = input('>>> ').strip()
        if len(user) > 0:
            return lst[int(user)]
    
    # If no choice was made, keep asking for input
    while len(user) == 0:
        user = input('>>> ').strip()
        if len(user) > 0:
            return lst[int(user)]
    
    # Turn back! You should not be here!
    raise Exception('You got somewhere illegal!')
