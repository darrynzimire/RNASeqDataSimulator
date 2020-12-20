#!/usr/bin/python3
# encoding = UTF-8

import random
import numpy as np
import scipy as sp

SEED = 12345


def negative_binomial():
    
    np.random.seed(SEED)
    NB = np.round(sp.random.negative_binomial(n=1, p=0.1, size=100000))
    NB_model = [1 if x == 0 else x for x in NB]
    
    return NB_model


def Uniform_dist(low, high=None, size=None):
    """
    Description
    This function generates random numbers according to a specified distribution
    Parameters
    low (int): The lowest value to sample
    high (int) the upper bound for sampling
    Return: Returns a vector of numbers drawn from a user-specified distribution
    """
    np.random.seed(SEED)
    counts = np.round(np.random.uniform(low=0, high=high, size=size))

    return counts



def Poisson():
    pass


