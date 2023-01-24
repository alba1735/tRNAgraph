#!/usr/bin/env python3

import os

def builder(directory):
    '''
    Function to create output directory if it does not already exist
    '''
    if not os.path.exists(directory):
        print('Creating output directory: {}'.format(directory))
        os.makedirs(directory, exist_ok=True)
    else:
        print('Output directory already exists: {}'.format(directory))

def path_cleaner(directory):
    '''
    Function to clean up input path
    '''
    if directory[-1] == '/':
        directory = directory[:-1]
    return directory