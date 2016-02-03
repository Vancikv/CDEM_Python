'''
Created on 18. 3. 2014

@author: Vancikv
'''

import os

def get_outfile(folder_path_list, file_name):
    '''Returns a file in the specified folder using the home
    directory as root.
    '''
    HOME_DIR = os.path.expanduser("~")
    out_dir = os.path.join(HOME_DIR, *folder_path_list)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    outfile = os.path.join(out_dir, file_name)
    return outfile
