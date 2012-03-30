'''
Created on 30 Mar 2012

@author: garyt




'''
import sys
import os
import glob

def check_args():
    
    args = []
    if len(sys.argv) > 1:
        args = sys.argv[1:]
    if len(args) != 1:
        template = "Camshift_table_builder requires a single directory as its argument but got '%s'"
        print >> sys.stderr, template % `args`
        print >> sys.stderr, "exiting..."
        sys.exit(-1)
        
def parameter_file_iterator(table_dir, template = '1[A-Z]*orgpairs.par*'):
    if not os.path.isabs(table_dir):
        table_dir = os.path.join(os.getcwd(), table_dir)
    query = os.path.join(table_dir, template)
    for filename in glob.glob(query):
        yield filename

class Camshift_table_builder(object):
    '''
    classdocs
    '''


    def __init__(self,params):
        '''
        Constructor
        '''
        



if __name__ == '__main__':
    check_args()
    
    table_dir = sys.argv[1]

    for par_filename in parameter_file_iterator(table_dir):
        print par_filename
    