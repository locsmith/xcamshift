#!/usr/bin/python

import sys

for line in sys.stdin:
     fields =  line.split()
     for field in fields:
          if field.startswith('test_data'):
                print field
