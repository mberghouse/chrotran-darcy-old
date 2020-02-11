#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:40:03 2020

@author: rleone
"""

# STEPS TO CREATE KEYWORD DATABASE
# 1. Run regression tests:
#   > cd $PFLOTRAN_DIR/regression_tests
#   > make test
# 2. Create a folder called 'docs' in keyword_database folder
#   > cd $PFLOTRAN_DIR/src/python/keyword_database
#   > mkdir docs
# 3. Run script
#   > python createdatabase.py
    
    
# STEPS TO CREATE HTML FILE
# 1. Go into docs directory
#   > cd docs
# 2. Setup sphinx (will give error that there is already an index file but 
#    ignore this)
#   > sphinx-quickstart
# 3. Make the html files
#   > make clean
#   > make html
# 4. Open up html files
#   > cd _build/html
#     open index.html 

import os
import numpy as np
import re

base_dir = os.getcwd()

os.chdir('../../../regression_tests')

files = []
for root, dirnames, filenames in os.walk(os.getcwd()):
    for filename in filenames:
        if filename.endswith('.stdout'):
            files.append(os.path.join(root,filename))
            

keywords_dict = {}

for file in files:
#    keywords_dict[file] = []
    with open(file,'r') as f:
        for line in f:
            if line.startswith(' KEYWORD'):
 #               print(re.split('(:|,|\s)\s*',line))
                words = re.split(': |,|\n',line.strip())
                
                
                while not words[-1].isupper(): ##dangerous???
                    words.pop(-1)
                
                keyword = words[-1]
                try:
                    previous_keyword = words[-2]
                except:
                    previous_keyword = None
                    
                try:
                    second_previous_keyword = words[-3]
                except:
                    second_previous_keyword = None
                    
                #####remove what isn't actually a keyword
                if previous_keyword == 'CHARACTERISTIC_CURVES' or previous_keyword == 'SATURATION_FUNCTIONS':
                    keyword = previous_keyword
                if previous_keyword == 'CHEMISTRY':# or 'CYBERNETIC':
                    if re.search(r'\d',keyword):
                        keyword = previous_keyword
                    if re.search('-',keyword):
                        keyword = previous_keyword
                    if re.search('\+',keyword):
                        keyword = previous_keyword
                if previous_keyword == 'CYBERNETIC':
                    if re.search(r'\d',keyword):
                        keyword = previous_keyword
                if previous_keyword == 'MINERALS':
                    keyword = previous_keyword
                if previous_keyword == 'CONCENTRATIONS':
                    keyword = previous_keyword
                if previous_keyword == 'FREE_ION_GUESS':
                    keyword = previous_keyword
                if previous_keyword == 'PREFACTOR_SPECIES':
                    keyword = previous_keyword
                if second_previous_keyword=='CHEMISTRY_2ND_PASS':
                    regex = re.compile('/+')
                    if re.search(r'\d',keyword):
                        keyword = second_previous_keyword
                    if re.search('-',keyword):
                        keyword = second_previous_keyword
                    if re.search('\+',keyword):
                        keyword = second_previous_keyword

                    
                #####Store keywords
                if keyword in keywords_dict.keys():
                    if file not in keywords_dict[keyword]:
                        keywords_dict[keyword].append(file)
                else:
                    keywords_dict[keyword] = [file]
#                if len(words)>1 and words[-1] not in keywords_dict[file]:
#                    keywords_dict[file].append(words[-1])


with open('{}/docs/index.rst'.format(base_dir),'w') as f:

    intro = """
*************************
PFLOTRAN Keyword Database
*************************

.. toctree::
   :maxdepth: 2
""" 
        
    f.write(intro)

    for keyword,filepaths in sorted(keywords_dict.items()):
        keyword_rst = keyword + '.rst'
        f.write("""
   {}""".format(keyword_rst))
    
        with open('{}/docs/{}'.format(base_dir,keyword_rst),'w') as fin:

          fin.write('.. _{}:'.format(keyword))
          fin.write('\n\n')
          fin.write('{}'.format(keyword))
          fin.write('\n')
          fin.write('{}\n'.format('='*len(keyword)))
          fin.write('\n')
          for filepath in filepaths:
            file = re.split('/|\.',filepath.strip())[-2]
            fin.write("""
{}
            """.format(file))

            
