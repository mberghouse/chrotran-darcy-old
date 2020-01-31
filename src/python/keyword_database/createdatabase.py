#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:40:03 2020

@author: rleone
"""

#######MUST RUN REGRESSIONS FIRST BEFORE RUNNING PYTHON PROGRAM#######
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
                if keyword in keywords_dict.keys():
                    if file not in keywords_dict[keyword]:
                        keywords_dict[keyword].append(file)
                else:
                    keywords_dict[keyword] = [file]
#                if len(words)>1 and words[-1] not in keywords_dict[file]:
#                    keywords_dict[file].append(words[-1])

###WORKING 
#with open('{}/keyword_database/index.rst'.format(root_dir),'w') as f:
#
#    intro = """
#*************************
#PFLOTRAN Keyword Database
#*************************
#
#.. toctree::
#   :maxdepth: 2
#""" 
#        
#    f.write(intro)
#
#    for keyword,filepaths in keywords_dict.items():
#        keyword_rst = keyword + '.rst'
#        f.write("""
#   {}""".format(keyword_rst))
#    
#    
#        with open('{}/keyword_database/{}'.format(root_dir,keyword_rst),'w') as fin:
#          fin.write('.. _{}:'.format(keyword))
#          fin.write('\n\n')
#          fin.write('{}'.format(keyword))
#          fin.write('\n')
#          fin.write('{}\n'.format('='*len(keyword)))
#          fin.write('\n')
#          for filepath in filepaths:
#            file = re.split('/|\.',filepath.strip())[-2]
#            fin.write("""
#{}
#            """.format(file))
            
##TEST SORT
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
#f = open('{}/keyword_database/index.rst'.format(root_dir),'w')  
#
#intro = """
#*************************
#PFLOTRAN Keyword Database
#*************************
#
#.. toctree::
#   :maxdepth: 2
#""" 
#        
#f.write(intro)
            
#for keyword,filepaths in keywords_dict.items():
#    keyword_rst = keyword + '.rst'
#    f.write("""
#   {}""".format(keyword_rst))
#    
#    
#    with open('{}/keyword_database/{}'.format(root_dir,keyword_rst),'w') as fin:
#        fin.write('.. _{}:'.format(keyword))
#        fin.write('\n\n')
#        fin.write('{}'.format(keyword))
#        fin.write('\n')
#        fin.write('{}\n'.format('='*len(keyword)))
#        fin.write('\n')
#        for filepath in filepaths:
#            file = re.split('/|\.',filepath.strip())[-2]
#            fin.write("""
#{}
#            """.format(file))
    
                    
                    
                    
#################################
#Dictionary based on regress test
#################################            
#keywords_dict = {}
#
#for file in files:
#    keywords_dict[file] = []
#    with open(file,'r') as f:
#        for line in f:
#            if line.startswith(' KEYWORD'):
# #               print(re.split('(:|,|\s)\s*',line))
#                words = re.split(': |,|\n',line.strip())
#                if not words[-1].isupper():
#                    words.pop(-1)
#                
#                if len(words)>1 and words[-1] not in keywords_dict[file]:
#                    keywords_dict[file].append(words[-1])
                    
                    
                
#                print(line.strip().split(','))#[0].split(':'))
   
####Write index file etc     
#f = open('{}/keyword_database/index.rst'.format(root_dir),'w')  
#
#intro = """
#*************************
#PFLOTRAN Keyword Database        
#*************************
#
#.. toctree::
#   :maxdepth: 2
#""" 
#        
#f.write(intro)
#
#for filepath,keywords in keywords_dict.items():
#    file = re.split('/|\.',filepath.strip())[-2] #filename.strip().split('/')[-1]
#    file_rst = re.split('/|\.',filepath.strip())[-2] + '.rst'
#    f.write('{}\n'.format(file_rst))
#    f.write("""
#   {}""".format(file_rst))
    
    
#    with open('{}/keyword_database/{}'.format(root_dir,file_rst),'w') as fin:
#        fin.write('.. _{}:'.format(file))
#        fin.write('\n\n')
#        fin.write('{}'.format(file))
#        fin.write('\n')
#        fin.write('{}\n'.format('='*len(file)))
#        fin.write('\n')
#        for keyword in keywords:
#            fin.write("""
#{}
#            """.format(keyword))
#    print(keywords_dict[file][3])
    
            