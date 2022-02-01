# -*- coding: utf-8 -*-

import sys
import os
import re

# Prefer to use the version of the theme in this repo
# and not the installed version of the theme.
sys.path.insert(0, os.path.abspath('.'))

import sphinx_rtd_theme

project = u'FEPX'
version = u'1.2.2-90-devel'
release = u'1.2.2-90-devel'
author = u'Matthew Kasemer'
copyright = u''
language = 'en'

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': False,
    'navigation_depth': 5,
}
