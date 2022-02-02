# -*- coding: utf-8 -*-

import sys
import os
import re


# Prefer to use the version of the theme in this repo
# and not the installed version of the theme.
sys.path.insert(0, os.path.abspath('.'))

import sphinx_rtd_theme

project = u'FEPX'
version = u'1.3.0'
release = version
author = u'ACME Lab'
copyright = u'ACME Lab'
language = 'en'

master_doc = 'index'

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': False,
    'navigation_depth': 5,
}
