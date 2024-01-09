# -*- coding: utf-8 -*-

import sys
import os
import re


# Prefer to use the version of the theme in this repo
# and not the installed version of the theme.
sys.path.insert(0, os.path.abspath('.'))

import sphinx_rtd_theme

project = u'FEPX'
version = u'2.0.0'
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

def setup(app):
    app.add_css_file('my_theme.css')

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
