# -*- coding: utf-8 -*-

import sys
import os
import re

# The manual's settings (docs/manual/conf.py of neperfepx/fepx3, step 9.21): mathematics,
# numbered figures, no table-of-contents entries for objects, its style sheet
extensions = ['sphinx.ext.mathjax', 'sphinx_design', 'sphinx_copybutton', 'sphinxext.opengraph']
numfig = True
toc_object_entries = False
exclude_patterns = ['build']

# Copy commands without the prompt, and keep the lines of multi-line commands
copybutton_prompt_text = '$ '
copybutton_line_continuation_character = '\\'

ogp_site_url = 'https://fepx.info/'
ogp_site_name = 'FEPX'
ogp_image = 'https://fepx.info/imgs/home-fepx.png'
ogp_social_cards = {'enable': False}

project = u'FEPX'
version = u'3.0.0'
release = version
author = u'ACME Lab, CNRS'
copyright = u'ACME Lab, CNRS'
language = 'en'

master_doc = 'index'

templates_path = ['_templates']

html_theme = 'pydata_sphinx_theme'
html_title = 'FEPX'
html_logo = '_static/fepx-logo-text.svg'
html_favicon = '_static/fepx-favicon.png'
html_theme_options = {
    'logo': {'alt_text': 'FEPX'},
    'navbar_align': 'left',
    'navbar_end': [],
    'header_links_before_dropdown': 8,
    'show_nav_level': 1,
    'navigation_depth': 5,
    'show_toc_level': 2,
    'secondary_sidebar_items': ['page-toc'],
    'navbar_start': ['navbar-logo', 'fepx-version'],
    'footer_start': ['copyright'],
    'footer_end': ['theme-version'],
}

# Light mode only: the images have white backgrounds
html_context = {'default_mode': 'light'}

# These pages have no section navigation to show
html_sidebars = {page: [] for page in ('index', 'gallery', 'papers', 'applications',
                                       'downloads', 'whatsnew', 'community')}

html_css_files = ['custom.css']

# Field tags of the application papers (applications.rst), written as :tag:`key`
APPLICATION_TAGS = {
    'xrd': 'diffraction experiments',
    'meth': 'methods & formulations',
    'proc': 'processing & texture',
    'ti': 'hexagonal metals',
    'frac': 'fatigue, cyclic loading & damage',
    'other': 'other',
}

def tag_role(name, rawtext, text, lineno, inliner, options=None, content=None):
    from docutils import nodes
    key = text.strip()
    if key not in APPLICATION_TAGS:
        msg = inliner.reporter.error('Unknown application tag "%s" (see APPLICATION_TAGS in conf.py)'
                                     % key, line=lineno)
        return [inliner.problematic(rawtext, rawtext, msg)], [msg]
    html = '<span class="neper-tag" data-tag="%s">%s</span>' % (key, APPLICATION_TAGS[key])
    return [nodes.raw('', html, format='html')], []

def setup(app):
    app.add_role('tag', tag_role)
    app.add_js_file('applications.js')
    app.add_js_file('applications-chart.js')
    app.add_css_file('my_theme.css')
    # Light mode only, whatever a previous visit saved in the browser
    app.add_js_file(None, body="""
      try {
        localStorage.setItem("mode", "light");
        localStorage.setItem("theme", "light");
      } catch (e) {}
      document.documentElement.dataset.mode = "light";
      document.documentElement.dataset.theme = "light";
    """)

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
