# This file is part of the FEPX software package.
# Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
# See the COPYING file in the top-level directory.

# Sphinx configuration of the FEPX documentation (ROADMAP task M.1).
#
# Build:  make html (in doc/), or sphinx-build -W -b html doc <output directory>
# Needs:  Sphinx (python3-sphinx) and, for the look of the website, pydata-sphinx-theme
#         (python3-pydata-sphinx-theme); without the theme the built-in "alabaster" theme is used.

project = 'FEPX'
author = 'DPLab, ACME Lab, CNRS'
copyright = 'DPLab, ACME Lab, CNRS'
version = u'3.0.0'
release = u'3.0.0'
language = 'en'

# reStructuredText sources: every page is listed in a toctree (no orphans).
source_suffix = {'.rst': 'restructuredtext'}
master_doc = 'index'
exclude_patterns = ['_build', 'build']

# Equations are written in LaTeX (:math:) and rendered with MathJax.
extensions = ['sphinx.ext.mathjax']
# Every reference must resolve: a keyword (confval) with no entry, or a broken link, is an
# error under -W (step 9.9, D104).
nitpicky = True
# The keyword entries (confval) are not listed in the table of contents.
toc_object_entries = False

# Numbered figures and tables, so pages can refer to them.
numfig = True

# Theme: the website's (pydata-sphinx-theme, light mode only) when installed, otherwise
# Sphinx's default.
try:
    import pydata_sphinx_theme  # noqa: F401
    html_theme = 'pydata_sphinx_theme'
    html_logo = '_static/fepx-logo-text.svg'
    html_theme_options = {
        'logo': {'alt_text': 'FEPX'},
        'navbar_align': 'left',
        'navbar_center': [],
        'navbar_end': [],
        'show_nav_level': 1,
        'navigation_depth': 5,
        'show_toc_level': 2,
        'secondary_sidebar_items': ['page-toc'],
        'navbar_start': ['navbar-logo', 'fepx-version'],
        'footer_start': ['copyright'],
        'footer_end': ['theme-version'],
    }
    html_context = {'default_mode': 'light'}
    # Every page of the documentation is a top-level page, which the theme would list in the header
    # and leave the section navigation empty: list them in the left sidebar instead, as on the
    # website (where they are under "Documentation").
    html_sidebars = {'**': ['fepx-doc-nav.html']}
    templates_path = ['_templates']
except ImportError:
    html_theme = 'alabaster'

html_favicon = '_static/fepx-favicon.png'
html_static_path = ['_static']
html_css_files = ['custom.css']
html_title = 'FEPX documentation'


def setup(app):
    # Light mode only, whatever a previous visit saved in the browser
    app.add_js_file(None, body="""
      try {
        localStorage.setItem("mode", "light");
        localStorage.setItem("theme", "light");
      } catch (e) {}
      document.documentElement.dataset.mode = "light";
      document.documentElement.dataset.theme = "light";
    """)
