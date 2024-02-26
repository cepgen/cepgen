# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = u'CepGen'
copyright = u'2013-@CURRENT_YEAR@, the CepGen Collaboration'
author = u'Laurent Forthomme'

# The short X.Y version
version = u'@CEPGEN_VERSION@'
# The full version, including alpha/beta/rc tags
release = u'@CEPGEN_VERSION@'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'breathe',
    'changelog',
    'sphinx_git',
    #'sphinxemoji.sphinxemoji',
    'sphinxcontrib.bibtex',
    'sphinx_togglebutton',
    'sphinx_toolbox.collapse',
    'sphinx_math_dollar',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
]

bibtex_bibfiles = ['_static/bibliography.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [u'_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
#pygments_style = 'colorful'
pygments_style = None


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_logo = 'small-cepgen-logo.png'

#import pietroalbini_sphinx_themes
#html_theme_path = [pietroalbini_sphinx_themes.themes_path()]
#html_theme = "pietroalbini"

#import sphinx_redactor_theme
#html_theme = 'sphinx_redactor_theme'
#html_theme_path = [sphinx_redactor_theme.get_html_theme_path()]

#html_theme = 'classic'
#html_theme = 'sphinx_rtd_theme'
#html_theme = 'pyramid'
#html_theme = 'haiku'
#html_theme = 'traditional'
#html_theme = 'alabaster'
#html_theme = 'karma_sphinx_theme'
html_theme = 'furo'

#import solar_theme
#html_theme = 'solar_theme'
#html_theme_path = [solar_theme.theme_path]

#import hachibee_sphinx_theme
#html_theme = 'hachibee'
#html_theme_path = [hachibee_sphinx_theme.get_html_themes_path()]

#import kotti_docs_theme
#html_theme = 'kotti_docs_theme'
#html_theme_path = [kotti_docs_theme.get_theme_dir()]

#import guzzle_sphinx_theme
#html_theme = 'guzzle_sphinx_theme'
#html_theme_path = guzzle_sphinx_theme.html_theme_path()
#extensions.append("guzzle_sphinx_theme")

html_extra_path = ['.htaccess']

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/cepgen/cepgen",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
    #'display_version': False,
    #'sidebarwidth': '300px',
    #'style_external_links': True,
    #'logo_only': True,
}
#html_use_index = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'CepGendoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': r'',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'CepGen.tex', u'CepGen Documentation',
     u'Laurent Forthomme', 'manual'),
]

mathjax3_config = {
  'tex': {
    'inlineMath': [["\\(","\\)"]],
    'displayMath': [["\\[","\\]"]],
    'processEscapes': True,
    'macros': {
      'Pom': "{\\rm I\\!P}",
      'Reg': "{\\rm I\\!R}",
      'gg': ["{\\gamma\\gamma\\rightarrow #1}", 1],
      'ggx': "{\\gg{X}}",
      'ggll': "{\\gg{\ell^+\\ell^-}}",
      'ggff': "{\\gg{f\\bar f}}",
      'ggww': "{\\gg{W^+W^-}}",
      'kt': "{k_{\\rm T}}",
      'pt': "{p_{\\rm T}}",
      'vecqt': "{\\bf q_{\\rm T}}",
      'xbj': "{x_{\\rm Bj}}",
    }
  },
}


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'cepgen', u'CepGen Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'CepGen', u'CepGen Documentation',
     author, 'CepGen', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

# Breathe Configuration
breathe_default_project = "CepGen"
breathe_implementation_filename_extensions = ['.cxx', '.C', '.f']

# Changelog configuration
#changelog_render_changeset = "https://phab.hepforge.org/rCEPGEN%s"
#changelog_render_pullreq = "https://phab.hepforge.org/D%s"
changelog_render_changeset = "https://github.com/cepgen/cepgen/commit/%s"
changelog_render_pullreq = "https://gitlab.cern.ch/lforthom/cepgen/-/merge_requests/%s"

def setup(app):
    app.add_css_file('hacks.css')
#    app.add_js_file('mathconf.js', type='text/x-mathjax-config')
#    app.add_js_file('https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_HTMLorMML', async=True)