# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os
import sys

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
morphman_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "morphman")

sys.path.insert(0, morphman_root)
sys.path.insert(0, os.path.join(morphman_root, "common"))
sys.path.insert(0, os.path.join(morphman_root, "misc"))

# -- Project information -----------------------------------------------------

project = 'morphman'
copyright = '2022, Aslak W. Bergersen & Henrik A. Kjeldsberg'
author = 'Aslak W. Bergersen & Henrik A. Kjeldsberg'

# The short X.Y version.
version = '1.3'
# The full version, including alpha/beta/rc tags
release = '1.3'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.githubpages'
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

autodoc_mock_imports = ["numpy", "scipy", "vtk", "vmtk", "morphman"]

# -- Options for HTML output -------------------------------------------------

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Output file base name for HTML help builder.
htmlhelp_basename = 'morphMan_docs_help'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# -- Options for LATEX output ------------------------------------------------
# latex_engine = 'pdflatex'
#
# latex_elements = {}
#
# latex_documents = [
#     (master_doc, 'morphMan_docs.tex', 'morphman\\_docs Documentation',
#      'Aslak W. Bergersen & Henrik A. Kjeldsberg', 'manual'),
# ]
#
# man_pages = [
#     (master_doc, 'morphMan_docs', u'morphMan_docs Documentation',
#      [author], 1)
# ]
# texinfo_documents = [
#     (master_doc, 'morphMan_docs', u'morphman Documentation',
#      author, 'morphMan_docs', 'One line description of project.',
#      'Miscellaneous'),
# ]

# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/': None}
