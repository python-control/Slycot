# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../slycot'))

project = 'Slycot'
copyright = '2023, Slycot Developers'
author = 'Slycot Developers'
release = '0.6'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.todo', 'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx', 'sphinx.ext.imgmath',
    'sphinx.ext.autosummary', 'nbsphinx', 'numpydoc',
    'sphinx.ext.doctest'
]
# scan documents for autosummary directives and generate stub pages for each.
autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
