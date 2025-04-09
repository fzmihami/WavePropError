# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0, os.path.abspath('../../')) 

project = 'WavePropError'
copyright = '2025, Fatima-Zahra Mihami'
author = 'Fatima-Zahra Mihami'
release = 'v1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosummary',
              'sphinx.ext.graphviz',
              'sphinxcontrib.bibtex',
              'sphinx.ext.mathjax',
              'sphinx_autodoc_typehints',
              'nbsphinx',]


# C++ Documentation
extensions += ['breathe']
breathe_projects = {
    "WavePropError": "../../docs/doxygen_output/xml"
}
breathe_default_project = "WavePropError"

templates_path = ['_templates']
exclude_patterns = []
bibtex_bibfiles = ["references.bib"]  

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Options for autodoc ----------------------------------------------------
autodoc_default_options = {
    'members': True,  # Include all members of a class
    'undoc-members': True,  # Include undocumented members
    'show-inheritance': True,  # Show inheritance in class documentation
}