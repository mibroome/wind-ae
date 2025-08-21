# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# sys.path.insert(0, os.path.abspath('../wind_ae/'))
# sys.path.insert(0, os.path.abspath('../wind_ae/wrapper/'))
# sys.path.insert(0, os.path.abspath('../wind_ae/wrapper/wrapper_utils/'))


# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Wind-AE'
copyright = '2025, Madelyn Broome, John McCann, Ruth Murray-Clay'
author = 'Madelyn Broome, John McCann, Ruth Murray-Clay'
release = '0.2.0'
root_doc = 'index'
autodoc_member_order = 'bysource'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "nbsphinx"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
