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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------


project = 'torch-adata'
copyright = '2022, Michael E. Vinyard'
author = 'Michael E. Vinyard'

# The full version, including alpha/beta/rc tags


release = '0.0.16'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = ["myst_parser"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Generate API documentation ----------------------------------------------------------

extensions = ['autoapi.extension']
autoapi_dirs = ['../../torch_adata']
autosummary_generate = True

html_static_path = ['_static']
templates_path = ["_templates"]


# -- Options for HTML output -------------------------------------------------------------

html_show_sourcelink = True
html_theme = 'pydata_sphinx_theme'
html_favicon = '../imgs/torch-adata.logo.small.svg'

html_context = dict(
    github_user="mvinyard",   # Username
    github_repo="torch-adata",   # Repo name
    github_version="master",  # Version
    doc_path="docs/",  # Path in the checkout to the docs root
)

# Set link name generated in the top bar.
html_title = "torch-adata"
html_logo = "../imgs/torch-adata.logo.small.svg"

html_theme_options = {
    "github_url": "https://github.com/mvinyard/torch-adata",
    "twitter_url": "https://twitter.com/vinyard_m",
    "codetextcolor": "ee4c2c",
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

