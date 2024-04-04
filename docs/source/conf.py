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
import inspect
sys.path.insert(0, os.path.abspath('.'))

__location__ = os.path.join(os.getcwd(), os.path.dirname(
    inspect.getfile(inspect.currentframe())))

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# add 'ray-optics/src' to path
sys.path.insert(0, os.path.join(__location__, '../../src'))

# -- Project information -----------------------------------------------------

project = u'ray-optics'
copyright = u'2017-2024, Michael J. Hayford'
author = u'Michael J. Hayford'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = ''
# The full version, including alpha/beta/rc tags.
release = ''

try:
    from rayoptics import __version__ as version
except ImportError:
    pass
else:
    release = version

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.intersphinx', 
              'sphinx.ext.todo', 'sphinx.ext.autosummary', 
              'sphinx.ext.viewcode', 'sphinx.ext.coverage',
              'sphinx.ext.doctest', 'sphinx.ext.ifconfig', 
              'sphinx.ext.mathjax', 'sphinx.ext.napoleon']

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['Thumbs.db', '.DS_Store', 'tests']

# The master toctree document.
master_doc = 'index'

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False

autodoc_member_order = 'bysource'

# The name of the Pygments (syntax highlighting) style to use.
#pygments_style = 'sphinx'
pygments_style = 'friendly'
# pygments_style = 'xcode'
# pygments_style = 'solarize-light'

rst_prolog = """
.. |ybar| replace:: :math:`y-\overline{y}`
.. |nubar| replace:: :math:`\omega-\overline{\omega}`
.. |minimum_python_version| replace:: 3.10
.. |minimum_numpy_version| replace:: 1.24.4
.. |Series| replace:: :class:`pandas.Series`
.. |DataFrame| replace:: :class:`pandas.DataFrame`
"""

# A list of ignored prefixes for module index sorting.
modindex_common_prefix = ['rayoptics.']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'agogo'
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 5,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- External mapping ------------------------------------------------------------
python_version = '.'.join(map(str, sys.version_info[0:2]))
intersphinx_mapping = {
    'sphinx': ('https://www.sphinx-doc.org/en/master', None),
    'python': ('https://docs.python.org/' + python_version, None),
    'matplotlib': ('https://matplotlib.org/stable', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'opticalglass': ('https://opticalglass.readthedocs.io/en/latest', None),
}
