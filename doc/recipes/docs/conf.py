#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Fortnet Recipes documentation build configuration file, created by
# sphinx-quickstart on Thu Jun 25 20:16:11 2021.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Fortnet'
copyright = '2020, T. W. van der Heide'
author = 'T. W. van der Heide'

# The full version, including alpha/beta/rc tags
release = '0.5'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
needs_sphinx = '1.8'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.mathjax', 'sphinxcontrib.bibtex']
bibtex_bibfiles = ['refs.bib']
bibtex_encoding = 'utf-8-sig'
bibtex_default_style = 'unsrt'
bibtex_reference_style = 'label'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'Fortnet Recipes'
copyright = '2022, T. W. van der Heide'
author = 'T. W. van der Heide'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '0.5'
# The full version, including alpha/beta/rc tags.
release = '0.5'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Figures are enumerated and can be reference by the :numref: directive
numfig = True

numfig_format = {
    'figure': 'Figure %s',
    'table': 'Table %s',
    'code-block': 'Listing %s',
    'section': 'Section'
}


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_logo = '../../../utils/art/logo_doc.svg'
html_theme_options = {
    'logo_only': False,
    'display_version': True,
}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
#html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'FortnetRecipesdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'FortnetRecipes.tex', 'Fortnet Recipes Documentation',
     'T. W. van der Heide', 'manual'),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'fortnetrecipes', 'Fortnet Recipes Documentation',
     [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'FortnetRecipes', 'Fortnet Recipes Documentation',
     author, 'FortnetRecipes', 'One line description of project.',
     'Miscellaneous'),
]



# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Additional functionality to archive files in _downloads

def setup(app):
    app.connect('builder-inited', create_archives)


def create_archives(app):
    '''Creates a .tar.bz2 archive in _downloads for each entry in _archives.'''

    import os
    import tarfile
    from sphinx.util import logging

    ARCHIVE_DIR = '_archives'
    DOWNLOAD_DIR = '_downloads/archives'
    IGNORE_FILE = '_ignore'
    COMPRESSION = 'bz2'

    logger = logging.getLogger(__name__)
    cwd = os.getcwd()
    archivedir = os.path.join(cwd, ARCHIVE_DIR)
    downloaddir = os.path.join(cwd, DOWNLOAD_DIR)
    if not os.path.exists(downloaddir):
        os.makedirs(downloaddir)
    archives = [fname for fname in os.listdir(archivedir) if fname[0] != '_']
    tarfilter = _get_tar_filter(os.path.join(archivedir, IGNORE_FILE))
    for archive in archives:
        tarname = os.path.join(DOWNLOAD_DIR, archive + '.tar.' + COMPRESSION)
        tarpath = os.path.join(cwd, tarname)
        filename = os.path.join(ARCHIVE_DIR, archive)
        filepath = os.path.join(cwd, filename)
        tartime = os.path.getmtime(tarpath) if os.path.exists(tarpath) else 0.0
        filetime = _newest_modification_time(filepath)
        if tartime <= filetime:
            with tarfile.open(tarpath, 'w:' + COMPRESSION) as tar:
                tar.add(filepath, arcname=archive, filter=tarfilter)
            logger.info('Created archive {} from {}'.format(tarname, filename))
        else:
            logger.info('Skipped unchanged archive {} from {}'\
                        .format(tarname, filename))


def _newest_modification_time(filepath):
    '''Returns the newest modification time among all files within a folder.'''
    import os
    mtime = os.lstat(filepath).st_mtime
    for root, dirs, entries in os.walk(filepath):
        for entry in dirs + entries:
            entrypath = os.path.join(root, entry)
            mtime = max(mtime, os.lstat(entrypath).st_mtime)
    return mtime


def _get_tar_filter(ignorefile):
    '''Filters out all files matching any pattern in an ignore file'''
    import os
    import fnmatch

    if not os.path.exists(ignorefile):
        return None
    
    with open(ignorefile, 'r') as fobj:
        patterns = [line.strip() for line in fobj.readlines()]

    def tar_filter(tarinfo):
        for pattern in patterns:
            if fnmatch.fnmatch(tarinfo.name, pattern):
                return None
        return tarinfo

    return tar_filter
