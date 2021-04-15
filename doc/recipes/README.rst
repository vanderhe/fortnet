***************
Fortnet Recipes
***************

The `docs/` folder contains the source of the Fortnet Recipes book. The content
is currently rendered on `fortnet-recipes.readthedocs.io
<http://fortnet-recipes.readthedocs.io>`_. You can also render it locally by
using the Sphinx documentation system.

This material is licensed under the CC-BY-SA license (see `LICENSE <LICENSE>`_).

Contributions are very welcome!


Technical notes
===============

* The book is written in the restructured text format and is contained in the
  `docs/` folder. In order to process it, you need Sphinx version 1.8 or above.

* Before the text of the book is built, the building system creates a tar.bz2
  archive for each top level file/directory in the `docs/_archives/`
  folder. Created archives will be then copied to the `docs/_downloads/archives`
  folder. You may then reference those .tar.bz2 archives from within the text of
  the book.
