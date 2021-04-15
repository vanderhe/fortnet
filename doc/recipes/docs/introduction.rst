.. _sec-introduction:

************
Introduction
************

This documentation is a symbiosis of examples demonstrating the usage of the
Behler-Parrinello-Neural-Network :cite:`bpnn` implementation `Fortnet`,
together with more detailed descriptions of the features provided.

Before You Start
================

The examples assume that you have the `latest stable version
<https://github.com/vanderhe/fortnet/releases/>`_ of Fortnet installed.
Although, many of the recipes may also work with older versions of the code.
Please consult the
`INSTALL.rst <https://github.com/vanderhe/fortnet/blob/master/INSTALL.rst>`_
file for detailed building instructions and troubleshooting.

The recipes in this document often only show the relevant parts of the input. In
order to obtain the full input files and in order to run the examples yourself,
please download the archive, containing all the inputs of the individual
recipes.

.. only :: builder_html or readthedocs

   Download :download:`archive with all inputs
   <_downloads/archives/recipes.tar.bz2>`.

.. only :: not (builder_html or readthedocs)

   This can be downloaded from the `online version of the Fortnet recipes
   <https://fortnet.readthedocs.io/>`_ at https://fortnet.readthedocs.io/.
   
In each recipe you will get an indication where to find the corresponding
directories in the archive with square brackets after the section title (e.g.
[Input: `recipes/basics/firsttrain/`]).

Where to start
==============

The individual chapters are more or less independent from each other, so you may
go directly to the one relevant to your interests. However, if you are new to
Fortnet, please make sure to work through the relevant introductory examples in
the :ref:`sec-basics` chapters first.

Please note that the example outputs in the recipes may have been created with
older versions of Fortnet and therefore could differ slightly in format from
output of the most recent code. The corresponding inputs in the archive should
work, without any changes, with the last stable release of Fortnet.
