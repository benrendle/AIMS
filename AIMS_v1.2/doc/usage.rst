Usage
=====

There are three different ways of using AIMS:

  1. generating a binary file with the grid of models (including names,
     global parameters, and pulsation frequencies).

     .. note::
        This step must be carried out before the following two steps
        as these require the above binary file to function correctly.

  2. carrying out tests to evaluate the accuracy of the interpolation
     for a given grid of models.
  3. finding the properties of an observed star thanks to its classic
     and seismic parameters.

The way AIMS is used is decided by the values given in the ``AIMS_configure.py``
file, which also contains a number of other control parameters.  Extensive
comments are included in this file to help the user know how to set the various
parameters.

Generating a binary grid
------------------------

  Requirements:
    * a grid of models, including the pulsation frequencies; the format for
      the files with the pulsation frequencies is described in
      :py:meth:`model.Model.read_file`.
    * a list with the paths and a set of global parameters for each model in
      the grid; the format this file is described in
      :py:meth:`model.Model_grid.read_model_list`.

  Relevant parameters in ``AIMS_configure.py``:
    * ``write_data``: set this to ``True`` so that AIMS will write binary grid.
    * ``list_grid``: set this to the filename of the file with the list of
      paths and global parameters.
    * ``binary_grid``: set this to the filename of the file which will contain
      the binary data.
    * ``grid_params``: specify the parameters relevant to the grid (excluding
      age, which is dealt with separately).  Different options can be found
      in the source to :py:func:`model.Model.string_to_param`.
    * ``npositive``: set this to ``True`` to only save modes with :math:`n \ge 0`
      in the binary file.
    
  To run AIMS in this configuration, just type the following in a terminal
  window::

  ./AIMS.py

Testing the accuracy of the interpolation
-----------------------------------------

  Requirements:
    * a binary grid of models as produced by AIMS

  Relevant parameters in ``AIMS_configure.py``:
    * ``write_data``: set this to ``False`` otherwise a binary grid will
      be produced, the interpolations tests will not be carried out.
    * ``test_interpolation``: set this to ``True`` so that AIMS
      will carry out the interpolation tests.
    * ``interpolation_file``: specify the name of the file in which to
      write the results from the interpolation test in binary format.
      These results can be plotted using ``plot_interpolation_test.py``.

  To run AIMS in this configuration, just type the following in a terminal
  window::

  ./AIMS.py

Characterising an observed star
-------------------------------

  Requirements:
    * a binary grid of models as produced by AIMS
    * a file with the observational data; the format for this file is
      similar to the format used for the
      `Asteroseismic Modeling Portal (AMP) <https://amp.phys.au.dk/>`_ with
      some simplifications and is described below.  It will be read by
      :py:meth:`AIMS.Likelihood.read_constraints`

  Relevant parameters in ``AIMS_configure.py``:
    * ``write_data``: set this to ``False``
    * ``test_interpolation``: set this to ``False``
    * most of the parameters in this file - see comments for details

  To run AIMS in this configuration, just type the following in a terminal
  window::

  ./AIMS.py file_with_constraints

  where ``file_with_constraints`` is the file with the observational
  constraints.
