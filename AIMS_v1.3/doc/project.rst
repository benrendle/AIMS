Project Summary
===============

Description
-----------

**Name**: "Asteroseismic Inference on a Massive Scale" (*AIMS*)

**Goals**:
  * estimate stellar parameters and credible intervals/error bars
  * chose a representative set or sample of reference models
  * be computationally efficient

**Inputs**:
  * classic constraints and error bars (Teff, L, ...)
  * seismic constraints and error bars (individual frequencies)

**Requirements**:
  * a *precalculated* grid of models including:

    - the models themselves
    - parameters for the model (M, R, Teff, age, ...)
    - theoretical frequency spectra for the models

**Methodology**:
  * applies an MCMC algorithm based on the python package `emcee <http://dan.iel.fm/emcee/current/>`_.
    Relevant articles include:

    - `Bazot et al. (2012, MNRAS 427, 1847) <http://ukads.nottingham.ac.uk/abs/2012MNRAS.427.1847B>`_
    - `Gruberbauer et al. (2012, ApJ 749, 109) <http://ukads.nottingham.ac.uk/abs/2012ApJ...749..109G>`_

  * interpolates within the grid of models using Delaunay tessellation
    (from the `scipy.spatial <http://docs.scipy.org/doc/scipy/reference/spatial.html>`_
    package which is based on the `Qhull <http://www.qhull.org/>`_ library)
  * modular approach: facilitates including contributions from different
    people

Contributors
------------

**Author**:

  * Daniel R. Reese

**Comments, corrections, suggestions, and contributions**:

  * Diego Bossini
  * Gael Buldgen
  * Tiago L. Campante
  * William J. Chaplin
  * Hugo R. Coelho
  * Guy R. Davies
  * Benoît D. C. P. Herbert
  * James S. Kuszlewicz
  * Yveline Lebreton
  * Martin W. Long
  * Mikkel N. Lund
  * Andrea Miglio
  * Ben Rendle

Supplementary material
----------------------

  * a more technical :download:`overview <./files/Overview.pdf>` of AIMS
  * a PDF version of this documentation may be downloaded
    :download:`here <./files/AIMS.pdf>`

Copyright information
---------------------

  * the AIMS project is distributed under the terms of the
    `GNU General Public License, version 3 <http://www.gnu.org/licenses/gpl-3.0.en.html>`_
  * a copy of of this license may be downloaded :download:`here <../COPYING>`
    and should also be included in ``AIMS.tgz``
