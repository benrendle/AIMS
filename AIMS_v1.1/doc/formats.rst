File formats
============

Format of a file with a list of models and properties:
------------------------------------------------------

  Description:
    * The first line is a header.   It contains the root folder
      (including the final slash) with the grid of models and
      optionally, a suffix for the names of the files with the
      theoretical pulsation frequencies.  For example::

        /home/dreese/models_inversions/Grid_mesa_MS/  .freq

    * Each of the following lines correspond to one model in
      the grid.  They are composed of 8 or more columns with the
      following information:

      1. The second part of the path for the given model.
         When concatenated with the prefix on the first line,
         this should give the full path to the model.  If,
         furthermore, the suffix from the first line is
         appended to it, it gives the name of the file
         with the frequencies.
      2. The stellar mass in :math:`\mathrm{g}`
      3. The stellar radius in :math:`\mathrm{cm}`
      4. The stellar luminosity in :math:`\mathrm{g.cm^2.s^{-3}}`
      5. The metallicity
      6. The hydrogen content
      7. The stellar age in :math:`\mathrm{Myrs}`
      8. The effective temperature in :math:`\mathrm{K}`
      9. (user-defined) This and the following columns correspond
         to the parameters specified in the ``user_params``
         variable given in ``AIMS_configure.py``.

    * Except for the first line, the order of the lines does
      not matter.  AIMS will construct evolutionary tracks based
      on the parameters selected in the ``grid_params`` variable
      given in ``AIMS_configure.py``, and sort them according to
      age.

  Example:
    Here's an example of a file read by AIMS (via the
    :py:meth:`model.Model_grid.read_model_list` method)::

      /home/dreese/models_inversions/Grid_mesa_MS/  .freq
      M0.80/LOGS_M0.80/M0.80Z0.0028Y0.2536/m0.80Y0.2536Z0.0028a1.8ovh0.2ovhe0_n1.profile.FGONG    1.59136E+33   5.02248266E+10   2.3309799E+33  0.0028  0.7436   1.0000000E-04   6000.94326   7.432106E-01
      M0.80/LOGS_M0.80/M0.80Z0.0028Y0.2536/m0.80Y0.2536Z0.0028a1.8ovh0.2ovhe0_n10.profile.FGONG   1.59136E+33   5.02990358E+10   2.1724140E+33  0.0028  0.7436   2.0974688E+03   5891.82623   6.146083E-01
      M0.80/LOGS_M0.80/M0.80Z0.0028Y0.2536/m0.80Y0.2536Z0.0028a1.8ovh0.2ovhe0_n11.profile.FGONG   1.59136E+33   5.04940406E+10   2.2011824E+33  0.0028  0.7436   2.3237113E+03   5899.81040   6.001537E-01

    It contains three models.  The structure of the first model can
    be found in the following file::

      /home/dreese/models_inversions/Grid_mesa_MS/M0.80/LOGS_M0.80/M0.80Z0.0028Y0.2536/m0.80Y0.2536Z0.0028a1.8ovh0.2ovhe0_n1.profile.FGONG

    and its frequencies in this file::
  
      /home/dreese/models_inversions/Grid_mesa_MS/M0.80/LOGS_M0.80/M0.80Z0.0028Y0.2536/m0.80Y0.2536Z0.0028a1.8ovh0.2ovhe0_n1.profile.FGONG.freq

    The ninth column corresponds to the central hydrogen
    content, as specified by the contents of the ``user_params``
    variable from ``AIMS_configure.py``::

      user_params = (("Xc", r'Central hydrogen, $%sX_c%s$'),)

Format of a file with theoretical frequencies:
----------------------------------------------

  Description:
    * the first line is a header (and is skipped)
    * the following lines contain five columns which correspond
      to l, n, frequency, a_value, inertia

      - the a_value column is ignored, so it could contain anything.
        ``InversionKit`` will typically put the difference between
        the numerical and variational frequencies in that column.

  Example:
    Here's an example of a file with theoretical pulsation
    frequencies which can be read by AIMS (via the
    :py:meth:`model.Model.read_file` method)::

          #l  n         nu_theo (muHz)  nu_var-nu_theo (muHz)                Inertia
          0  15  3.225852209451052e+03  1.312960435370769e-03  3.233628965187502e-09
          0  16  3.421699035498995e+03 -2.482639610207116e-03  2.229252226305757e-09
          0  17  3.615805033992529e+03  3.993051574070705e-03  1.618154348529283e-09
          0  18  3.809740380503104e+03  9.650666734160040e-04  1.250359548964621e-09
          0  19  4.003716857281849e+03 -7.991676880010345e-03  1.033914933206195e-09
          0  20  4.198691419457581e+03  1.742711681799847e-03  8.866985261874711e-10
          1  15  3.316007619955153e+03  5.056100344972947e-03  2.715966891128009e-09
          1  16  3.511258977705781e+03  1.855844971032639e-04  1.902147334986236e-09
          1  17  3.705576731149742e+03 -2.505276897409203e-03  1.424266453221534e-09
          1  18  3.899485457373566e+03  5.212276555539575e-03  1.134594720287415e-09
          1  19  4.094401244305849e+03  6.020260397235688e-03  9.579611596023003e-10
          1  20  4.289716814475406e+03 -1.019475706561934e-02  8.344804874142957e-10
          2  15  3.399280335063532e+03 -8.466318249702454e-04  2.315947651745295e-09
          2  16  3.594141943503532e+03  4.712417365681176e-03  1.665322627996223e-09
          2  17  3.788792185755381e+03 -1.167229517704982e-03  1.277569745555387e-09
          2  18  3.983271067684743e+03 -6.187409578615188e-03  1.048757367028520e-09
          2  19  4.178866833517976e+03  6.893199766636826e-03  8.963691946280509e-10
          2  20  4.374959711016754e+03  3.274638356742798e-03  7.911508926344487e-10
          3  15  3.476224140192640e+03 -2.524210208321165e-03  2.009476926536794e-09
          3  16  3.671438520072859e+03  2.351724720028869e-04  1.485336526791650e-09
          3  17  3.866350877376991e+03  5.643782460992952e-03  1.167619144668003e-09
          3  18  4.061929209725198e+03 -1.552865011490212e-03  9.789648655155361e-10
          3  19  4.258077196700047e+03 -8.629839649984206e-03  8.472972126693386e-10
          3  20  4.455063887754256e+03  1.484804296796938e-02  7.528069568152023e-10

Format of a file with observational constraints:
------------------------------------------------

  Description:
    * a collection of lines with frequency data with either
      (l, freq, error_bar) or (l, n, freq, error_bar) (depending on
      the value of ``read_n`` in the ``AIMS_configure.py``
      file).  For example::

        0 1503.5 0.16

      or the following if specifying the radial order::

        0 15 1503.5 0.16

    * a collection of lines with classical constraints.  These
      start with the name of the relevant parameter (see
      possible options in :py:func:`model.Model.string_to_param`)
      followed by a description of its probability distribution
      function.  This probability distribution function is
      specified in two possible ways:

      - it is implicitly assumed to be Gaussian.  In this situation
        it is only necessary to specify the mean value and the
        one sigma error bar.  For example::

          Teff 6100 80

      - it is explicitly specified (different options are given
        in :py:class:`AIMS.Distribution`)::

          Teff Uniform 6000 6200

    * anything following a ``#`` is a comment

    * the order of the lines does not matter

  Examples:
    * example of a file where n is *not* specified::

          0 1582.20 0.13  # this is a (useless) comment
          0 1684.02 0.16
          0 1785.57 0.15
          1 1526.55 0.29
          1 1628.90 0.30
          1 1730.45 0.17
          2 1575.49 0.82
          2 1676.25 0.51
          2 1777.62 0.27
          Teff 6060.00 84.00
          Fe_H -0.20 0.09

    * example of a file where n is specified::

          0 15 1582.20 0.13
          0 16 1684.02 0.16
          Teff 6060.00 84.00 # AIMS doesn't worry about the order of the lines
          0 17 1785.57 0.15
          1 14 1526.55 0.29
          1 15 1628.90 0.30 
          1 16 1730.45 0.17
          2 14 1575.49 0.82
          2 15 1676.25 0.51
          2 16 1777.62 0.27
          Fe_H -0.20 0.09

  Differences with `AMP <https://amp.phys.au.dk/>`_:
    * the number of frequencies does not need to be specified
      (if this line contains supplementary parameters, than
      ``AIMS.py`` may confuse it with frequency data)
    * there are no flags (one should adjust the parameters in
      ``AIMS_configure.py`` instead)
    * the order of the lines is not important (one can mix
      the classic and seismic observables)
    * it is possible to specify radial orders (depending on
      the value of ``read_n`` in the ``AIMS_configure.py``
      file)
    * the treatment of non-seismic constraints is more flexible

      - a larger variety of non-seismic constraints can be included
        (see possible options in :py:func:`model.Model.string_to_param`)
      - full parameter names are allowed (and preferred); for compatibility
        with `AMP <https://amp.phys.au.dk/>`_, the same one letter
        abbreviations are also allowed
      - it is possible to specify the probability distribution function
