!------------------------------------------------------------------------------
! This method takes two sets of modes and produces a third set by combining
! the two first sets using the interpolation weights coef1 and coef2. 
! Only modes which are common to both sets will be present in the output set.
!------------------------------------------------------------------------------
! INPUT:
! coef1 = first interpolation coefficient
! narr1 = n values of the first set of modes
! larr1 = l values of the first set of modes
! farr1 = frequencies of the first set of modes
! iarr1 = inertias of the first set of modes
! n1    = number of modes in the first set of modes
! coef2 = second interpolation coefficient
! narr2 = n values of the second set of modes
! larr2 = l values of the second set of modes
! farr2 = frequencies of the second set of modes
! iarr2 = inertias of the second set of modes
! n2    = number of modes in the second set of modes
!
! OUTPUT: 
! narr3 = n values of the resultant set of modes
! larr3 = l values of the resultant set of modes
! farr3 = frequencies of the resultant set of modes
! iarr3 = inertias of the resultant set of modes
! n3    = number of modes in the resultant set of modes
!------------------------------------------------------------------------------
      subroutine combine_modes(coef1,narr1,larr1,farr1,iarr1, n1, &
                               coef2,narr2,larr2,farr2,iarr2, n2, &
                               narr3,larr3,farr3,iarr3, n3)

      implicit none

      ! input arguments
      integer, intent(in) :: n1, n2
      integer(kind=2), intent(in) :: narr1(0:n1-1), narr2(0:n2-1)  ! np.int16
      integer(kind=1), intent(in) :: larr1(0:n1-1), larr2(0:n2-1)  ! np.int8
      real(kind=8), intent(in) :: coef1, coef2
      real(kind=8), intent(in) :: farr1(0:n1-1),iarr1(0:n1-1)
      real(kind=8), intent(in) :: farr2(0:n2-1),iarr2(0:n2-1)

      ! output arguments
      integer, intent(inout) :: n3
      !f2py intent(in,out) n3
      integer(kind=2), intent(inout) :: narr3(0:n3-1)
      !f2py intent(in,out) narr3
      integer(kind=1), intent(inout) :: larr3(0:n3-1)
      !f2py intent(in,out) larr3
      real(kind=8), intent(inout) :: farr3(0:n3-1),iarr3(0:n3-1)
      !f2py intent(in,out) farr3, iarr3

      integer i1, i2, i3

  
      i1 = 0
      i2 = 0
      i3 =-1 
      do while ((i1.lt.n1).and.(i2.lt.n2))
        if (larr1(i1).lt.larr2(i2)) then
          i1 = i1+1
          cycle
        endif
        if (larr1(i1).gt.larr2(i2)) then
          i2 = i2+1
          cycle
        endif
        if (narr1(i1).lt.narr2(i2)) then
          i1 = i1+1
          cycle
        endif
        if (narr1(i1).gt.narr2(i2)) then
          i2 = i2+1
          cycle
        endif
        ! now the two modes have the same n and l values:

        ! first increment i3
        i3 = i3 + 1
        !! sanity check (to be removed eventually)
        !if (i3.gt.n3) then
        !  print*,i3,n3
        !  stop "array3 too small in aims_fortran.combine"
        !endif

        narr3(i3) = narr1(i1)
        larr3(i3) = larr1(i1)
        farr3(i3) = coef1*farr1(i1) + coef2*farr2(i2) 
        iarr3(i3) = coef1*iarr1(i1) + coef2*iarr2(i2) 

        i1 = i1+1
        i2 = i2+1
      enddo
      n3 = i3+1

      end subroutine combine_modes

!------------------------------------------------------------------------------
! This method finds a mapping between an observed set of modes and a
! theoretical set of modes from a model.  The mapping is obtained from
! matching (l,n) values.
!------------------------------------------------------------------------------
! INPUT:
! nobs     = n values of the observed modes
! lobs     = l values of the observed modes
! size_obs = number of observed modes
! nmod     = n values of the theoretical modes
! lmod     = l values of the theoretical modes
! size_mod = number of theoretical modes
!
! OUTPUT: 
! mode_map = mapping between observed and theoretical modes
! nmissing = number of unmatched observed modes
!------------------------------------------------------------------------------
      subroutine find_map_n(nobs, lobs, size_obs, nmod, lmod, &
                            size_mod, mode_map, nmissing)

      implicit none

      ! input arguments
      integer, intent(in) :: size_obs, size_mod
      integer(kind=2), intent(in) :: nmod(0:size_mod-1), nobs(0:size_obs-1)
      integer(kind=1), intent(in) :: lmod(0:size_mod-1), lobs(0:size_obs-1)

      ! output arguments
      integer, intent(out) :: nmissing
      integer, intent(inout) :: mode_map(0:size_obs-1)
      !f2py intent(in,out) mode_map

      integer iobs, imod

      !initialisation
      iobs = 0
      imod = 0
      nmissing = 0
      mode_map = -1

      ! NOTE: this assumes the observed modes are sorted according to (l,n)
      do while((iobs.lt.size_obs).and.(imod.lt.size_mod))

        if (lmod(imod).lt.lobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        if (lmod(imod).gt.lobs(iobs)) then
           iobs = iobs + 1
           nmissing = nmissing + 1
           cycle
        endif

        if (nmod(imod).lt.nobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        if (nmod(imod).gt.nobs(iobs)) then
           iobs = iobs + 1
           nmissing = nmissing + 1
           cycle
        endif

        mode_map(iobs) = imod

        imod = imod + 1
        iobs = iobs + 1
      enddo

      ! keep track of the number of missing modes:
      nmissing = nmissing + (size_obs-iobs)
      end subroutine find_map_n
        
!------------------------------------------------------------------------------
! This method finds a mapping between an observed set of modes and a
! theoretical set of modes from a model.  The mapping is obtained from
! l values and frequency proximity.
!------------------------------------------------------------------------------
! INPUT:
! fobs     = frequencies of the observed modes (units = muHz)
! lobs     = l values of the observed modes
! size_obs = number of observed modes
! fmod     = frequencies of the theoretical modes (units = muHz)
! lmod     = l values of the theoretical modes
! ind      = index array for which the theoretical modes are sorted according
!            to (l,freq)
! size_mod = number of theoretical modes
!
! OUTPUT: 
! mode_map = mapping between observed and theoretical modes
! nmissing = number of unmatched observed modes
!------------------------------------------------------------------------------
      subroutine find_map_freq(fobs, lobs, size_obs, fmod, lmod, ind, &
                            size_mod, mode_map, nmissing) 

      implicit none

      ! input arguments
      integer, intent(in) :: size_obs, size_mod
      integer, intent(in) :: ind(0:size_mod-1)
      integer(kind=1), intent(in) :: lmod(0:size_mod-1), lobs(0:size_obs-1)
      real(kind=8), intent(in) ::    fmod(0:size_mod-1), fobs(0:size_obs-1)

      ! output arguments
      integer, intent(out) :: nmissing
      integer, intent(inout) :: mode_map(0:size_obs-1)
      !f2py intent(in,out) mode_map

      real(kind=8) :: diff0, diff1, diffmin
      real(kind=8), allocatable :: diffs(:)
      integer :: iobs, imod, imod0, imod1, val, j, jmin

      !initialisation
      iobs = 0
      imod = 0
      nmissing = 0
      mode_map = -1
      allocate(diffs(0:size_obs-1))

      ! This is a two step process:
      !   1. The first step finds the nearest theoretical mode to
      !      each observed mode.
      !   2. When a theoretical mode has several observed modes that
      !      correspond to it, only the nearest observed mode is kept.
      !      The other observed modes are no longer mapped to the
      !      theoretical mode and become "missing".

      ! Step 1: find nearest theoretical modes:
      do while((iobs.lt.size_obs).and.(imod.lt.size_mod))
        imod1 = ind(imod)

        if (lmod(imod1).lt.lobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        if (lmod(imod1).gt.lobs(iobs)) then
           iobs = iobs + 1
           nmissing = nmissing + 1
           cycle
        endif
        
        if (fmod(imod1).lt.fobs(iobs)) then
           imod = imod + 1
           cycle
        endif

        diff1 = fmod(imod1) - fobs(iobs)
        if (imod.gt.0) then
          imod0 = ind(imod-1)
          if (lmod(imod0).eq.lobs(iobs)) then
            diff0 = fobs(iobs) - fmod(imod0)
            if (diff0 < diff1) then
              mode_map(iobs) = imod0
              diffs(iobs) = diff0
            else
              mode_map(iobs) = imod1
              diffs(iobs) = diff1
            endif
          else
            mode_map(iobs) = imod1
            diffs(iobs) = diff1
          endif
        else
          mode_map(iobs) = imod1
          diffs(iobs) = diff1
        endif

        iobs = iobs + 1
      enddo

      ! keep track of the number of missing modes:
      nmissing = nmissing + (size_obs-iobs)

      ! Step 2: filter out cases where more than one observed mode 
      !         corresponds to the same theoretical mode.
      do iobs = 0, size_obs-1
        if (mode_map(iobs).eq.-1) cycle

        val = mode_map(iobs)
        jmin = iobs
        diffmin = diffs(iobs)
        do j=iobs+1,size_obs-1
          if (mode_map(j).ne.val) exit
          if (diffs(j).lt.diffmin) then
            jmin = j
            diffmin = diffs(j)
          endif
        enddo
        mode_map(iobs:j-1) = -1
        mode_map(jmin) = val
        nmissing = nmissing + j-iobs-1
      enddo

      deallocate(diffs)

      end subroutine find_map_freq

!------------------------------------------------------------------------------
! This method calculates differences between observed and theoretical frequency
! combinations.
!------------------------------------------------------------------------------
! INPUT:
! freq           = theoretical frequencies
! mode_map       = mapping between observed and theoretical frequencies
! values         = values of observed frequency combinations
! ncoeff         = number of terms in the frequency combinations
! coeff          = coefficients in the frequency combinations
! indices        = mode indices in the frequency combinations
! nmodes         = number of theoretical modes
! nobs           = number of observed modes
! ncomb          = number of frequency combinations
! nmax           = maximum number of terms in a frequency combination
!
! OUTPUT:
! dvalues        = differences between observed and theoretical frequency
!                  combinations
!------------------------------------------------------------------------------
      subroutine compare_frequency_combinations(freq,mode_map,values,&
                    ncoeff,coeff,indices,nmodes,nobs,ncomb,nmax,dvalues)

      implicit none

      integer, intent(in) :: nmodes, nobs, ncomb, nmax
      real(kind=8), intent(in) :: freq(0:nmodes-1),values(0:ncomb-1),&
                                  coeff(0:nmax-1,0:1,0:ncomb-1)
      integer, intent(in) :: mode_map(0:nobs-1), ncoeff(0:1,0:ncomb-1),&
                             indices(0:nmax-1,0:1,0:ncomb-1)
      real(kind=8), intent(inout) :: dvalues(0:ncomb-1)
      !f2py intent(in,out) dvalues

      real(kind=8) :: den
      integer i, j

      ! NOTE: dvalues has already been initialised to 0 in AIMS.py
      do i=0, ncomb-1
        do j=0,ncoeff(0,i)-1
          dvalues(i) = dvalues(i) + coeff(j,0,i)*freq(mode_map(indices(j,0,i)))
        enddo
        if (ncoeff(1,i).gt.0) then
          den = 0d0
          do j=0,ncoeff(1,i)-1
            den = den + coeff(j,1,i)*freq(mode_map(indices(j,1,i)))
          enddo
          dvalues(i) = dvalues(i)/den
        endif
        dvalues(i) = dvalues(i) - values(i)
      enddo

      end subroutine compare_frequency_combinations

!------------------------------------------------------------------------------
! This method reads a set of modes from a file in "agsm" format, a fortran
! binary format used by ADIPLS.
!------------------------------------------------------------------------------
! INPUT:
! filename       = name of the agsm file
! npositive      = specifies whether to retain only modes with n >= 0
! below_cutoff   = specifies whether to only retain modes below the cutoff frequency
! freqlim        = upper frequency bound. Mode with frequencies above this bound
!                  are discarded
!
! OUTPUT:
! narr           = n values of the modes which are read
! larr           = l values of the modes which are read
! farr           = frequencies of the modes which are read
! iarr           = inertias of the modes which are read
! nn             = number of modes which are read
! exceed_freqlim = specifies whether a mode exceeded the frequency limit
!------------------------------------------------------------------------------
      subroutine read_file_agsm(filename,npositive,below_cutoff, &
                                freqlim,narr,larr,farr,iarr,nn,  &
                                exceed_freqlim)

      implicit none
      integer, parameter :: nmax = 1000
      character(len=*), intent(in) :: filename
      logical, intent(in) :: npositive, below_cutoff
      real(kind=8), intent(in) :: freqlim
      real(kind=8), intent(out) :: farr(nmax), iarr(nmax)
      integer(kind=2), intent(out) :: narr(nmax)
      integer(kind=1), intent(out) :: larr(nmax)
      logical, intent(out) :: exceed_freqlim
      integer, intent(out) :: nn
      real(kind=8) :: cs(50)
      integer :: ics(8), ntot, i, j

      ! transform the last part of the cs array into integers:
      equivalence(ics(1), cs(39))

      ntot = 0
      open(unit=31,file=filename,status='old',form='unformatted')
      ! note: the convert keyword can be used to change the endian
      do
        read(31,end=10) (cs(i),i=1,50)
        ntot = ntot + 1
      enddo
10    close(31)
      
      ! sanity check:
      if (ntot.gt.nmax) then
        stop "Please increase nmax in aims_fortran.read_file_agsm"
      endif

      exceed_freqlim = .False.
      open(unit=31,file=filename,status='old',form='unformatted')
      nn = 1 
      do i=1,ntot
        read(31) (cs(j),j=1,50)
        if (below_cutoff.and.(ics(5).ne.10010)) cycle
        if (npositive.and.(cs(19).lt.-0.1d0)) cycle
        if (cs(37)*1d3.gt.freqlim) then
          exceed_freqlim = .True.
          cycle
        endif
        larr(nn) = nint(cs(18))
        narr(nn) = nint(cs(19))
        farr(nn) = cs(37)*1d3 ! Richardson frequency in muHz
        !farr(nn) = cs(27)*1d3 ! variational frequency in muHz
        iarr(nn) = cs(24)
        nn = nn + 1
      enddo
      close(31)
      nn = nn - 1

      end subroutine read_file_agsm

