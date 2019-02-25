!********************************************************************************************************************************!      
!                                                                                                                                !
!                 Purpose: Calculate free energy of activation and reaction over the entire reaction coordinate from		 !
!		           perturbation free energy differences in each window calculated by BAR	                         !
!                  Author: Mohammad Shamsus Saleheen, Department of Chemical Engineering,USC                                     !
!                    Date: 06.20.17                                                                                              !
!            Modification:                                                                                                       !
! Reasons of Modification:                                                                                                       !
!                   Notes: I do not recommend using this script to anyone else since it's too much customized for a specific     !
!                          purpose.                                                                                              !
!                                                                                                                                !
!********************************************************************************************************************************!

program BARsum
  use constants	
  use BARsumroutines
  implicit none

  ! Input and output files
  character(sl)			:: Ifile = 'BARList'
  character(sl)			:: Ofile = 'DeltaF-BAR'

  ! Declaration of variable
  integer			:: elements, GsumIndex, TSIndex=68
  real(dp), allocatable         :: Gbar(:)
  real(dp), allocatable         :: Gsum(:)
  real(dp)			:: delGRxn, delGFor
  integer			:: nline, alloc_err

  call countline(Ifile, nline)
  elements  = nline
  GsumIndex = nline + 1
  allocate (Gbar(elements))
  allocate (Gsum(GsumIndex))

  call readarray (Ifile, elements, Gbar)
  call cumulativeGsum (elements, Gbar, GsumIndex, Gsum)
  call getdelG (GsumIndex, TSIndex, Gsum, delGRxn, delGFor)
  call writedelG (Ofile, GsumIndex, Gsum, delGRxn, delGFor)

  deallocate (Gbar, stat = alloc_err)
  deallocate (Gsum, stat = alloc_err)
      
end program BARsum
