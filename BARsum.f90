!********************************************************************************************************************************!      
!                                                                                                                                !
!                 Purpose: Calculate free energy of activation and reaction over the entire reaction coordinate from		 !
!		           perturbation free energy differences in each window calculated by BAR, OPTimized IS, TS and Final
!                          states	                                                                                         !
!                  Author: Mehdi Zare, Department of Chemical Engineering,USC                                                    !
!                    Date: 02.12.17                                                                                              !
!            Modification: This code is the modification of Saleheen's BARAggregate Code                                                                                                      !
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
  character(sl)			:: IfileIS  = 'BARList-GL-IS'
  character(sl)                 :: IfileTS  = 'BARList-GL-TS'
  character(sl)                 :: IfileFS  = 'BARList-GL-FS'
  character(sl)                 :: IfileFEP = 'BARList-FEP'
  character(sl)			:: Ofile    = 'FEP-OPT'

  ! Declaration of variable
  integer			:: elementsIS,elementsTS,elementsFS,elementsFEP, TSIndexOPT, TSIndex=68
  integer                       :: GbarIndex, GsumIndex, i, ii, jj, kk
  real(dp), allocatable         :: GbarIS(:), GbarTS(:), GbarFS(:), GbarFEP(:)
  real(dp), allocatable         :: Gsum(:), Gbar(:)
  real(dp)			:: delGRxn, delGFor
  integer			:: nline, alloc_err
  
  call countline(IfileIS, nline)
  elementsIS  = nline
  write(*,*) "elementsIS=", elementsIS          !!FLAG
  allocate (GbarIS(elementsIS))
  call readarray (IfileIS, elementsIS, GbarIS)

  call countline(IfileFEP, nline)
  elementsFEP  = nline
   write(*,*) "elementsFEP=", elementsFEP          !!FLAG
  allocate (GbarFEP(elementsFEP))
  call readarray (IfileFEP, elementsFEP, GbarFEP)

  call countline(IfileTS, nline)
  elementsTS  = nline
 write(*,*) "elementsTS=", elementsTS          !!FLAG
  allocate (GbarTS(elementsTS))
  call readarray (IfileTS, elementsTS, GbarTS)

  call countline(IfileFS, nline)
  elementsFS  = nline
 write(*,*) "elementsFS=", elementsFS          !!FLAG
  allocate (GbarFS(elementsFS))
  call readarray (IfileFS, elementsFS, GbarFS)
  
  GbarIndex = elementsIS + elementsFEP + elementsTS * 2 + elementsFS
  GsumIndex = GbarIndex + 1
   write(*,*) "GbarIndex=", GbarIndex          !!FLAG
   write(*,*) "GsumIndex=", GsumIndex          !!FLAG
  allocate (Gsum(GsumIndex))
  allocate (Gbar(GbarIndex))
    
    ii = 1
    kk = elementsIS
    jj = kk
    do i = ii, kk         
      Gbar(i) = -GbarIS(jj)
      jj = jj - 1
    end do
   
    ii = kk + 1
    kk = ii + (TSIndex -1) - 1
    jj = 1
    do i = ii, kk
       Gbar(i) = GbarFEP(jj)
       jj = jj + 1
    end do
  
    ii = kk + 1
    kk = ii + elementsTS - 1
    jj = 1
    do i = ii, kk
       Gbar(i) = GbarTS(jj)
       jj = jj + 1
    end do
    
    ii = kk + 1
    kk = ii + elementsTS - 1
    jj = elementsTS
    do i = ii, kk
       Gbar(i) = -GbarTS(jj)
       jj = jj - 1
    end do
    
    ii = kk + 1
    kk = ii + (elementsFEP-TSIndex)
    jj = TSIndex
    do i = ii, kk
       Gbar(i) = GbarFEP(jj)
       jj = jj + 1
    end do
 
   ii = kk + 1
   kk = ii + elementsFS - 1
   jj = 1
   do i = ii, kk
       Gbar(i) = GbarFS(jj)
       jj = jj + 1
    end do

  call cumulativeGsum (GbarIndex, Gbar, GsumIndex, Gsum)
     
  write(*,*) "Gsum(145)= " , Gsum(145)
  write(*,*) "Gsum(73)= " , Gsum(73)
  write(*,*) "Gsum(77)= " , Gsum(77)
  write(*,*) "Gsum(81)= " , Gsum(81)    
  
  TSIndexOPT = elementsIS +  elementsTS + TSIndex
  write(*,*) "delESCF= " , Gsum(GsumIndex)
  write(*,*) "delEFORSCF= " , Gsum(TSIndexOPT)
   
  call getdelG (GsumIndex, TSIndexOPT, Gsum, delGRxn, delGFor)
  call writedelG (Ofile, GsumIndex, Gsum, delGRxn, delGFor)

  deallocate (Gbar, stat = alloc_err)
  deallocate (Gsum, stat = alloc_err)
  deallocate (GbarIS, stat = alloc_err)
  deallocate (GbarTS, stat = alloc_err)
  deallocate (GbarFS, stat = alloc_err)
  deallocate (GbarFEP, stat = alloc_err)

end program BARsum
