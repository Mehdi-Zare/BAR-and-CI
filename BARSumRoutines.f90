module BARsumroutines
use constants
implicit none

   ! Selective parameters for usage in this specific script
   real(dp), parameter :: Temp = 500.00_dp
   real(dp), parameter :: kBT  = kB * Temp
   real(dp), parameter :: beta = one / kBT

   ! Gas phase partition functions and free eenrgies for COH cleavage of Ethylene Glycol calculated at 500K
   real(dp), parameter :: ZPCRxn        =  -0.06993552395_dp
   real(dp), parameter :: ZPCFor        =  -0.12930749364_dp
   real(dp), parameter :: QvibR         =   34351.480629083_dp
   real(dp), parameter :: QvibTS        =   21147.2573379462_dp
   real(dp), parameter :: QvibP         =   10818.22645197180_dp  
 
contains

!****************************************************************************************
!	                   Procedure for counting lines of a file 							!
!****************************************************************************************
    
    subroutine countline (filename, nline)			        				
    implicit none

	character(sl), intent(in) 	:: filename
	integer,       intent(out)	:: nline
	integer			        :: ierror
	
	nline=0
     	OPEN(UNIT = 10, file = filename, status = 'old',&
             action = 'read', iostat = ierror)
             do
           	read(10, *, iostat = ierror)
           	if (ierror /= 0) exit
		nline = nline + 1
             end do
     	CLOSE(UNIT = 10)
    
    end subroutine

!****************************************************************************************
!			 Procedure for reading energy values reported in input files    			!
!****************************************************************************************
        
    subroutine readarray (filename, Narr, wf)
    implicit none
    character(sl), intent(in) 	:: filename
    integer,       intent(in) 	:: Narr	
    real(dp),      intent(out) 	:: wf(Narr)
    real(dp)                    :: temp2
    integer			:: temp, i, ierror
        
        OPEN(Unit = 10, file = filename, status = 'old', &
             action = 'read', iostat = ierror)
             do i = 1, Narr
            	read(10, *, iostat = ierror) temp, temp2, wf(i)
             if  (ierror /= 0) exit
	     end do
    	CLOSE(10)  
    
    end subroutine

!****************************************************************************************
!			 Procedure for calulating forward Gsum				!
!****************************************************************************************
    subroutine cumulativeGsum(Narr, Gbar, GsumIndex, Gsum)
    implicit none
    integer,  intent(in)  :: Narr
    integer,  intent(in)  :: GsumIndex
    real(dp), intent(in)  :: Gbar(Narr)
    real(dp), intent(out) :: Gsum(GsumIndex)
    integer	          :: i
    Gsum = zero
    do i = 1, Narr
      Gsum(i+1) = Gsum(i) + Ha2eV * Gbar(i) 
    WRITE(*,*) Gsum(i+1)                                        !!!!!FLAG
    end do
    end subroutine

!****************************************************************************************
!	Procedure for calulating delG using gas phase partition functions		!
!****************************************************************************************
    subroutine getdelG(GsumIndex, TSIndex, Gsum, delGRxn, delGFor)
    implicit none
    integer,  intent(in)  :: GsumIndex
    integer,  intent(in)  :: TSIndex
    real(dp), intent(in)  :: Gsum(GsumIndex)
    real(dp), intent(out) :: delGRxn
    real(dp), intent(out) :: delGFor
    real(dp)		  :: delE, delEFor
    real(dp)		  :: QelecRxn, QelecFor
    
    delE     = ZPCRxn + Gsum(GsumIndex)
    delEFor  = ZPCFor + Gsum(TSIndex)
    QelecRxn = exp ( - delE / kBT )
    QelecFor = exp ( - delEFor / kBT )
    delGRxn  = - kBT * log ((QelecRxn * QvibP) / QvibR )
    delGFor  = - kBT * log ((QelecFor * QvibTS) / QvibR )

    end subroutine

!************************************************************************************************!  
!   Procedure for writing cumulative free energies and free energies of activation and reaction  !
!************************************************************************************************!
    subroutine writedelG(filename, GsumIndex, Gsum, delGRxn, delGFor)
    implicit none
    character(sl), intent(in)  	:: filename
    integer, intent(in)   	:: GsumIndex
    real(dp), intent(in)  	:: Gsum(GsumIndex)
    real(dp), intent(in) 	:: delGRxn
    real(dp), intent(in) 	:: delGFor
    integer			:: i, ierror

	OPEN (UNIT = 20, file = filename, access= 'sequential', &
    	  status = 'replace', action = 'write', position = 'append', &
          iostat = ierror)
          write(20,1020)
          1020 format (10X, ' deltaGRxn',15X, 'deltaGactivation')
          write(20,1030)
          1030 format (6X, '==================',7X, '=========================') 
          write(20,1040) delGRxn, delGFor
          1040 format (3X, F20.13, 8X, F20.13)
          write(20,1050)
          1050 format (2X,'Index', 3X,'Cumulative Free Energies')
          write(20,1060)
          1060 format (2X,'=====',3X, '========================')          
          do i = 1, GsumIndex
            write(20, 1070) i, Gsum(i)
            1070 format (4X, I4, 8X, F20.13)
          end do  
        CLOSE(10)      
    
    end subroutine    

end module BARsumroutines
