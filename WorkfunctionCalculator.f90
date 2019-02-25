!********************************************************************************************************************************
!							                                                                        !														   			
!		  	   Purpose: Caclulate workfunctions [ U(perturbed)-U(sampled) ] from lof files				!													   			
!	               	    Author: Mohammad Shamsus Saleheen, Department of Chemical Engineering,USC                           !               
!	       	              Date: 06.18.17											!					   
!   	              Modification:												!					 
!          Reasons of Modification:												!					   
!                                                               								!					   
!********************************************************************************************************************************


program workfunction
implicit none


        !*******************************************    Member Elements ******************************************!
        !                                                                                                         !
        !       'infile'-- input file name for configurational energy of each conformations                       !
        !       'uo(:)'--- conformational energies (array) for states over which MM space has been sampled        !
        !       'u1(:)'--- conformational energies (array) for perturbed states                                   !
        !       'nline'--- number of lines in input file                                                          !
        !       'datapoins'--number of conformational energies reported in input file                             !
        !       'upoints'--number of conformational energies for initial and final states                         !
        !       'delta'--- change in conformational energy between initial and final states                       !
        !       'DeltaF'-- output filename                                                                        !
        !       'countline(filename,nline)'-- subroutine for counting lines 	                                  !
        !                                                                                                         !
        !*********************************************************************************************************!


! Setting data type 
integer, parameter:: dp = selected_real_kind(15, 307)


! Input variables and filenames
real(dp), allocatable:: uo(:), u1(:)
character(20)::infile='log.MMS_ENERGY'


! Temporary variables 
integer:: nline							
integer:: datapoints
integer:: upoints
integer:: alloc_err,ierror
integer:: i,j,k,temp


! Output variables
real(dp), allocatable:: delta(:)
character(20)::outfile='DeltaF'

! Count number of lines and data points
call countline(infile,nline)
datapoints=nline-4
upoints=datapoints/2

! Allocate representative arrays
allocate(uo(upoints))
allocate(u1(upoints))
allocate(delta(upoints))

! Open input file to read values

	OPEN(UNIT=10,file=infile,status='old',action='read',&
         iostat=ierror)

         	read(10,*,iostat=ierror)
            read(10,*,iostat=ierror)
         		do i=1,upoints
					read(10,*,iostat=ierror)temp, uo(i)
                    if (ierror /= 0) exit
                end do
			read(10,*,iostat=ierror)
            read(10,*,iostat=ierror)
            	do j=1,upoints
					read(10,*,iostat=ierror)temp, u1(j)
                    if (ierror /= 0) exit
                end do
    CLOSE(UNIT=10)
    		                   
! Open output file, calculate delta values and write         
	OPEN(UNIT=20,file=outfile,status='replace',&
    	access='sequential',action='write',&
        position='append',iostat=ierror)
        write(20,100)
        100 format(3X,"Ufinal-Uinitial")
        write(20,110)
        110 format(3X,"===============")
 
        		do k = 1, upoints
                	delta(k)=u1(k)-uo(k)
                	write(20,120) delta(k)
                	120 format(3X,F18.13)
               end do
    CLOSE(20)
    
! Deallocate arrays    
deallocate(uo,stat=alloc_err)
deallocate(u1,stat=alloc_err)                
deallocate(delta,stat=alloc_err) 
                
end program workfunction

!**********Subroutine for counting lines**************!

Subroutine countline(infile,nline)
implicit none
character(20),intent(in):: infile
integer,intent(out):: nline
integer::ierror
	
	nline=0
    OPEN(UNIT=10,file=infile,status='old',action='read',&
         iostat=ierror)
         do
           read(10,*,iostat=ierror)
           if (ierror /= 0) exit
		   nline=nline+1
         end do
    CLOSE(UNIT=10)

End subroutine countline

           

