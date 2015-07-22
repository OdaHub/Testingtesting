!-------------------------------------------------------
!Module to read variables from the input file
!-------------------------------------------------------
!BM---------------------------------------
Module read_input_file
implicit none

  !INTERFACE ASSIGMENTS: define operators on the objects(types)
  !read_inputparameter is a generic subroutine to read parameters
  !from an input file. The location of the parameter is traced by
  !its keyword. The subroutine can assign the value to strings, integer
  !values, double precision, logicals, vector of characters, vector of 
  !double precisions, and double precision matrices.

  !BI------------------------------------------------
  interface read_inputparameter 
    module procedure get_string, get_integer, get_dblprec, get_logical, &
                     get_dp_vector, get_int_vector, get_dp_matrix, &
                     get_char_vector, get_next_full_line
  end interface read_inputparameter
  !EI-----------------------------------------------


contains
!BS--------------------------------------------
subroutine get_string(string,filenm,keyword) 
use stringlengths
implicit none
character(LEN=*), intent(in)::filenm,keyword
character(LEN=*), intent(out)::string
   call search_keyword(filenm,keyword,keyname=string)
end subroutine get_string
!ES--------------------------------------------


!BS--------------------------------------------
subroutine get_integer(integ,filenm,keyword)
implicit none
character(LEN=*), intent(in)::filenm,keyword
integer, intent(out)::integ
   call search_keyword(filenm,keyword,keyinteger=integ)
end subroutine get_integer
!ES--------------------------------------------


!BS--------------------------------------------
subroutine get_dblprec(dbl,filenm,keyword)
implicit none
character(LEN=*), intent(in)::filenm,keyword
double precision, intent(out)::dbl
   call search_keyword(filenm,keyword,keydblprec=dbl)
end subroutine get_dblprec
!ES--------------------------------------------


!BS--------------------------------------------
subroutine get_logical(logic,filenm,keyword)
implicit none
character(LEN=*), intent(in)::filenm,keyword
logical, intent(out)::logic
   call search_keyword(filenm,keyword,keylog=logic)
end subroutine get_logical
!ES--------------------------------------------


!BS--------------------------------------------
subroutine get_char_vector(vector,filenm,keyword)
implicit none
character(LEN=*), intent(in) ::filenm,keyword
character(LEN=*), intent(out)::vector(:)
integer::Nrow

   Nrow=size(vector)
   call search_keyword(filenm,keyword,char_vector=vector,Nrow=Nrow)
end subroutine get_char_vector
!ES--------------------------------------------

!BS--------------------------------------------
subroutine get_dp_vector(vector,filenm,keyword)
implicit none
character(LEN=*), intent(in) ::filenm,keyword
double precision, intent(out)::vector(:)
integer::Nrow

   Nrow=size(vector)
   call search_keyword(filenm,keyword,dp_vector=vector,Nrow=Nrow)
end subroutine get_dp_vector
!ES--------------------------------------------

!BS--------------------------------------------
subroutine get_int_vector(vector,filenm,keyword)
implicit none
character(LEN=*), intent(in) ::filenm,keyword
integer, intent(out)::vector(:)
integer::Nrow

   Nrow=size(vector)
   call search_keyword(filenm,keyword,int_vector=vector,Nrow=Nrow)
end subroutine get_int_vector
!ES--------------------------------------------


!BS--------------------------------------------
subroutine get_dp_matrix(matrix,filenm,keyword)
implicit none
character(LEN=*), intent(in)::filenm,keyword
double precision, intent(out)::matrix(:,:)
integer::Nrow
Nrow=size(matrix(:,1))
   call search_keyword(filenm,keyword,dp_matrix=matrix,Nrow=Nrow)
end subroutine get_dp_matrix
!ES--------------------------------------------

!BF-----------------------------------------------
subroutine get_next_full_line(string,filenm,keyword,full_line)
implicit none
character(LEN=*), intent(in)::filenm,keyword
logical, intent(in)::full_line 
!only there to distuinguish between subroutine get_character
character(LEN=*), intent(out)::string
   call search_keyword(filenm,keyword,next_full_line=string)
end subroutine get_next_full_line
!EF-----------------------------------------------

!BS------------------------------------------------
subroutine search_keyword(filenm,keyword,keyname,keyinteger,keydblprec,keylog,&
char_vector,dp_vector,int_vector,dp_matrix,Nrow,next_full_line)
use defaultfile
implicit none
character(LEN=*), intent(in)::filenm,keyword
character(LEN=*), optional:: keyname
integer,optional::keyinteger
double precision, optional::keydblprec
logical, optional::keylog
integer, optional::Nrow
character(LEN=*), optional::char_vector(:)
double precision, optional::dp_vector(:),dp_matrix(:,:)
integer, optional::int_vector(:)
character(LEN=*), optional::next_full_line
character(len=LEN(keyword)+1)::word
integer::status
logical::keywordfound
character(LEN=max(LEN(filenm),LEN(inputdefault)))::file_array(2)
integer::i,j



file_array=(/filenm, inputdefault /)

do i=1,2

  keywordfound=.false.
  open(1,file=file_array(i),status="old")
    do

      read(1,*,iostat=status) word 
      if (status /=0 ) exit
      if (word==keyword) then

        keywordfound=.true.

        if (present(keyname)) then
          backspace 1
          read(1,*) word, keyname
          print *,trim(file_array(i))//":"//word//"="//trim(keyname)
        endif
  
        if (present(keyinteger)) then
          backspace 1
          read(1,*) word, keyinteger
          print *,trim(file_array(i))//":"//word//"=",keyinteger
        endif
  
        if (present(keydblprec)) then
          backspace 1
          read(1,*) word, keydblprec
          print *,trim(file_array(i))//":"//word//"=",keydblprec
        endif

        if (present(keylog)) then
          backspace 1
          read(1,*) word, keylog
          print *,trim(file_array(i))//":"//word//"=",keylog
        endif

        if (present(char_vector)) then
          print *,trim(file_array(i))//":"//word
          do j=1,Nrow
            read(1,*,iostat=status) char_vector(j)
            if (status /=0 ) then 
              print *,"ERROR: character-vector" 
              stop
            endif
            print *,j,char_vector(j) 
          enddo 
        endif 

        if (present(dp_vector)) then
          print *,trim(file_array(i))//":"//word
          do j=1,Nrow
            read(1,*,iostat=status) dp_vector(j)
            if (status /=0 ) then
              print *,"ERROR: dp-vector"
              stop
            endif
            print *,j,dp_vector(j)
          enddo
        endif

       if (present(int_vector)) then
          print *,trim(file_array(i))//":"//word
          do j=1,Nrow
            read(1,*,iostat=status) int_vector(j)
            if (status /=0 ) then
              print *,"ERROR: int-vector"
              stop
            endif
            print *,j,int_vector(j)
          enddo
        endif

        if (present(dp_matrix)) then
          print *,trim(file_array(i))//":"//word
          do j=1,Nrow
            read(1,*,iostat=status) dp_matrix(j,:)
            if (status /=0 ) then
              print *,"ERROR: dp-matrix"
              stop
            endif
            print *,j,dp_matrix(j,:)
          enddo
        endif

        if (present(next_full_line)) then
          print *,trim(file_array(i))//":"//word
          read(1,'(A200)',iostat=status) next_full_line 
          if (status /=0 ) then
            print *,"ERROR: next_full_line"
            stop
          endif
          if (next_full_line(200:200)/=" ") then
            print *,"directory name too long > 200 characters"
            stop
          endif
          print *,trim(next_full_line)
        endif
       
        exit  !keyword was found

      endif 

    enddo
  close(1)
  if (keywordfound) exit !no need to search in defaul-inputfile
enddo  

if (.not.keywordfound) then
  print *,"ERROR"
  print *,"keyword:",keyword," was not found in input or default file"
  stop
endif

end subroutine search_keyword
!ES------------------------------------------------


end Module read_input_file
!EM-----------------------------------------
