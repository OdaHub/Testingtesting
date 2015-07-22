module output
use shell
implicit none

contains

!-----------------------------------------
subroutine makeoutput
use var 
implicit none
integer::i
character(LEN=LSTR)::TISMOLtex,TISMOLdvi,TISMOLps,TISMOLpdf,string
  TISMOLtex="TISMOL.tex"
  TISMOLdvi="TISMOL.dvi"
  TISMOLps="TISMOL.ps"
  TISMOLpdf="TISMOL.pdf"
  open(1,file=TISMOLtex)
  call make_header(1) 
  
  string="\includegraphics[angle=-90,width=5.6cm]{"
  do i=1,NPSF
    write(1,*) trim(string)//trim(psfiles(i))//"}"
  enddo
  write(1,*) "\clearpage"
  write(1,*) "\input{summary}"
  write(1,*) "\end{document}"
  close(1)
  
  write(11,*)  "latex "//trim(TISMOLtex)
  write(11,*)  "dvips -o -t a4 "//trim(TISMOLdvi)
  write(11,*)  "ps2pdf "//trim(TISMOLps)
  write(11,*) "gv "//trim(TISMOLps)//" & "
  close(11)
  !print *,"make ""fexc.sh"" executable by chmod 770 ""fexc.sh"""
  !print *,"followed by ./fexc.sh"
  call do_shell("chmod 770 fexc.sh")
  call do_shell("./fexc.sh")
 
end subroutine makeoutput
!ES-----------------------------------------------

!BS---------------------------------------------
subroutine make_header(iU)
implicit none
integer, intent(in)::iU
write(iU,*) "\documentclass[12pt,english]{article}"
write(iU,*) "\usepackage[dvips]{graphicx}"
write(iU,*) "\setlength{\textwidth}{17cm}"
write(iU,*) "\setlength{\textheight}{31cm}"
write(iU,*) "\addtolength{\evensidemargin}{-3cm}"
write(iU,*) "\addtolength{\oddsidemargin}{-3cm}"
write(iU,*) "\addtolength{\topmargin}{-3.0cm}"
write(iU,*) "\pagestyle{empty}"
write(iU,*) "\begin{document}"
write(iU,*) "\noindent"
end subroutine make_header
!ES------------------------------------------------
end module output
