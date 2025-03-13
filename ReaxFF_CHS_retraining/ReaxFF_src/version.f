********************************************************************** 
********************************************************************** 
*                                                                    *
*     REAXFF Reactive force field program                            *
*                                                                    *
*     Developed and written by Adri van Duin, duin@wag.caltech.edu   *
*                                                                    *
*     Copyright (c) 2001-2010 California Institute of Technology     *
*                                                                    *
*     This is an open-source program. Feel free to modify its        *
*     contents. Please keep me informed of any useful modification   *
*     or addition that you made. Please do not distribute this       *
*     program to others; if people are interested in obtaining       *
*     a copy of this program let them contact me first.              *
*                                                                    *
********************************************************************** 
      subroutine version   
      include 'cbka.blk'
      character*20 qhulp
      qhulp='Compiled on:'
      write (*,100)
      write (*,110)qhulp
      return
  100 format ('ReaxFF version 2.0')
  110 format(a20,'Thu Jul 17 11:46:50 PDT 2008')
      end
