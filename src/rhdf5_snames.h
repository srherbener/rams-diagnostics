/*
!
! Copyright (C) 1991-2004  ; All Rights Reserved ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 2.2.0
!###########################################################################
*/

#if defined(SUN) || defined(ALPHA) || defined(SGI) || defined (PC_LINUX1) || defined(NEC_SX)

#define rh5f_open            rh5f_open_
#define rh5f_close           rh5f_close_
#define rh5f_create          rh5f_create_
#define rh5l_exists          rh5l_exists_
#define rh5g_create          rh5g_create_
#define rh5g_open            rh5g_open_
#define rh5g_close           rh5g_close_
#define rh5d_open            rh5d_open_
#define rh5d_write           rh5d_write_
#define rh5d_setup_and_write rh5d_setup_and_write_
#define rh5d_setup_and_read  rh5d_setup_and_read_
#define rh5d_close           rh5d_close_
#define rh5t_close           rh5t_close_
#define rh5s_close           rh5s_close_
#define rh5p_close           rh5p_close_
#define rh5a_write_anyscalar rh5a_write_anyscalar_
#define rh5a_read_anyscalar  rh5a_read_anyscalar_

#endif

#if defined(CRAY)

#define rh5f_open            RH5F_OPEN
#define rh5f_close           RH5F_CLOSE
#define rh5f_create          RH5F_CREATE
#define rh5l_exists          RH5L_EXISTS
#define rh5g_create          RH5G_CREATE
#define rh5g_open            RH5G_OPEN
#define rh5g_close           RH5G_CLOSE
#define rh5d_open            RH5D_OPEN
#define rh5d_write           RH5D_WRITE
#define rh5d_setup_and_write RH5D_SETUP_AND_WRITE
#define rh5d_setup_and_read  RH5D_SETUP_AND_READ
#define rh5d_close           RH5D_CLOSE
#define rh5t_close           RH5T_CLOSE
#define rh5s_close           RH5S_CLOSE
#define rh5p_close           RH5P_CLOSE
#define rh5a_write_anyscalar RH5A_WRITE_ANYSCALAR
#define rh5a_read_anyscalar  RH5A_READ_ANYSCALAR

#endif

