
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_physcon.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> Contains physical constants required by the ice model.
module glimmer_physcon

  use glimmer_global, only : dp,sp

  implicit none
  
  save

  !TODO - Verify that all physical constants are mks. (I think they are.)
  !       Add a d0 to scyr, rhoi and grav

  real(dp),parameter :: scyr = 31556926.0        ! Number of seconds in a year (s). 
                                                 ! Note that this is for a 365.242 day year, and might need changing.

  real(dp),parameter :: pi = 3.1415926535897d0   !< Value of \f$\pi\f$.

  real(dp),parameter :: rhoi = 910.0             !< The density of ice (kg m<SUP>-3</SUP>)
  real(dp),parameter :: rhom = 3300.0d0          !< The density of magma(?) (kg m<SUP>-3</SUP>)

  real(dp),parameter :: rhoo = 1028.0d0          !< The density of the ocean (kg m<SUP>-3</SUP>)
  real(dp),parameter :: rhow = 1000.0d0          !< The density of fresh water (kg m<SUP>-3</SUP>)
  real(dp),parameter :: rhos = 2600.0d0          !*FD The density of solid till (kg m$^{-3}$)
  real(dp),parameter :: f = - rhoo / rhoi

  real(dp),parameter :: grav = 9.81              !< The acceleration due to gravity (m s<SUP>-2</SUP>)

  integer, parameter :: gn = 3                   !< The power dependency of Glenn's flow law.

  real(dp),parameter :: arrmlh = 1.733d3         !< Constant of proportionality in Arrhenius relation
                                                 !< in \texttt{patebudd}, for \f$T^{*}\geq263\f$K.
                                                 !< (Pa<SUP>-3</SUP> s<SUP>-1</SUP>) 
  real(dp),parameter :: arrmll = 3.613d-13       !< Constant of proportionality in Arrhenius relation
                                                 !< in \texttt{patebudd}, for \f$T^{*}<263\f$K.
                                                 !< (Pa<SUP>-3</SUP> s<SUP>-1</SUP>) 
  real(dp),parameter :: gascon = 8.314d0         !< The gas ideal constant \f$R\f$ (J mol<SUP>-1</SUP> K<SUP>-1</SUP>)
  real(dp),parameter :: actenh = 139.0d3         !< Activation energy in Glenn's flow law for \f$T^{*}\geq263\f$K. (J mol<SUP>-1</SUP>)
  real(dp),parameter :: actenl = 60.0d3          !< Activation energy in Glenn's flow law for \f$T^{*}<263\f$K. (J mol<SUP>-1</SUP>)

  real(dp),parameter :: shci = 2009.0d0          !< Specific heat capacity of ice (J kg<SUP>-1</SUP> K<SUP>-1</SUP>)
  real(dp),parameter :: lhci = 335.0d3           !< Latent heat of melting of ice (J kg<SUP>-1</SUP>) 
  real(dp),parameter :: coni = 2.1d0             !< Thermal conductivity of ice (W m<SUP>-1</SUP> K<SUP>-1</SUP>)

  real(dp),parameter :: pmlt = 9.7456d-8         !< Factor for dependence of melting point on pressure (K Pa<SUP>-1</SUP>)
  real(dp),parameter :: trpt = 273.15d0          !< Triple point of water (K)

end module glimmer_physcon
