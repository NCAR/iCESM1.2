Dear user,

Please note that line 3807 in ice_therm_vertical.F90:

melts(i,j) = melts(i,j) - dhs

has been uncommented in this particular version of iCESM.  The file itself can
be found here:

<iCESM_tag_directory>/models/ice/cice/src/source/ice_therm_vertical.F90

Having this line of code active results in a better water isotope simulation.
However, it also can change the actual physical climate, and thus result in
climates that don't match non-isotopic CICE4 results.  Thus if one needs the
isotopic simulations to match a standard CICE4 or CESM1 simulation, then this
line of code should be commented out.

Any questions or concerns with this particular issue can be directed to Jesse Nusbaumer <jesse.nusbaumer@nasa.gov>.

Good luck with the simulations, and have a great day!

