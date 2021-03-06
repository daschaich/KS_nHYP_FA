// -----------------------------------------------------------------
// Application-dependent routine for writing gauge info file
// For SU(3) lattice doubling

// This file is an ASCII companion to the gauge configuration file
// and contains information about the action used to generate it.
// This information is consistently written in the pattern
//     keyword value
// or
//     keyword[n] value1 value2 ... valuen
// where n is an integer.
//
// Possible keywords are listed in io_lat.h
#include "double_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Write optional entries in the ASCII info file
// Call from one of the lattice output routines in io_lat4.c
// File has already been opened and the required magic number,
// time stamp and lattice dimensions have already been written
void write_appl_gauge_info(FILE *fp) {
  write_gauge_info_item(fp, "gauge.description", "%s",
                        "\"Doubled in all four directions\"",
                        0, 0);
}
// -----------------------------------------------------------------
