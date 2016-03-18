// -----------------------------------------------------------------
// How does this relate to ../generic/gauge_stuff.c?
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void char_num(int dig[], int *chr, int *ch, int length) {
  int j, bdig[10], tenl = 1, newv, old;

  *ch = 0;
  for (j = 0; j < length - 1; j++)
    tenl *= 10;

  *chr = dig[length - 1];
  for (j = length - 2; j >= 0; j--)
    *chr = *chr * 10 + dig[j];

  // Forward
  old = *chr;
  for (j = length - 1; j >= 1; j--) {
    newv = old - tenl * dig[j];
    newv = newv * 10 + dig[j];
    if (newv < *chr)
      *chr = newv;

    old = newv;
  }

  // Backward
  for (j = 0; j < length; j++)
    bdig[j] = 7 - dig[length - j - 1];

  old = bdig[length - 1];
  for (j = length - 2; j >= 0; j--)
    old = old * 10 + bdig[j];

  if (old < *chr)
    *chr = old;

  for (j = length - 1; j >= 1; j--) {
    newv = old - tenl * bdig[j];
    newv = newv * 10 + bdig[j];
    if (newv < *chr)
      *chr = newv;

    old = newv;
  }

  if (*chr < tenl)
    *ch = 1;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void make_loop_term() {
  int dirs[10], dir = 0, ln, iloop, l_num = 0, length, k, irep;

  // Suppress annoying compiler warnings
  for (k = 0; k < 10; k++)
    dirs[k] = 0;

  // Loop over all the different contributions to plaquettes
  for (iloop = 0; iloop < nloop; iloop++) {
    // node0_printf("iloop = %d\n", iloop);
    length = loop_length[iloop];
    // Loop over rotations and reflections
    for (ln = 0; ln < loop_num[iloop]; ln++) {
      if (loop_ch[iloop][ln] == 1) {
        // Set up dirs and sign
        for (k = 0; k < length; k++) {
          if (loop_table[iloop][ln][k] < 4)
            dirs[k] = (dir + loop_table[iloop][ln][k]) % 4;
          else
            dirs[k] = (7 + dir - loop_table[iloop][ln][k]) % 4;
        }

        // Loop over position in path
        for (k = 0; k < length; k++) {
          if (dirs[k] == dir) {
            for (irep = 0; irep < nreps; irep++) {
              loop_term[l_num][irep] = loop_coeff[iloop][irep];
//              node0_printf("%d %d %.4g\n",
//                           l_num, irep, loop_term[l_num][irep]);
            }
            l_num++;
          }
        } // End loop over location in path
      }
    } // End loop over rotations and reflections
  } // End loop over different contributions to plaquettes
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void make_loop_table2() {
  int perm[8], pp[8], ir[4], vec[10];
  int length, iloop, i, j, chr, count, flag;
  Real u0;
  void make_loop_term();

  static int loop_length_in[nloop] = {4, 6, 6, 6, 8, 8};
  static int loop_ind[nloop][10] = {{0, 1, 7, 6, 0, 0, 0, 0, 0, 0},
                                    {0, 1, 1, 7, 6, 6, 0, 0, 0, 0},
                                    {0, 1, 2, 7, 6, 5, 0, 0, 0, 0},
                                    {0, 1, 2, 7, 5, 6, 0, 0, 0, 0},
                                    {0, 1, 0, 1, 7, 6, 7, 6, 0, 0},
                                    {0, 0, 1, 1, 7, 7, 6, 6, 0, 0}};

  for (j = 0; j < nloop; j++) {
    loop_num[j] = 0;
    loop_length[j] = loop_length_in[j];
  }

  // Set up the loop coefficients
  for (i = 0; i < nreps; i++) {
    for (j = 0; j < nloop; j++)
      loop_coeff[j][i] = 1;
  }

  // Plaquette
  loop_coeff[0][0]= 1;

  // a**2 improved Symanzik with funny couplings
  u0 = 0.868;
  loop_coeff[0][0] = 1;
  loop_coeff[1][0] = (0.6264 * log(u0) - 1) / (20.0 * u0 * u0);
  loop_coeff[2][0] = 0.04335 * log(u0) / (u0 * u0);

  for (iloop = 0; iloop < nloop; iloop++) {
    length = loop_length[iloop];
    count=0;
    // Permutations
    for (perm[0] = 0; perm[0] < 4; perm[0]++) {
      for (perm[1] = 0; perm[1] < 4; perm[1]++) {
        for (perm[2] = 0; perm[2] < 4; perm[2]++) {
          for (perm[3] = 0; perm[3] < 4; perm[3]++) {
            if (perm[0] != perm[1]
                && perm[0] != perm[2]
                && perm[0] != perm[3]
                && perm[1] != perm[2]
                && perm[1] != perm[3]
                && perm[2] != perm[3]) {

              // Reflections

              for (ir[0] = 0; ir[0] < 2; ir[0]++) {
                for (ir[1] = 0; ir[1] < 2; ir[1]++) {
                  for (ir[2] = 0; ir[2] < 2; ir[2]++) {
                    for (ir[3] = 0; ir[3] < 2; ir[3]++) {
                      for (j = 0; j < 4; j++) {
                        pp[j] = perm[j];
                        if (ir[j] == 1)
                          pp[j] = 7 - pp[j];

                        pp[7 - j] = 7 - pp[j];
                      }

                      // Create new vector
                      for (j = 0; j < length; j++)
                        vec[j] = pp[loop_ind[iloop][j]];

                      char_num(vec,&chr,&ch,length);
                      flag = 0;

                      // Check if it's a new set
                      for (j = 0; j < count; j++) {
                        if (chr == loop_char[j])
                          flag = 1;
                      }
                      if (flag == 0) {
                        loop_char[count] = chr;
                        loop_ch[iloop][count] = ch;
                        for (j = 0; j < length; j++)
                          loop_table[iloop][count][j] = vec[j];

                        count++;
                      }
                      loop_num[iloop] = count;
                    } // End loop over ir[3] reflections
                  } // End loop over ir[2] reflections
                } // End loop over ir[1] reflections
              } // End loop over ir[0] reflection
            }
          } // End loop over perm[3] permutations
        } // End loop over perm[2] permutations
      } // End loop over perm[1] permutations
    } // End loop over perm[0] permutations
//    node0_printf("iloop, loop_num %d %d\n", iloop, loop_num[iloop]);
  } // End loop over iloop

//  for (iloop = 0; iloop < nloop; iloop++) {
//    for (count = 0; count < loop_num[iloop]; count++) {
//      node0_printf(" %d %d %d %d ",
//                   iloop, count, loop_char[count], loop_ch[iloop][count]);
//      for (j=0;j<length;j++)
//        node0_printf(" %d",loop_table[iloop][count][j]);
//
//      node0_printf("\n");
//    }
//  }

  // Print out the loop coefficients
  node0_printf("Loop coefficients: nloop, rep, loop_coeff, multiplicity\n");
  for (i = 0; i < nreps; i++) {
    for (j = 0; j < nloop; j++) {
      node0_printf("%d %d %.4g %d\n",
                   j, i, loop_coeff[j][i], loop_num[j]);
    }
  }
  node0_printf("\n");
  make_loop_term();
}
// -----------------------------------------------------------------
