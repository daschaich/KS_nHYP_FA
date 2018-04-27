// -----------------------------------------------------------------
// Compute SU(3) field strength tensor
// Need to define six field strength indices
// Use tempmat and tempmat2 for temporary storage

// This cartoon shows how the plaquettes are calculated
// The path begins and ends at the 'O' at the corner
// F_munu is the sum of these plaquettes minus their adjoints

//  ^         --------<--------       --------<--------
//  |dir2     |               |       |               |
//  |         |               |       |               |
//  |         |               |       |               |
//  |         |               ^       |               ^
//  ------>   |               |       |               |
//    dir     |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            -------->-------O       O------->--------
//
//            --------<-------O       O-------<--------
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            |               ^       |               ^
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            -------->--------       -------->--------

// Convention: try to use gen_pt[0] and mtag for links in direction dir
// These are gathered from +/- dir2

#include "generic_includes.h"

#define LINK_OFFSET(dir) link_src + sizeof(matrix) * (dir)
#define LINK(dir) (((matrix *)F_PT(s, link_src))[dir])
#define FS(component) (((matrix *)F_PT(s, field_dest))[component])

// link_src is offset for matrix link[4] in site struct
// field_dest is offset for matrix fieldstrength[6] in site struct
void make_field_strength(field_offset link_src, field_offset field_dest) {
  register int i, component, dir = -99, dir2 = -99;
  register site *s;
  int j;
  complex cc;
  matrix tmat, tmat2;
  msg_tag *mtag, *mtag2;

  for (component = FS_XY; component <= FS_ZT; component++) {
    switch(component) {
      case FS_XY: dir = XUP; dir2 = YUP; break;
      case FS_XZ: dir = XUP; dir2 = ZUP; break;
      case FS_YZ: dir = YUP; dir2 = ZUP; break;
      case FS_XT: dir = XUP; dir2 = TUP; break;
      case FS_YT: dir = YUP; dir2 = TUP; break;
      case FS_ZT: dir = ZUP; dir2 = TUP; break;
    }

    // +dir +dir2 plaquette
    mtag  = start_gather_site(LINK_OFFSET(dir), sizeof(matrix),
                              dir2, EVENANDODD, gen_pt[0]);
    mtag2 = start_gather_site(LINK_OFFSET(dir2), sizeof(matrix),
                              dir, EVENANDODD, gen_pt[1]);

    wait_gather(mtag);
    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_nn(&LINK(dir), (matrix *)(gen_pt[1][i]), &tmat);
      mult_na(&tmat, (matrix *)(gen_pt[0][i]), &tmat2);
      mult_na(&tmat2, &LINK(dir2), &tmat);
      adjoint(&tmat, &tmat2);
      sub_matrix(&tmat, &tmat2, &FS(component));
    }
    cleanup_gather(mtag2);

    // -dir +dir2 plaquette
    // Reuse link[dir] gather from dir2 corresponding to mtag
    FORALLSITES(i, s) {
      mult_an(&LINK(dir2), &LINK(dir), &tmat);
      mult_an((matrix *)(gen_pt[0][i]), &tmat, &(tempmat[i]));
    }
    mtag2 = start_gather_field(tempmat, sizeof(matrix),
                               OPP_DIR(dir), EVENANDODD, gen_pt[1]);

    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_nn(&LINK(dir2), (matrix *)(gen_pt[1][i]), &tmat);
      adjoint(&tmat, &tmat2);
      sum_matrix(&tmat, &FS(component));
      dif_matrix(&tmat2, &FS(component));
    }
    cleanup_gather(mtag);
    cleanup_gather(mtag2);

    // -dir -dir2 plaquette
    mtag = start_gather_site(LINK_OFFSET(dir), sizeof(matrix),
                             OPP_DIR(dir), EVENANDODD, gen_pt[0]);
    mtag2 = start_gather_site(LINK_OFFSET(dir2), sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

    wait_gather(mtag);
    wait_gather(mtag2);
    FORALLSITES(i,s){
      mult_nn((matrix *)(gen_pt[0][i]), &LINK(dir2), &(tempmat[i]));
      mult_nn((matrix *)(gen_pt[1][i]), &LINK(dir), &(tempmat2[i]));
    }
    cleanup_gather(mtag);
    cleanup_gather(mtag2);

    mtag = start_gather_field(tempmat, sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[0]);
    mtag2 = start_gather_field(tempmat2, sizeof(matrix),
                               OPP_DIR(dir), EVENANDODD, gen_pt[1]);

    wait_gather(mtag);
    wait_gather(mtag2);
    FORALLSITES(i,s){
      mult_an((matrix *)(gen_pt[1][i]), (matrix *)(gen_pt[0][i]), &tmat);
      adjoint(&tmat, &tmat2);
      sum_matrix(&tmat, &FS(component));
      dif_matrix(&tmat2, &FS(component));
    }
    cleanup_gather(mtag);
    cleanup_gather(mtag2);

    // +dir -dir2 plaquette
    mtag2 = start_gather_site(LINK_OFFSET(dir2), sizeof(matrix),
                              dir, EVENANDODD, gen_pt[1]);

    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_an(&LINK(dir2), &LINK(dir), &tmat);
      mult_nn(&tmat, (matrix *)(gen_pt[1][i]), &tempmat[i]);
    }
    cleanup_gather(mtag2);

    mtag = start_gather_field(tempmat, sizeof(matrix),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[0]);
    wait_gather(mtag);
    FORALLSITES(i,s){
      mult_na((matrix *)(gen_pt[0][i]), &LINK(dir), &tmat);
      adjoint(&tmat, &tmat2);
      sum_matrix(&tmat, &FS(component));
      dif_matrix(&tmat2, &FS(component));
    }
    cleanup_gather(mtag);

    // Make traceless
    FORALLSITES(i, s) {
      cc = trace(&FS(component));
      CDIVREAL(cc, 3.0, cc);
      for (j = 0; j < 3; j++)
        CSUB(FS(component).e[j][j], cc, FS(component).e[j][j]);
    }
  }
}
// -----------------------------------------------------------------
