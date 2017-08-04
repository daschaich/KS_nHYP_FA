/******************** ks_utilities.c *************************************/
/* Miscellaneous (mostly) KS utilities for staggered code */
#include "generic_ks_includes.h"

/*--------------------------------------------------------------------*/
vector *create_v_field(void){
  vector *cv;
  
  cv = (vector *)
    malloc(sizeof(vector)*sites_on_node);
  
  if(cv == NULL){
    printf("create_cv_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }
  return cv;
}

/*--------------------------------------------------------------------*/
void clear_v_field(vector *v){

  int i;
  site *s;

  FORALLSITES(i,s){
      clearvec( &v[i] );
  }
}

/*--------------------------------------------------------------------*/
void destroy_v_field(vector *v){
  free(v);
}

/*--------------------------------------------------------------------*/
ks_prop_field create_ksp_field(void){
  ks_prop_field ksp;
  int color;
  
  ksp = (ks_prop_field) malloc(3*sizeof(vector *));
  if(ksp == NULL){
    printf("create_ksp_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }

  for(color= 0; color < 3; color++)
    ksp[color] = create_v_field();
  
  return ksp;
}

/*--------------------------------------------------------------------*/
void clear_ksp_field(ks_prop_field ksp){
  int color;
  
  for(color= 0; color < 3; color++)
    clear_v_field(ksp[color]);
}

/*--------------------------------------------------------------------*/
void copy_ksp_field(ks_prop_field kspcopy, ks_prop_field ksp){
  int color, i;
  site *s;

  for(color = 0; color < 3; color++){
    FORALLSITES(i,s){
      kspcopy[color][i] = ksp[color][i];
    }
  }
}

/*--------------------------------------------------------------------*/
ks_prop_field create_ksp_field_copy(ks_prop_field k){
  ks_prop_field ksp;

  ksp = create_ksp_field();
  copy_ksp_field(ksp, k);

  return ksp;
}

/*--------------------------------------------------------------------*/
void destroy_ksp_field(ks_prop_field ksp){
  int color;

  if(ksp == NULL)return;
  for(color = 0; color < 3; color++){
    destroy_v_field(ksp[color]);
  }
  free(ksp);
}
