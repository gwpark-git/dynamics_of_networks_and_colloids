
#ifndef CONNECTIVITY_H
#define CONNECTIVITY_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <mkl.h>
#include "matrix.h"
#include "read_file_condition.h"

class CONNECTIVITY
{
 public:
  MATRIX *HASH;
  MKL_LONG *TOKEN;
  
  MKL_LONG Np; // number of particles (particle == bead == micelle)
  MKL_LONG Mc; // maximally avaliable connections
  
  CONNECTIVITY()
    {
      std::cout << "ERR:: Class CONNECTIVITY must have argument\n";
    }
  CONNECTIVITY(MKL_LONG number_of_particles, MKL_LONG maximum_connections);
  CONNECTIVITY(COND& given_condition);
  virtual ~CONNECTIVITY()
    {
      mkl_free(HASH);
      mkl_free(TOKEN);
    }

  MKL_LONG read_exist_hash(const char* fn_hash);
  
  MKL_LONG check_valid(MKL_LONG index_particle, MKL_LONG index_target);
};

class CHAIN_NODE
{
  /* 
     This class is designed for suggesting node for connectivity information based on single chain subject.
     Basically, its data structure is given by linked-list, and the index itself is related with the array of object.
   */
 public:
  /* MKL_LONG &HEAD; */
  /* MKL_LONG &TAIL; */
  MKL_LONG HEAD, TAIL; // the design is changed according to the design interface.

  MKL_LONG& index(MKL_LONG flag_HEAD_TAIL)
    {
      /*
	It prevent the complicate interface design which contains duplicated definition.
       */
      if (flag_HEAD_TAIL)
	return TAIL;
      return HEAD;
      // note that the reference variable can omit the & symbol when return its values.
    }
  /* MKL_LONG &HEAD() */
  /*   { */
  /*     return index[0]; */
  /*   } */
  /* MKL_LONG &TAIL() */
  /*   { */
  /*     return index[1]; */
  /*   } */
  /* MKL_LONG **index; */
  CHAIN_NODE()
    {
      /* index = (MKL_LONG**) mkl_malloc(2*sizeof(MKL_LONG*), BIT); */
      /* index[0] = &HEAD; */
      /* index[1] = &TAIL; */
      /* HEAD = index[0]; */
    /* TAIL = index[1]; */
    }
  virtual ~CHAIN_NODE(){
    /* mkl_free(index); */
    /* mkl_free(index); */
    /* HEAD = NULL; */
    /* TAIL = NULL; */
  }
};

/* class TRACKING_CHAIN */
/* { */
/*  public: */
/*   /\* Originally, the following scheme is related with the HEAD and TAIL*\/ */
/*   /\* However, the new scheme is used united approach based on index for chain ends rather than chain *\/ */
/*   /\* The rule is simple: using 2*Ntot as index and 0~Ntot-1 is related with the head and left are tail *\/ */
/*   /\* /\\* the HEAD and TAIL is enough for producing all the information *\\/ *\/ */
/*   /\* /\\* However, it makes large overhead during simulation *\\/ *\/ */
/*   /\* /\\* Therefore, the new index table is applied based on the number of particles *\\/ *\/ */
/*   /\* MKL_LONG *HEAD; *\/ */
/*   /\* MKL_LONG *TAIL; *\/ */
/*   MKL_LONG *CE_ATTACH; */
    
/*   MKL_LONG **C_INDEX; */
/*   MKL_LONG *C_TOKEN; */
  
/*   MKL_LONG initialization; */
/*   MKL_LONG N_tot; */
/*   MKL_LONG Np; */
/*   TRACKING_CHAIN() */
/*     { */
/*       std::cout << "ERR: Class TRACKING_CHAIN must have argument\n"; */
/*     } */
/*   TRACKING_CHAIN(MKL_LONG number_of_chains, MKL_LONG number_of_particles) */
/*     { */
/*       N_tot = number_of_chains; */
/*       Np = number_of_particles; */
/*       MKL_LONG Nc = N_tot/Np; */
/*       C_INDEX = (MKL_LONG**)mkl_malloc(Np*sizeof(MKL_LONG*), BIT); */
/*       C_TOKEN = (MKL_LONG*)mkl_malloc(Np*sizeof(MKL_LONG), BIT); */
/*       for(MKL_LONG i=0; i<Np; i++) */
/*         { */
/*           C_INDEX[i] = (MKL_LONG*)mkl_malloc(4*Nc*sizeof(MKL_LONG*), BIT); */
/*           C_TOKEN[i] = 0; */
/*         } */
      
/*       /\* HEAD = (MKL_LONG*)mkl_malloc(N_tot*sizeof(MKL_LONG), BIT); *\/ */
/*       /\* TAIL = (MKL_LONG*)mkl_malloc(N_tot*sizeof(MKL_LONG), BIT); *\/ */
/*       CE_ATTACH = (MKL_LONG*)mkl_malloc(2*N_tot*sizeof(MKL_LONG), BIT); */
/*       // 2*N_tot is of importance since the both of head and tail of subject chain are described by one approach */
/*       /\* for(MKL_LONG i=0; i<N_tot; i++) *\/ */
/*       /\*   { *\/ */
/*       /\*     HEAD[i] = i%N_p; // allocate the index based on the number of particles *\/ */
/*       /\*     TAIL[i] = HEAD[i]; *\/ */
/*       /\*   } *\/ */
/*       for(MKL_LONG i=0; i<2*N_tot; i++) */
/*         { */
/*           CE_ATTACH[i] = i%Np; */
/*           C_INDEX[CE_ATTACH[i], TOKEN[CE_ATTACH[i]]++] = i; // it allocates the index table (2Ntot*2Ntot) for the different chain ends */
/*         } */

/*       initialization = TRUE; */
/*     } */
/*   ~TRACKING_CHAIN() */
/*     { */
/*       if(initialization) */
/*         { */
/*           /\* mkl_free(HEAD); *\/ */
/*           /\* mkl_free(TAIL); *\/ */
/*           mkl_free(CE_ATTACH); */
/*           for(MKL_LONG i=0; i<Np; i++) */
/*             mkl_free(C_INDEX[i]); */
/*           mkl_free(C_INDEX); */
/*         } */
/*     } */


/*   MOV_CONNECTION(MKL_LONG index_chain_end, MKL_LONG index_target_particle) */
/*     { */
/*       for(MKL_LONG i=0; i<TOKEN[CE_ATTACH[index_chain_end]]; i++) */
        
/*       CE_ATTACH[index_chain_end] = index_target_particle; */
/*     } */
/*   /\* C_INDEX_ADD(MKL_LONG given_index_particle, MKL_LONG  *\/ */
  
/*   /\* SELECTING_CHAIN(MKL_LONG given_index) *\/ */
/*   /\*   { *\/ */
/*   /\*     for(MKL_LONG i=0; i<N_tot; i++) *\/ */
/*   /\*       { *\/ */
          
/*   /\*       } *\/ */
/*   /\*   } *\/ */
/* }; */


#endif
