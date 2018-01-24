/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:11:45
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-04 16:49:42
 */

#include "pfHome.h"

MPI_Datatype MPI_VECTOR;
MPI_Datatype MPI_ATOM;

#define MAX_MPI_COMPONENTS 30
#define CHECK_RETURN(a)                                            \
  do {                                                             \
    int r = a;                                                     \
    if (r != MPI_SUCCESS) {                                        \
      printf("Error in %s on line %d: %d", __FILE__, __LINE__, r); \
      return PFERROR;                                            \
    }                                                              \
  } while (0);

void pfHome::pfMPIinit(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &gsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mrank);
  createMPIdataTypes();
  return;
}

void pfHome::pfMPIfinalize() { MPI_Finalize(); }

int pfHome::createMPIdataTypes() {
  /***************************************************************
   *  MPI vector
   ***************************************************************/
  CHECK_RETURN(MPI_Type_contiguous(3, MPI_DOUBLE, &(MPI_VECTOR)));
  CHECK_RETURN(MPI_Type_commit(&MPI_VECTOR));

  int size_a = 0;
  int size_b = 0;

  int data_len[MAX_MPI_COMPONENTS];
  MPI_Datatype data_type[MAX_MPI_COMPONENTS];
  MPI_Aint data_size[MAX_MPI_COMPONENTS];

  /***************************************************************
   *  MPI atomMPI
   ***************************************************************/
  data_len[size_a] = 1;
  data_type[size_a++] = MPI_INT;

  data_len[size_a] = 1;
  data_type[size_a++] = MPI_VECTOR;

  data_len[size_a] = 1;
  data_type[size_a++] = MPI_VECTOR;

  data_len[size_a] = 1;
  data_type[size_a++] = MPI_DOUBLE;

  AtomMPI oneAtom;
  CHECK_RETURN(MPI_Get_address(&oneAtom.id, &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&oneAtom.pst, &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&oneAtom.frc, &data_size[size_b++]));
  CHECK_RETURN(MPI_Get_address(&oneAtom.absfrc, &data_size[size_b++]));

  if (size_a != size_b) {
    fprintf(stderr,
            "Sizes of item array and address array differ"
            "in %s:%d",
            __FILE__, __LINE__);
    MPI_Finalize();
    exit(PFERROR);
  }
  /* calculate the displacements from address to get size */
  for (int i = 1; i < size_a; ++i) data_size[i] -= data_size[0];
  data_size[0] = 0;

  CHECK_RETURN(MPI_Type_create_struct(size_a, data_len, data_size, data_type,
                                      &MPI_ATOM));
  CHECK_RETURN(MPI_Type_commit(&MPI_ATOM));
  return PFSUCCESS;
}
