/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
//
// GOAL: Implement a 1D communication scheme along
//       X axis with blocking communications.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex1(lbm_comm_t * comm, int total_width, int total_height)
{
	int rank;
	int comm_size;
	int local_width_without_ghost;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	comm->nb_x = comm_size;
	comm->nb_y = 1;

	comm->rank_x = rank;
	comm->rank_y = 0;

	if (total_width % comm->nb_x != 0)
		MPI_Abort(MPI_COMM_WORLD, -1);

	local_width_without_ghost = total_width / comm->nb_x;

	comm->width = local_width_without_ghost + 2;
	comm->height = total_height + 2;

	comm->x = comm->rank_x * local_width_without_ghost;
	comm->y = 0;

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex1(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 1D communication with blocking MPI functions (MPI_Send & MPI_Recv)
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[DIRECTIONS] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);
}
