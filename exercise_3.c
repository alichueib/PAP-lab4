/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement non-blocking 1D communication scheme
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
// NEW:
//     - >>> Non-blocking communications <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex3(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex3(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int left_rank;
	int right_rank;
	int req_count = 0;
	MPI_Request requests[4];
	const int column_size = comm->height * DIRECTIONS;

	left_rank = (comm->rank_x == 0) ? MPI_PROC_NULL : comm->rank_x - 1;
	right_rank = (comm->rank_x == comm->nb_x - 1) ? MPI_PROC_NULL : comm->rank_x + 1;

	// Receive the left neighbor interior column into our left ghost column (x = 0).
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, 0, 0),
		column_size,
		MPI_DOUBLE,
		left_rank,
		1,
		MPI_COMM_WORLD,
		&requests[req_count++]
	);

	// Receive the right neighbor interior column into our right ghost column (x = width - 1).
	MPI_Irecv(
		lbm_mesh_get_cell(mesh, comm->width - 1, 0),
		column_size,
		MPI_DOUBLE,
		right_rank,
		0,
		MPI_COMM_WORLD,
		&requests[req_count++]
	);

	// Send our first interior column (x = 1) to the left neighbor.
	MPI_Isend(
		lbm_mesh_get_cell(mesh, 1, 0),
		column_size,
		MPI_DOUBLE,
		left_rank,
		0,
		MPI_COMM_WORLD,
		&requests[req_count++]
	);

	// Send our last interior column (x = width - 2) to the right neighbor.
	MPI_Isend(
		lbm_mesh_get_cell(mesh, comm->width - 2, 0),
		column_size,
		MPI_DOUBLE,
		right_rank,
		1,
		MPI_COMM_WORLD,
		&requests[req_count++]
	);

	MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
}
