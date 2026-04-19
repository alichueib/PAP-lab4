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
	int left;
	int right;
	int req_count = 0;
	MPI_Request requests[4]; const int column_size = comm->height *DIRECTIONS;
	if(comm->rank_x == 0){
		left = MPI_PROC_NULL;}
	else{
		left = comm->rank_x - 1;}

	if(comm->rank_x == comm->nb_x - 1){
		right = MPI_PROC_NULL;}
	else{
		right = comm->rank_x + 1;
	}

	//receive the left neighbor interior column into our left ghost column (x = 0)
	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, 0), column_size,MPI_DOUBLE,left,1, MPI_COMM_WORLD, &requests[req_count++]);
	//receive the right neighbor interior column into our right ghost column (x = width - 1)
	MPI_Irecv(lbm_mesh_get_cell(mesh, comm->width - 1, 0), column_size,MPI_DOUBLE,right,0,MPI_COMM_WORLD, &requests[req_count++]);
	//send our first interior column (x = 1) to the left neighbor
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, 0),column_size, MPI_DOUBLE, left,0,MPI_COMM_WORLD, &requests[req_count++]);
	//send our last interior column (x = width - 2) to the right neighbor.
	MPI_Isend(lbm_mesh_get_cell(mesh, comm->width - 2, 0), column_size,MPI_DOUBLE,right, 1, MPI_COMM_WORLD, &requests[req_count++]);
	MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
}
