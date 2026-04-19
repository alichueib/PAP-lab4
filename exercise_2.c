/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement odd/even 1D blocking communication scheme 
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
// NEW:
//     - >>> Odd/even communication ordering <<<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex2(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex2(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int rank;
	int left;
	int right;
	const int column_size = comm->height *DIRECTIONS;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(comm->rank_x == 0){
		left = MPI_PROC_NULL;}
	else{
		left = comm->rank_x - 1;}

	if(comm->rank_x == comm->nb_x - 1){
		right = MPI_PROC_NULL;}
	else{
		right = comm->rank_x + 1;
	}
	if ((rank % 2) == 0)
	{
		//in case of even rank, they first receive the left-going column from the right neighbor
		MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, 0),column_size, MPI_DOUBLE, right,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		//then they receive the right-going column from the left neighbor into their left ghost column
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, 0),column_size, MPI_DOUBLE,right,1,MPI_COMM_WORLD);

		//then they send their first interior column to the left neighbor, which will be received into its right ghost column
		MPI_Send(lbm_mesh_get_cell(mesh, 1, 0), column_size,MPI_DOUBLE, left,0,MPI_COMM_WORLD);
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, 0),column_size, MPI_DOUBLE, left,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else
	{
		//in case of odd rank, they first send their first interior column to the left neighbor, which will be received into its right ghost column
		MPI_Send(lbm_mesh_get_cell(mesh, 1, 0),column_size, MPI_DOUBLE, left,0,MPI_COMM_WORLD);

		//then they receive the right-going neighbot
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, 0),column_size, MPI_DOUBLE, left,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		//then receive the left-going column from the right neighbor
		MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, 0),column_size, MPI_DOUBLE,right,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE
		);
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, 0), column_size, MPI_DOUBLE,right,1,MPI_COMM_WORLD);
	}
}
