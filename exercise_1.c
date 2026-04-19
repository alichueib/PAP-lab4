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
	int size;
	int local_width_without_ghost;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	comm->nb_x = size;
	comm->nb_y = 1;
	comm->rank_x = rank;
	comm->rank_y = 0;

	//the global width must be divisible by the number of X sub-domains
	if (total_width % size != 0)
		MPI_Abort(MPI_COMM_WORLD, -1); //we tried using return but it does not work, the program just hangs (deadlock) instead of exiting

	local_width_without_ghost = total_width /size;
	comm->width = local_width_without_ghost + 2;
	comm->height = total_height + 2;
	comm->y = 0;
	comm->x = rank * local_width_without_ghost;

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex1(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	int left;
	int right;
	const int column_size = comm->height * DIRECTIONS;
	left = (comm->rank_x == 0) ? MPI_PROC_NULL : comm->rank_x - 1;
	right = (comm->rank_x == comm->nb_x - 1) ? MPI_PROC_NULL : comm->rank_x + 1;

	if(comm->rank_x == 0){
		left = MPI_PROC_NULL;}
	else{
		left = comm->rank_x - 1;}

	if(comm->rank_x == comm->nb_x - 1){
		right = MPI_PROC_NULL;}
	else{
		right = comm->rank_x + 1;
	}

	if (left != MPI_PROC_NULL)
		MPI_Recv(lbm_mesh_get_cell(mesh, 0, 0), column_size, MPI_DOUBLE, left, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	if (right != MPI_PROC_NULL)
		MPI_Send(lbm_mesh_get_cell(mesh, comm->width - 2, 0),column_size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);

	//now we should send to the left and receive from the right
	if (right != MPI_PROC_NULL)MPI_Recv(lbm_mesh_get_cell(mesh, comm->width - 1, 0), column_size, MPI_DOUBLE,right,1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	if (left != MPI_PROC_NULL)MPI_Send(lbm_mesh_get_cell(mesh, 1, 0),column_size, MPI_DOUBLE, left,1,MPI_COMM_WORLD);
}
