/* =============================================
 * Kostas Evangelou
 * mpi_heat2D_opt.c
 * Parallel Systems
 * MPI Optimization of heat2D problem
 * =============================================*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NXPROB      20                 /* x dimension of problem grid */
#define NYPROB      20                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */
#define DTAG        0                  /* message tag */
#define UTAG        1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE       -1                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */
#define UP          0
#define DOWN        1
#define LEFT        2
#define RIGHT       3

struct Parms {
    float cx;
    float cy;
} parms = {0.1, 0.1};

int main(int argc, char *argv[]) {
    void    inidat(), prtdat(), update_in(), update_out();
    float   ***u;                       /* array for grid */
    float   *temp[2];                   /* temp 1d array for fast memory allocation */
    int     taskid,                     /* this task's unique id */
            numworkers,                 /* number of worker processes */
            numtasks,                   /* number of tasks */
            averow, rows, offset, extra,/* for sending rows of data */
            dest, source,               /* to - from for message send-receive */
            left, right,                /* neighbor tasks */
            msgtype,                    /* for message types */
            rc, start, end,             /* misc */
            i, ix, iy, iz, it,          /* loop variables */
            dims[2];                    /* cart argument*/
            reorder = 1,                /* Ranking may be reordered (true) or not (false) (logical) */
            periods[2] = {0, 0}         /* Logical  array of size ndims specifying whether the grid is periodic */
            neighbors[4];               /* Neighbors indicator */
    MPI_Status status;
    double start_time = 0.0,
            end_time = 0.0;

/* First, find out my taskid and how many tasks are running */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    numworkers = numtasks;

    if ((NYPROB * NXPROB) % numworkers != 0) {
        printf("ERROR: Number of tasks must have remainder 0 when divided with %d x %d \n", NXPROB, NYPROB);
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    /* Mpi Cartesian Grid */

    MPI_Comm cart_comm;
    size = sqrt(numworkers);
    dims[0] = size;
    dims[1] = size;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);

    MPI_Barrier(cartcomm);
    MPI_Cart_shift(cartcomm, 0, 1, &neighbors[UP], &neighbors[DOWN]);
    MPI_Cart_shift(cartcomm, 1, 1, &neighbors[LEFT], &neighbors[RIGHT]);

    sub_table_dim = sqrt(NXPROB*NYPROB/numtasks);
    sub_x = sub_table_dim + 2;
    sub_y = sub_table_dim + 2;

    /* Declare MPI Types */
    MPI_Type_contiguous(sub_x, MPI_FLOAT, &row);
    MPI_Type_commit(&MPI_row);
    MPI_Type_vector(sub_y, 1, sub_ x, MPI_FLOAT, &column);
    MPI_Type_commit(&MPI_column);

    MPI_Request req[8];

    /* Allocate Memory for Table */
    u=(float***)malloc(2 * sizeof(float**));
    for(i=0;i<2;i++){
      u[i] = (float**)malloc(sub_x * sizeof(float*));
      temp[i] = (float*)malloc(sub_y*sub_x*sizeof(float));
      for(j=0;j<sub_x;j++){
        u[i][j] = temp[i] + (j*sub_y);
      }
    }

    /* Init table with 0s */
    for (iz=0; iz<2; iz++)
        for (ix=0; ix<sub_x; ix++)
            for (iy=0; iy<sub_y; iy++)
                u[iz][ix][iy] = 0.0;

    /* Set Random values to sub_table */
    inidat(sub_x, sub_y, taskid, u[0]);

    MPI_Barrier(cart_comm);
    start_time = MPI_Wtime();

#if (CONVERGENCE == 1)
    task_convergence = reduced_convergence = 0;
#endif

    for(step=0; step<=STEPS;steps++){
        /* Send and Receive asynchronous the shared values of neighbor */
        if (neighbors[UP] >= 0){
            MPI_Isend(table_u + iz*sub_x*sub_y + sub_y + 1, sub_table_dimention, MPI_FLOAT, neighbors[UP], DTAG, cartcomm, &req[0]);
            MPI_Irecv(table_u + iz*sub_x*sub_y + 1, sub_table_dimention, MPI_FLOAT, neighbors[UP], UTAG, cartcomm, &req[1]);
        }

        if (neighbors[DOWN] >= 0){
            MPI_Isend(table_u + iz*sub_x*sub_y + sub_table_dimention*sub_y + 1, sub_table_dimention , MPI_FLOAT, neighbors[DOWN], UTAG, cartcomm, &req[2]);
            MPI_Irecv(table_u + iz*sub_x*sub_y + (sub_table_dimention+1)*sub_y + 1, sub_table_dimention , MPI_FLOAT, neighbors[DOWN], DTAG, cartcomm, &req[3]);
        }

        if (neighbors[LEFT] >= 0){
            MPI_Isend(table_u + iz*sub_x*sub_y + sub_y + 1, 1, COL_INT, neighbors[LEFT], RTAG, cartcomm,&req[4]);
            MPI_Irecv(table_u + iz*sub_x*sub_y + sub_y, 1, COL_INT, neighbors[LEFT], LTAG, cartcomm, &req[5]);
        }

        if (neighbors[RIGHT] >= 0 ){
            MPI_Isend(table_u + iz*sub_x*sub_y + sub_y + sub_table_dimention, 1, COL_INT, neighbors[RIGHT], LTAG, cartcomm,&req[6]);
            MPI_Irecv(table_u + iz*sub_x*sub_y + sub_y + sub_table_dimention + 1, 1, COL_INT , neighbors[RIGHT], RTAG, cartcomm,&req[7]);
        }

        /* Update inside table while the process wait for neighbor values */
        update_inside_table(sub_table_dimention - 2, table_u + iz*sub_x*sub_y, table_u + (1-iz)*sub_x*sub_y);

        /* Wait for neighbor values */
        if(neighbors[UP] >= 0){
            MPI_Wait(&req[0],&status[0]);
            MPI_Wait(&req[1],&status[1]);
        }
        if(neighbors[DOWN] >= 0){
            MPI_Wait(&req[2],&status[2]);
            MPI_Wait(&req[3],&status[3]);
        }
        if(neighbors[LEFT] >= 0){
            MPI_Wait(&req[4],&status[4]);
            MPI_Wait(&req[5],&status[5]);
        }
        if(neighbors[RIGHT] >= 0){
            MPI_Wait(&req[6],&status[6]);
            MPI_Wait(&req[7],&status[7]);
        }

        /* Update outside table with neighboor values */
        update_outside_table(sub_table_dimention, table_u + iz*sub_x*sub_y, table_u + (1-iz)*sub_x*sub_y);

        /* Next loop with have to deal with the other table */
        iz = 1 - iz;

#if (CONVERGENCE == 1)
        if (step % CONVERGENCE_n == 0) {
            task_convergence = check_convergence(1, task_X - 2, 1, task_Y - 2, task_Y, *u[iz], *u[1 - iz]);
            MPI_Barrier(MPI_COMM_CARTESIAN);
            MPI_Allreduce(&task_convergence, &reduced_convergence, 1, MPI_INT, MPI_LAND, MPI_COMM_CARTESIAN);
            if (reduced_convergence == 1) {
                break;
            }
        }
#endif
    }
    /* Stop timer, calculate duration, reduce */
    end_time = MPI_Wtime();
    task_time = end_time - start_time;
    MPI_Barrier(MPI_COMM_CARTESIAN);
    MPI_Reduce(&task_time, &reduced_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_CARTESIAN);
    if (rank == 0) {
        printf("Convergence: %d\n", CONVERGENCE);
        printf("u size: [%d][%d]\n", NXPROB, NYPROB);
        printf("tasks: %d\n", number_of_tasks);
        printf("Time elapsed: %f seconds\n", reduced_time);
    }
    /****************/

    /* Cleanup everything */
    free(temp[0]);
    free(temp[1]);
    free(u[0]);
    free(u[1]);
    free(u);
    MPI_Type_free(&MPI_row);
    MPI_Type_free(&MPI_column);
    MPI_Finalize();
}


/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update_calculation(int ix, int iy, int y, float *u1, float *u2)
{

    *(u2 + ix * y + iy) = *(u1 + ix * y + iy) +
        parms.cx * (*(u1 + (ix + 1) * y + iy) +
        *(u1 + (ix - 1) * y + iy) -
        2.0 * *(u1 + ix * y + iy)) +
        parms.cy * (*(u1 + ix * y + iy + 1) +
        *(u1 + ix * y + iy - 1) -
        2.0 * *(u1 + ix * y + iy));

}


void update_inside_table(int end, float *u1, float *u2)
{
    int i, j;
    for (i = 2; i <= end + 1; i++) {
        for (j = 2; j <= end + 1; j++) {
            update_calculation(i, j, end + 4, u1, u2);
        }
    }
}


void update_outside_table(int end, float *u1, float *u2)
{
    int i;
    for (i = 1; i <= end; i++) {
        update_calculation(1, i, end + 2, u1, u2);
        update_calculation(end, i, end + 2, u1, u2);
        update_calculation(i, 1, end + 2, u1, u2);
        update_calculation(i, end, end + 2, u1, u2);
    }
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, int id, float **u) {
int ix, iy;
int r = rand()%10;
for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     u[ix][iy] = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1))*(r+id);
}
