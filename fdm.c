/*****************************************************************
/* Author: Pablo Goitia <pablo.goitia@alumnos.unican.es>
/* Project: Finite-Difference Approximation to the Heat Equation
/* Date: Mar-2025
/*
/* Usage: ./fdm [deltaH]
/*****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

unsigned int mdf_heat(double ***__restrict__ u0,
					  double ***__restrict__ u1,
					  const unsigned int npX,
					  const unsigned int npY,
					  const unsigned int npZ,
					  const double deltaH,
					  const double deltaT,
					  const double inErr,
					  const double boundaries,
					  MPI_Comm compute_comm,
					  const int myrank,
					  const int size,
					  const unsigned int first_point,
					  const unsigned int last_point)
{

	register double alpha = deltaT / (deltaH * deltaH);
	int continued = 1;
	register unsigned int steps = 0;

	// Loop indexes
	unsigned int i, j, k, i_internal;
	unsigned int points_per_slice = last_point - first_point + 1;

	// Boundaries
	register double left;
	register double right;
	register double up;
	register double down;
	register double top;
	register double bottom;

	double ***ptr;
	double err, maxErr;

	// MPI comms handling variables
	MPI_Status status;
	MPI_Request request_left, request_right;

	// MPI Pack and Unpack variables
	int position = 0;

	// Buffer for MPI communication
	double *buffer = (double *)malloc(npY * npZ * sizeof(double));
	double t_block_start, t_block_end, total_time = 0.0f;

	#pragma omp parallel \
	shared(u0, u1, inErr, boundaries, compute_comm, myrank, size, first_point, last_point, alpha, continued, steps, points_per_slice, ptr, status, request_left, request_right, buffer) \
	private(i, j, k, i_internal, left, right, up, down, top, bottom, err, maxErr, position)
	{
		while (continued)
		{
			#pragma omp single
			{
				steps++;
				t_block_start = MPI_Wtime();
			}

			err = 0.0f;
			maxErr = 0.0f;

			#pragma omp for schedule(static) reduction(& : continued)
			for (i = first_point; i <= last_point; i++)
			{
				for (j = 0; j < npY; j++)
				{
					for (k = 0; k < npZ; k++)
					{
						i_internal = i - first_point + 1;
						left = boundaries;
						right = boundaries;
						up = boundaries;
						down = boundaries;
						top = boundaries;
						bottom = boundaries;

						if ((i > 0) && (i < (npX - 1)))
						{
							left = u0[i_internal - 1][j][k];
							right = u0[i_internal + 1][j][k];
						}
						else if (i == 0)
							right = u0[i_internal + 1][j][k];
						else
							left = u0[i_internal - 1][j][k];

						if ((j > 0) && (j < (npY - 1)))
						{
							up = u0[i_internal][j - 1][k];
							down = u0[i_internal][j + 1][k];
						}
						else if (j == 0)
							down = u0[i_internal][j + 1][k];
						else
							up = u0[i_internal][j - 1][k];

						if ((k > 0) && (k < (npZ - 1)))
						{
							top = u0[i_internal][j][k - 1];
							bottom = u0[i_internal][j][k + 1];
						}
						else if (k == 0)
							bottom = u0[i_internal][j][k + 1];
						else
							top = u0[i_internal][j][k - 1];

						u1[i_internal][j][k] = alpha * (top + bottom + up + down + left + right - (6.0f * u0[i_internal][j][k])) + u0[i_internal][j][k];
						err = fabs(u1[i_internal][j][k] - boundaries);
						if (err > inErr)
							maxErr = err;
						else
							continued = 0;
					}
				}
			}
			
			#pragma omp single
			{
				t_block_end = MPI_Wtime();
				total_time += t_block_end - t_block_start;

				// Use MPI_LAND to check if any process has set 'continued' to 0
				MPI_Allreduce(&continued, &continued, 1, MPI_INT, MPI_LAND, compute_comm);
			}
			
			if (continued == 0)
				break;

			#pragma omp single
			{
				// Swap the pointers (the next instant of time is now the current time)
				ptr = u0;
				u0 = u1;
				u1 = ptr;

				// Exchange first and last point values with neighbour MPI processes
				if (myrank > 0)
				{
					position = 0;

					// Send the content of the first point of the slice to the left neighbour
					for (i = 0; i < npY; i++)
						MPI_Pack(&u0[1][i][0], npZ, MPI_DOUBLE, buffer, npY * npZ * sizeof(double), &position, compute_comm);

					MPI_Isend(buffer, npY * npZ, MPI_DOUBLE, myrank - 1, steps, compute_comm, &request_left);

					position = 0;

					// Receive the content of the first point of the slice from the left neighbour
					MPI_Recv(buffer, npY * npZ, MPI_DOUBLE, myrank - 1, steps, compute_comm, &status);
					for (i = 0; i < npY; i++)
						MPI_Unpack(buffer, npY * npZ * sizeof(double), &position, &u0[0][i][0], npZ, MPI_DOUBLE, compute_comm);
				}

				if (myrank < (size - 1))
				{
					position = 0;

					// Send the content of the last point of the slice to the right neighbour
					for (i = 0; i < npY; i++)
						MPI_Pack(&u0[points_per_slice][i][0], npZ, MPI_DOUBLE, buffer, npY * npZ * sizeof(double), &position, compute_comm);

					MPI_Isend(buffer, npY * npZ, MPI_DOUBLE, myrank + 1, steps, compute_comm, &request_right);

					position = 0;

					// Receive the content of the last point of the slice from the right neighbour
					MPI_Recv(buffer, npY * npZ, MPI_DOUBLE, myrank + 1, steps, compute_comm, &status);
					for (i = 0; i < npY; i++)
						MPI_Unpack(buffer, npY * npZ * sizeof(double), &position, &u0[points_per_slice + 1][i][0], npZ, MPI_DOUBLE, compute_comm);
				}

				// Check if the messages have been sent and received
				if (myrank > 0)
					MPI_Wait(&request_left, &status);

				if (myrank < (size - 1))
					MPI_Wait(&request_right, &status);
			}
		}
	}

	printf("Process %d: Computation block time: %f seconds\n", myrank, total_time);

	return steps;
}

int main(int ac, char **av)
{
	// Data structures
	double ***u0;
	double ***u1;

	// Loop indexes
	unsigned int i, j;

	// Simulation parameters
	double deltaT = 0.01;
	double deltaH = 0.005f;
	double sizeX = 1.0f;
	double sizeY = 1.0f;
	double sizeZ = 1.0f;

	double *aux0, *aux1;

	unsigned int npX, npY, npZ;

	// MPI variables
	int my_global_rank;
	MPI_Status status;

	// Subcommunicator for the MPI processes that will compute the points
	MPI_Comm compute_comm;
	int myrank, size;
	int is_computing_proc = 0;

	// Variables for the MPI processes
	unsigned int points_per_slice, remaining_points;
	unsigned int first_point, last_point;

	// Variables for results
	unsigned int steps, max_steps;
	double start_time, end_time, execution_time, max_time;

	// Start MPI communication
	MPI_Init(&ac, &av);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_global_rank);

	// Input processing (common for all MPI processes)
	if (ac > 1)
	{
		deltaH = atof(av[1]);
		if (deltaH <= 0.0f || deltaH >= 1.0f)
		{
			if (my_global_rank == 0)
				fprintf(stderr, "Error: deltaH must be a value in the range (0.0, 1.0)\n");
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
	}

	// Calculate the number of points in each axis
	npX = (unsigned int)(sizeX / deltaH);
	npY = (unsigned int)(sizeY / deltaH);
	npZ = (unsigned int)(sizeZ / deltaH);

	// Discriminate the MPI processes that will compute the points and the remaining ones
	is_computing_proc = (my_global_rank < npX) ? 1 : 0;
	MPI_Comm_split(MPI_COMM_WORLD, is_computing_proc, my_global_rank, &compute_comm);
	MPI_Comm_rank(compute_comm, &myrank);
	MPI_Comm_size(compute_comm, &size);

	if (is_computing_proc)
	{
		/** The three-dimensional space will be divided into a number of slices across the X axis
		 * matching the number of MPI processes, and which are made up of a set of points. */
		points_per_slice = npX / size; // Provisional value

		/** Some slices might have more points than others if npX cannot be exactly divided by
		 * the number of MPI processes, so we will implement the simplest approach for a load
		 * balancing algorithm:
		 *  	Distribute the remaining points among the first MPI processes.
		 *  	Assuming equal workload for each point and equal hardware resources, this would
		 *  	be the fairest algorithm.
		 *
		 *  Obviously, those conditions are not frequently met, so this algorithm is not
		 *  the best one, but it will be enough for this first implementation.
		 */
		remaining_points = npX % size;

		// Calculate the definitive number of points for each MPI process
		if (myrank < remaining_points)
			points_per_slice++;

		/** Calculate the first and last point in the slice for each MPI process,
		 * considering that some processes will have one more point. */
		first_point = 0;
		for (i = 0; i < myrank; i++)
			first_point += (i < npX % size) ? (points_per_slice + 1) : points_per_slice;
		last_point = first_point + points_per_slice - 1;

		// Processes print the number of points in each axis
		printf("I am process %d of %d. points_per_slice=%d of %d. first_point=%d, last_point=%d\n",
			   myrank + 1, size, points_per_slice, npX, first_point, last_point);

		// Allocating memory for the tri-dimensional space
		/** Memory allocation for X axis.
		 * We will allocate memory for:
		 * - The points in each slice, that is, the value of points_per_slice.
		 * - Two more points in the X axis, because the algorithm will need to access the
		 *   left and right neighbour points.
		 */
		u0 = (double ***)malloc((points_per_slice + 2) * sizeof(double **));
		u1 = (double ***)malloc((points_per_slice + 2) * sizeof(double **));

		// Memory allocation for Y axis
		for (i = 0; i < (points_per_slice + 2); i++)
		{
			u0[i] = (double **)malloc(npY * sizeof(double *));
			u1[i] = (double **)malloc(npY * sizeof(double *));
		}

		// Memory allocation for Z axis
		for (i = 0; i < (points_per_slice + 2); i++)
		{
			for (j = 0; j < npY; j++)
			{
				aux0 = (double *)malloc(npZ * sizeof(double));
				aux1 = (double *)malloc(npZ * sizeof(double));
				// Initial condition: zero in all points
				memset(aux0, 0x01, npZ * sizeof(double));
				memset(aux1, 0x02, npZ * sizeof(double));
				u0[i][j] = aux0;
				u1[i][j] = aux1;
			}
		}

		start_time = MPI_Wtime();

		// Each MPI process will compute its assigned points with the finite difference method
		steps = mdf_heat(u0, u1, npX, npY, npZ, deltaH, deltaT, 1e-15, 100.0f, compute_comm, myrank, size, first_point, last_point);

		// Collect the number of steps from all MPI processes and get the maximum value
		MPI_Reduce(&steps, &max_steps, 1, MPI_UNSIGNED, MPI_MAX, 0, compute_comm);

		end_time = MPI_Wtime();
		execution_time = end_time - start_time;

		// Get max execution time across all processes
		MPI_Reduce(&execution_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, compute_comm);

		// Print the results
		if (myrank == 0)
			printf("Done! in %u steps\nExecution time: %f seconds\n", max_steps, max_time);

		// Free memory
		for (i = 0; i < (points_per_slice + 2); i++)
		{
			for (j = 0; j < npY; j++)
			{
				free(u0[i][j]);
				free(u1[i][j]);
			}
		}

		for (i = 0; i < (points_per_slice + 2); i++)
		{
			free(u0[i]);
			free(u1[i]);
		}

		free(u0);
		free(u1);
	}

	// Free the subcommunicators
	// Processes that do not compute points will also terminate their own subcommunicator
	MPI_Comm_free(&compute_comm);

	MPI_Finalize();

	return EXIT_SUCCESS;
}
