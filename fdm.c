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
#include <unistd.h>

#define TYPES_OF_NODES 	2
#define MAX_FREQ 		4.8
#define MIN_FREQ 		3.4

typedef struct _node {
	char root_name[10];
	int num_nodes;
	int cores_per_node;
	float frequency;
} node;

static node nodes[TYPES_OF_NODES] = {
	{"n16-8", 4, 6, 4.8},	// Fast nodes
	{"n16-9", 4, 4, 3.4}	// Slower nodes
};

unsigned int mdf_heat(double *__restrict__ u0,
					  double *__restrict__ u1,
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
	unsigned int i, j, k, idx;
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

	// Buffer for MPI communication
	double t_block_start, t_block_end, total_time = 0.0f;

	#pragma omp parallel \
	shared(u0, u1, inErr, boundaries, compute_comm, myrank, size, first_point, last_point, alpha, continued, steps, points_per_slice, ptr, status, request_left, request_right) \
	private(i, j, k, idx, left, right, up, down, top, bottom, err, maxErr)
	{
		while (continued)
		{
			#pragma omp single nowait
			{
				steps++;
				t_block_start = MPI_Wtime();
			}

			err = 0.0f;
			maxErr = 0.0f;

			#pragma omp for schedule(dynamic) reduction(& : continued)
			for (i = first_point; i <= last_point; i++)
			{
				for (j = 0; j < npY; j++)
				{
					for (k = 0; k < npZ; k++)
					{
						left = boundaries;
						right = boundaries;
						up = boundaries;
						down = boundaries;
						top = boundaries;
						bottom = boundaries;

						idx = (i - first_point + 1) * npY * npZ + j * npZ + k;

						if ((i > 0) && (i < (npX - 1)))
						{
							left = u0[idx - npY * npZ];
							right = u0[idx + npY * npZ]; 
						}
						else if (i == 0)
							right = u0[idx + npY * npZ];
						else
							left = u0[idx - npY * npZ];

						if ((j > 0) && (j < (npY - 1)))
						{
							up = u0[idx - npZ];
							down = u0[idx + npZ];
						}
						else if (j == 0)
							down = u0[idx + npZ];
						else
							up = u0[idx - npZ];

						if ((k > 0) && (k < (npZ - 1)))
						{
							top = u0[idx - 1];
							bottom = u0[idx + 1];
						}
						else if (k == 0)
							bottom = u0[idx + 1];
						else
							top = u0[idx - 1];

						u1[idx] = alpha * (top + bottom + up + down + left + right - (6.0f * u0[idx])) + u0[idx];

						err = fabs(u1[idx] - boundaries);
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
				// Exchange first and last point values with neighbour MPI processes
				if (myrank > 0)
				{
					// Send the content of the first point of the slice to the left neighbour
					MPI_Isend(&u1[npZ * npY], npY * npZ, MPI_DOUBLE, myrank - 1, steps, compute_comm, &request_left);

					// Receive the content of the first point of the slice from the left neighbour
					MPI_Irecv(&u1[0], npY * npZ, MPI_DOUBLE, myrank - 1, steps, compute_comm, &request_left);

				}

				if (myrank < (size - 1))
				{
					// Send the content of the last point of the slice to the right neighbour
					MPI_Isend(&u1[npZ * npY * points_per_slice], npY * npZ, MPI_DOUBLE, myrank + 1, steps, compute_comm, &request_right);

					// Receive the content of the last point of the slice from the right neighbour
					MPI_Irecv(&u1[npZ * npY * (points_per_slice + 1)], npY * npZ, MPI_DOUBLE, myrank + 1, steps, compute_comm, &request_right);
				}

				// Check if the messages have been received
				if (myrank > 0)
					MPI_Wait(&request_left, &status);

				if (myrank < (size - 1))
					MPI_Wait(&request_right, &status);

				// Swap the pointers (the next instant of time is now the current time)
				double *temp = u0;
				u0 = u1;
				u1 = temp;
			}
		}
	}

	printf("Process %d: Computation block time: %f seconds\n", myrank, total_time);

	return steps;
}

int main(int ac, char **av)
{
	// Data structures
	double *u0;
	double *u1;

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
	node mynode;

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

	// Determine the node where the process is running
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	// Search the node in the list of nodes
	for (i = 0; i < TYPES_OF_NODES; i++)
	{
		if (strncmp(processor_name, nodes[i].root_name, strlen(nodes[i].root_name)) == 0)
		{
			mynode = nodes[i];
			break;
		}
	}

	// Calculate the number of points in each axis
	npX = (unsigned int)(sizeX / deltaH);
	npY = (unsigned int)(sizeY / deltaH);
	npZ = (unsigned int)(sizeZ / deltaH);

	// The X axis is divided into a fair number of slices depending on the power of the cores
	double node_capacity = 0;
	double W = 0;	// Total capacity of the system (total frequency)

	for (int i = 0; i < TYPES_OF_NODES; i++)
		W += nodes[i].cores_per_node * nodes[i].frequency;

	double node_capacity = mynode.cores_per_node * mynode.frequency;

	// Workload fractions per node
	double fraction = node_capacity / W;

	// Point assignment per process
	points_per_slice = floor(npX * fraction);

	// Discriminate the MPI processes that will compute the points and the remaining ones
	is_computing_proc = (points_per_slice > 0) ? 1 : 0;
	MPI_Comm_split(MPI_COMM_WORLD, is_computing_proc, my_global_rank, &compute_comm);
	MPI_Comm_rank(compute_comm, &myrank);
	MPI_Comm_size(compute_comm, &size);

	if (is_computing_proc)
	{
		/* Calculate the first point for each MPI process based on their processing power.
		 * We need to sum the points assigned to previous processes. */
		int* global_points = malloc(size * sizeof(int));

		// Get the number of points assigned to each MPI process
		MPI_Allgather(&points_per_slice, 1, MPI_INT, global_points, 1, MPI_INT, compute_comm);

		// Calculate the first point for this MPI process
		first_point = 0;
		for (i = 0; i < myrank; i++)
			first_point += global_points[i];

		// Calculate the last point for this MPI process
		last_point = first_point + points_per_slice - 1;

		// Processes print the number of points to compute
		printf("I am process %d of %d. points_per_slice=%d of %d. first_point=%d, last_point=%d\n",
			   myrank + 1, size, points_per_slice, npX, first_point, last_point);

		/** Allocating memory for the tri-dimensional space, including:
		 * - The points in each slice, that is, the value of points_per_slice.
		 * - Two more points in the X axis, because the algorithm will need to access the
		 *   left and right neighbour points.
		 */
		u0 = malloc((points_per_slice + 2) * npY * npZ * sizeof(double *));
		u1 = malloc((points_per_slice + 2) * npY * npZ * sizeof(double *));

		// Initial condition: zero in all points
		memset(u0, 0x01, (points_per_slice + 2) * npY * npZ * sizeof(double *));
		memset(u1, 0x02, (points_per_slice + 2) * npY * npZ * sizeof(double *));

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

			if (myrank == 0)
			{
				printf("Frequency of cores in node: %d\n", get_frequency_of_cores_in_node());
			}

		// Free memory
		free(u0);
		free(u1);
	}

	// Free the subcommunicators
	// Processes that do not compute points will also terminate their own subcommunicator
	MPI_Comm_free(&compute_comm);

	MPI_Finalize();

	return EXIT_SUCCESS;
}
int get_frequency_of_cores_in_node() {
	FILE *fp;
	int freq = 0;
	
	fp = fopen("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq", "r");
	if (fp != NULL) {
		fscanf(fp, "%d", &freq);
		fclose(fp);
		return freq / 1000; // Convert from KHz to MHz
	}
	return -1; // Return -1 if unable to read frequency
}