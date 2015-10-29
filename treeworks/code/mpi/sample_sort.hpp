/***
 *  Project: TreeWorks
 *
 *  File: sample_sort.hpp
 *  Created: Mar 18, 2009
 *
 *  Author: Abhinav Sarje <abhinav.sarje@gmail.com>
 *  Adapted from Benjamin G. Jackson's ParallelLib.h
 */

#ifndef SAMPLE_SORT_HPP
#define SAMPLE_SORT_HPP

#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <string.h>

namespace twlib {

template <typename T>
void quick_sort(T* data, int num) {

	// special cases
	if(num <= 1) return;

	if(num == 2) {
		if(data[1] < data[0]) {
			T temp = data[0];
			data[0] = data[1];
			data[1] = temp;
		}
		return;
	}

	if(num == 3) {
		if(data[1] < data[0]) {
			T temp = data[0];
			data[0] = data[1];
			data[1] = temp;
		}
		if(data[2] < data[1]) {
			T temp = data[1];
			data[1] = data[2];
			data[2] = temp;

			if(data[1] < data[0]) {
				T temp = data[0];
				data[0] = data[1];
				data[1] = temp;
			}
		}
		return;
	}

	// general case

	// pick a random pivot
	srand((unsigned)time(0));
	int random_index = rand() % num;
	T pivot = data[random_index];

	// we will divide the array into 3 parts, that part less than the pivot, that part equal to
	// the pivot and that part greater than the pivot, using the following three indexes
	int index1 = 0;
	int index2 = num - 1;
	int index3;

	T temp;
	// divide the array into Part 1 and Part 2
	// Part 1:  < pivot
	// Part 2:  >= pivot
	do {
		while(index1 < num && data[index1] < pivot) index1++;
		while(index2 >= 0 && !(data[index2] < pivot)) index2--;

		if(index1 < index2) {
			temp = data[index2];
			data[index2] = data[index1];
			data[index1] = temp;
		}
	} while(index1 < index2);

	index3 = index1;
	index2 = num-1;

	// index1 now holds the start index of part 2 (and the size of part1)

	// divide Part 2 into parts 2a and 2b
	// Part 2a: == pivot
	// Part 2b: > pivot
	do {
		while(index3 < num && !(pivot < data[index3])) index3++;
		while(index2 >= index1 && pivot < data[index2]) index2--;

		if(index3 < index2) {
			temp = data[index3];
			data[index3] = data[index2];
			data[index2] = temp;
		}
	} while(index3 < index2);

	// index3 now holds the start location of part 2b
	// recursively sort parts 1 and 2b

	quick_sort(data, index1);
	quick_sort(data + index3, num - index3);
} // quick_sort()


static void set_displacements(int* counts, int* displacements, int proc) {
	displacements[0] = 0;
	for(int i = 1; i < proc; ++i)
		displacements[i] = displacements[i-1] + counts[i-1];
} // set_displacements()


template <typename T>
void all_gather(T* &data, int& n, MPI_Comm comm) {
	int total = 0, proc = 0, rank = 0;

	MPI_Comm_size(comm, &proc);
	MPI_Comm_rank(comm, &rank);

	MPI_Allreduce(&n, &total, 1, MPI_INT, MPI_SUM, comm);

	if(total == 0) return;

	T* all_data = new T[total];
	memset(all_data, 0, total*sizeof(T));

	n *= sizeof(T);

	int* counts = new int[proc];
	int* displacements = new int[proc];

	memset(counts, 0, sizeof(int)*proc);
	memset(displacements, 0, sizeof(int)*proc);

	MPI_Allgather(&n, 1, MPI_INT, counts, 1, MPI_INT, comm);

	set_displacements(counts, displacements, proc);

	MPI_Allgatherv(data, n, MPI_BYTE, all_data, counts, displacements,  MPI_BYTE, comm);

	delete[] counts;
	delete[] displacements;
	delete[] data;

	data = all_data;
	n = total;
} // all_gather()


// all_gather without throwing out the input
template <typename T>
void all_gather(T* &data, int n, T* &recv_data, int& total, MPI_Comm comm) {
	int proc;

	MPI_Comm_size(comm, &proc);

	MPI_Allreduce(&n, &total, 1, MPI_INT, MPI_SUM, comm);

	if(total == 0) return;

	T* all_data = new T[total];
	memset(all_data, 0, total*sizeof(T));

	n *= sizeof(T);

	int* counts = new int[proc];
	int* displacements = new int[proc];

	memset(counts, 0, sizeof(int)*proc);
	memset(displacements, 0, sizeof(int)*proc);

	MPI_Allgather(&n, 1, MPI_INT, counts, 1, MPI_INT, comm);

	set_displacements(counts, displacements, proc);

	MPI_Allgatherv(data, n, MPI_BYTE, all_data, counts, displacements,  MPI_BYTE, comm);

	delete[] counts;
	delete[] displacements;

	recv_data = all_data;
} // all_gather()


template <typename T>
void all_to_all(int* &counts, T* &data, int& n, MPI_Comm comm) {
	int proc;
	MPI_Comm_size(comm, &proc);

	int* receive_counts = new int[proc];
	int* receive_displacements = new int[proc];
	int* displacements = new int[proc];

	memset(receive_counts, 0, sizeof(int)*proc);
	memset(receive_displacements, 0, sizeof(int)*proc);
	memset(displacements, 0, sizeof(int)*proc);

	T* result = 0;
	int result_size;

	for(int i = 0; i < proc; ++i) counts[i] *= sizeof(T);

	set_displacements(counts, displacements, proc);

	// gather data info to prepare for all to all
	MPI_Alltoall(counts, 1, MPI_INT, receive_counts, 1, MPI_INT, comm);

	// create receive info
	set_displacements(receive_counts, receive_displacements, proc);
	result_size = receive_displacements[proc-1] + receive_counts[proc-1];
	result_size /= sizeof(T);
	result = new T[result_size];

	memset(result, 0, sizeof(T) * result_size);
	MPI_Alltoallv(data, counts, displacements, MPI_BYTE, result, receive_counts,
			receive_displacements, MPI_BYTE, comm);

	for(int i = 0; i < proc; ++i) receive_counts[i] /= sizeof(T);

	delete[] counts;
	delete[] data;

	counts = receive_counts;
	data = result;
	n = result_size;

	delete[] displacements;
	delete[] receive_displacements;
} // all_to_all()


template <typename T>
void distribute(T* &data, int& n, T* boundaries,MPI_Comm comm) {
	int proc;
	MPI_Comm_size(comm, &proc);

	int* counts = new int[proc];  // holds the number of elements to be sent to each processor
	memset(counts, 0, sizeof(int) * proc);

	int i = 0; //i = data array index, j = processor index
	for (int j = 0; j < proc; j++) {
		while (i < n && !(boundaries[j] < data[i])) { // data[i] <= boundaries[j]
			counts[j]++;
			i++;
		}
	}

	// some elements might be greater than the last boundary element.  This will be sent to the 
	// last processor.
	while (i < n) {
		counts[proc-1]++;
		i++;
	}

	// perform the all to all, data is replaced by the elements recieved by this processor
	// size is also updated.     
	all_to_all(counts, data, n, comm);

	delete[] counts;
} // distribute()




// a overloading where the send buffers are not replaced with recv, they are kept separate
template <typename T>
void all_to_all(int* &send_counts, T* &send_data,
		int* &recv_counts, T* &recv_data, int& n, MPI_Comm comm) {
	int proc;
	MPI_Comm_size(comm, &proc);

	int* receive_counts = new int[proc];
	int* temp_send_counts = new int[proc];
	int* receive_displacements = new int[proc];
	int* displacements = new int[proc];

	memset(receive_counts, 0, sizeof(int)*proc);
	memset(receive_displacements, 0, sizeof(int)*proc);
	memset(displacements, 0, sizeof(int)*proc);

	T* result = 0;
	int result_size;

	for(int i = 0; i < proc; ++i) temp_send_counts[i] = send_counts[i] * sizeof(T);

	set_displacements(temp_send_counts, displacements, proc);

	// gather data info to prepare for all to all
	MPI_Alltoall(temp_send_counts, 1, MPI_INT, receive_counts, 1, MPI_INT, comm);

	// create receive info
	set_displacements(receive_counts, receive_displacements, proc);
	result_size = receive_displacements[proc-1] + receive_counts[proc-1];
	result_size /= sizeof(T);
	result = new T[result_size];

	memset(result, 0, sizeof(T) * result_size);
	MPI_Alltoallv(send_data, temp_send_counts, displacements, MPI_BYTE, result, receive_counts,
			receive_displacements, MPI_BYTE, comm);

	for(int i = 0; i < proc; ++i) receive_counts[i] /= sizeof(T);

	recv_counts = receive_counts;
	recv_data = result;
	n = result_size;

	delete[] temp_send_counts;
	delete[] displacements;
	delete[] receive_displacements;
} // all_to_all()


template<typename T>
bool sample_sort(T* &data, int &n, const MPI_Comm comm) {
	// T has comparison operators implemented
	
	int rank = 0, proc = 0;
	
	MPI_Comm_size(comm, &proc);
	MPI_Comm_rank(comm, &rank);
	uint64_t n_global = 0;

	MPI_Allreduce(&n, &n_global, 1, MPI_LONG_LONG, MPI_SUM, comm);
	if(n_global == 0) return false;

	// sort the local data
	quick_sort<T>(data, n);

	// choose local samples
	double run = (double) n / sqrt(proc);
	int offset = (int) (run * (double) rank / (double) proc);
	int sample_size = 0;

	for(double k = offset; (int)k < n; k += run) ++ sample_size;

	T* samples = 0;
    samples = new T[sample_size];

	int i = 0;
	double k;
	for(k = offset, i = 0; (int)k < n && i < sample_size; k += run, ++ i)
		samples[i] = data[(int)k];

	// collect all samples
	all_gather(samples, sample_size, comm);

	// sort the samples
	quick_sort<T>(samples, sample_size);

	// choose the splitters and bucket the local data
	T *boundaries = new T[proc];

	double loc = run = (double)(sample_size - 1) / (double)proc;
	for (i = 0; i < proc; ++ i, loc += run)
		boundaries[i] = samples[(int)loc];

	// distribute the data among all processors
	distribute(data, n, boundaries, comm);

	// perform final local sort on the received data
	quick_sort<T>(data, n);

	delete[] boundaries;
	delete[] samples;

	return true;
} // sample_sort()

} // namespace twlib

#endif // SAMPLE_SORT_HPP
