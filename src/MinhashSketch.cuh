#ifndef MINHASHSKETCH_H
#define MINHASHSKETCH_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <cstdint>
#include <cmath>
#include <queue>
#include <thread>
#include <unordered_set>
#include <cuda_runtime.h>
#include "SpookyV2_d.h"
#include <cub/cub.cuh>

typedef vector<vector<uint64>> signature;

__device__ int base2int(char base);
__device__ uint64 getHashValue (uint64 *x, uint64 b, int k);
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__device__ void BlockGetHashValues (const int k, char *dna_d, int thread_offset, uint64 *thread_dataR, int hash_b);
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__device__ void BlockRemoveDupe (const int m, int *thread_dataS, int *thread_dataD,uint64 *thread_dataR, uint64 *sketch);
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__device__ void BlockGetBack (const int m, int threadID, uint64 *thread_dataR, uint64 *sketch);
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void getBlockSketch (const int k, const int m, char *dna_d, uint64 *input_d, int numElem_dna, int numElem_list, uint64 hash_b);
template<int BLOCK_THREADS, int ITEMS_PER_THREAD>
__device__ void BlockGetRank(int m, uint64 *thread_data64, int *thread_dataR, int *thread_dataD, uint64 *shared_B);
template<int BLOCKS_NUM, int BLOCK_THREADS, int ITEMS_PER_THREAD>
__global__ void getAllSketch (const int m, uint64 *input_d, int numElem_list);
void rMerge(const int m, uint64 *sketch_h, uint64 *output_h);
signature genSig(const int k, const int m, const int t, char *dnaList, int length, uint64 *hashes_b);

#endif