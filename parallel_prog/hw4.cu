//Udacity HW 4
//Radix Sorting

#include "reference_calc.cpp"
#include "utils.h"

#include <stdio.h>

/* Red Eye Removal
   ===============

   For this assignment we are implementing red eye removal.  This is
   accomplished by first creating a score for every pixel that tells us how
   likely it is to be a red eye pixel.  We have already done this for you - you
   are receiving the scores and need to sort them in ascending order so that we
   know which pixels to alter to remove the red eye.

   Note: ascending order == smallest to largest

   Each score is associated with a position, when you sort the scores, you must
   also move the positions accordingly.

   Implementing Parallel Radix Sort with CUDA
   ==========================================

   The basic idea is to construct a histogram on each pass of how many of each
   "digit" there are.   Then we scan this histogram so that we know where to put
   the output of each digit.  For example, the first 1 must come after all the
   0s so we have to know how many 0s there are to be able to start moving 1s
   into the correct position.

   1) Histogram of the number of occurrences of each digit
   2) Exclusive Prefix Sum of Histogram
   3) Determine relative offset of each digit
        For example [0 0 1 1 0 0 1]
                ->  [0 1 0 1 2 3 2]
   4) Combine the results of steps 2 & 3 to determine the final
      output location for each element and move it there

   LSB Radix sort is an out-of-place sort and you will need to ping-pong values
   between the input and output buffers we have provided.  Make sure the final
   sorted results end up in the output buffer!  Hint: You may need to do a copy
   at the end.

 */


__global__
void performParallelPrefixSum(unsigned int *d_input, int numElem)
{
    int tId = threadIdx.x;
    extern __shared__ int shdata[];
    int *sdataIn = (int *)shdata;
    int *sdataOut = (int *)(&shdata[blockDim.x]);
    int step;

    if (tId >= numElem)
        return;
    sdataIn[tId] = d_input[tId];
    sdataOut[tId] = sdataIn[tId];
    __syncthreads();
    for (unsigned int n = 0; ((pow(2.0, (double)n)) < blockDim.x); n++)
    {
        step = pow(2.0, (double)n);
        if (tId >= step)
        {
            sdataOut[tId] = sdataIn[tId] + sdataIn[tId-step];
        }
        __syncthreads();
        sdataIn[tId] = sdataOut[tId];
        __syncthreads();
    }
    d_input[tId] = sdataIn[tId];
}


__global__
void getCounts(unsigned int *d_inputVals, unsigned int * d_gPSum,
               unsigned int *d_cnt0s, unsigned int *d_cnt1s, unsigned int *d_cnt2s, 
               unsigned int *d_cnt3s,
               int numElem, int shift)
{

    int tId = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int mask, val, index;
    __shared__ int cnt0, cnt1, cnt2, cnt3; //operating on 2-bit. Hence 4 combinations

/*    extern __shared__ int s_data[];
    int *s_presum0 = (int *)s_data;
    int *s_presum1 = (int *)(&s_data[blockDim.x]);
    int *s_presum2 = (int *)(&s_data[2*blockDim.x]);
    int *s_presum3 = (int *)(&s_data[3*blockDim.x]);
*/
    if (tId >= numElem)
        return;

    if (threadIdx.x == 0)
      cnt0 = cnt1 = cnt2 = cnt3 = 0;

  /*  s_presum0[threadIdx.x] = 0;
    s_presum1[threadIdx.x] = 0;
    s_presum2[threadIdx.x] = 0;
    s_presum3[threadIdx.x] = 0;
    */
    __syncthreads();
    mask = 3 << shift;
    val = (mask & d_inputVals[tId]) >> shift;
    switch (val)
    {
        case 0:
            atomicAdd(&cnt0, 1);

        break;

        case 1:
            atomicAdd(&cnt1, 1);

        break;

        case 2:
            atomicAdd(&cnt2, 1);

        break;

        case 3:
            atomicAdd(&cnt3, 1);

    }
    __syncthreads();
    if (threadIdx.x == 0)
    {
        atomicAdd(&d_gPSum[0], cnt0);
        atomicAdd(&d_gPSum[1], cnt1);
        atomicAdd(&d_gPSum[2], cnt2);
        atomicAdd(&d_gPSum[3], cnt3);

        d_cnt0s[blockIdx.x] = cnt0;
        d_cnt1s[blockIdx.x] = cnt1;
        d_cnt2s[blockIdx.x] = cnt2;
        d_cnt3s[blockIdx.x] = cnt3;
    }
}

__global__
void scatter(unsigned int *d_inputVals, unsigned int *d_inputPos, unsigned int *d_globalPrefixSum,
             unsigned int *d_blockwise0s, unsigned int *d_blockwise1s, unsigned int *d_blockwise2s,
             unsigned int *d_blockwise3s,
             unsigned int *d_ovals, unsigned int *d_opos,
             int numElem, int shift)
{
    int tId = blockIdx.x * blockDim.x + threadIdx.x;
    int mask, val, new_index;
    __shared__ int cnt0, cnt1, cnt2, cnt3; //operating on 2-bit. Hence 4 combinations

    if (tId >= numElem)
        return;

    if (threadIdx.x > 4)
      return;

    if (threadIdx.x == 0)
      cnt0 = cnt1 = cnt2 = cnt3 = 0;

    __syncthreads();

    mask = 3 << shift;
    for (unsigned int i = 0; i < blockDim.x; i++) {
      unsigned int id = blockIdx.x * blockDim.x + i;
      val = (mask & d_inputVals[id]) >> shift;

      if (threadIdx.x != val) continue;

      if (val == 0)
        new_index = d_blockwise0s[blockIdx.x] + cnt0++;
      else if (val == 1)
        new_index = d_blockwise1s[blockIdx.x] + cnt1++;
      else if (val == 2)
        new_index = d_blockwise2s[blockIdx.x] + cnt2++;
      else if (val == 3)
        new_index = d_blockwise3s[blockIdx.x] + cnt3++;
      else
        printf("Something wrong\n");
    }

    new_index += d_globalPrefixSum[val];
    d_ovals[new_index] = d_inputVals[tId];
    d_opos[new_index] = d_inputPos[tId];

}

__global__
void mycopy(unsigned int *d_inputVals, unsigned int *d_inputPos,
            unsigned int *d_outputVals, unsigned int *d_outputPos, int numElem)
{
    int tId = blockIdx.x * blockDim.x + threadIdx.x;

    if (tId >= numElem)
        return;
    d_outputVals[tId] = d_inputVals[tId];
    d_outputPos[tId] = d_inputPos[tId];
}


void performSerialPrefixSum(int *h_globalPrefixSum, int numBins)
{
    int sum = 0;
    int *temp_buf = (int *)malloc(sizeof(int) * numBins);
    memset(temp_buf, 0, sizeof(int) * numBins);

    for (unsigned int i = 1; i < numBins; i++)
    {
        temp_buf[i] = sum;
        sum = sum + h_globalPrefixSum[i];
    }
    memcpy(h_globalPrefixSum, temp_buf, sizeof(int) * numBins);
    free(temp_buf);
}

void getCountsSerial(unsigned int *h_inputVals, unsigned int *cmp, int numblocks,
                     int numElems, int shift)
{
    int mask, val;
    int cnt0, cnt1, cnt2, cnt3;

    cnt0 = cnt1 = cnt2 = cnt3 = 0;
    mask = 3 << shift;
    for (unsigned int i = 0; i < numElems; i++) {
      val = (mask & h_inputVals[i]) >> shift;
      if (val == 0)
        cnt0++;
      else if (val == 1)
        cnt1++;
      else if (val == 2)
        cnt2++;
      else if (val == 3) {
        cnt3++;
      }
      else
        printf("Something wrong\n");
          
     // printf("%d %d \n", val, h_check[i]);
    }
    printf("shift:%d \n", shift);
    printf("Serial: %d %d %d %d\n", cnt0, cnt1, cnt2, cnt3);
    printf("Parallel: %d %d %d %d\n", cmp[0], cmp[1], cmp[2], cmp[3]);
         

}

void your_sort(unsigned int* const d_inputVals,
               unsigned int* const d_inputPos,
               unsigned int* const d_outputVals,
               unsigned int* const d_outputPos,
               const size_t numElems)
{
     unsigned int numBits = 2;
     unsigned int numBins = 1 << numBits;
     unsigned int *d_blockwise0s, *d_blockwise1s, *d_blockwise2s, *d_blockwise3s;
     unsigned int *d_globalPrefixSum;
     unsigned int *h_globalPrefixSum;
     unsigned int *d_ipVals, *d_ipPos, *d_opVals, *d_opPos;


     dim3 blockSize(1024, 1, 1);
     dim3 gridSize((numElems + blockSize.x - 1)/blockSize.x, 1, 1);
     printf("Elements: %d BlockDim: %d NumBlocks: %d\n", numElems, blockSize.x,
                                      (numElems + blockSize.x - 1)/blockSize.x);

     checkCudaErrors(cudaMalloc((void **)&d_blockwise0s, sizeof(unsigned int) * gridSize.x ));
     checkCudaErrors(cudaMalloc((void **)&d_blockwise1s, sizeof(unsigned int) * gridSize.x ));
     checkCudaErrors(cudaMalloc((void **)&d_blockwise2s, sizeof(unsigned int) * gridSize.x ));
     checkCudaErrors(cudaMalloc((void **)&d_blockwise3s, sizeof(unsigned int) * gridSize.x ));
     checkCudaErrors(cudaMalloc((void **)&d_globalPrefixSum, sizeof(unsigned int) * numBins));
     checkCudaErrors(cudaMalloc((void **)&d_ipVals, sizeof(unsigned int) * numElems));
     checkCudaErrors(cudaMalloc((void **)&d_ipPos, sizeof(unsigned int) * numElems));
     checkCudaErrors(cudaMalloc((void **)&d_opVals, sizeof(unsigned int) * numElems));
     checkCudaErrors(cudaMalloc((void **)&d_opPos, sizeof(unsigned int) * numElems));

     h_globalPrefixSum = (unsigned int *)malloc(sizeof(unsigned int) * numBins);
     unsigned int *h_inputVals = (unsigned int*)malloc(sizeof(unsigned int)*numElems);;
     cudaMemcpy(h_inputVals, d_inputVals, sizeof(unsigned int)*numElems, cudaMemcpyDeviceToHost);

      
     mycopy<<<gridSize, blockSize>>>(d_inputVals, d_inputPos, d_ipVals, d_ipPos, numElems);
     for (unsigned int i = 2; i < 3/*8 * sizeof(unsigned int)*/; i = i + numBits)
     {
         checkCudaErrors(cudaMemset(d_globalPrefixSum, 0,  sizeof(unsigned int) * numBins));
         checkCudaErrors(cudaMemset(d_blockwise0s, 0,  sizeof(unsigned int) * gridSize.x));
         checkCudaErrors(cudaMemset(d_blockwise1s, 0,  sizeof(unsigned int) * gridSize.x));
         checkCudaErrors(cudaMemset(d_blockwise2s, 0,  sizeof(unsigned int) * gridSize.x));
         checkCudaErrors(cudaMemset(d_blockwise3s, 0,  sizeof(unsigned int) * gridSize.x));

         //checkCudaErrors(cudaMemset(d_opVals, 0,  sizeof(unsigned int) * numElems));
         //checkCudaErrors(cudaMemset(d_opPos, 0,  sizeof(unsigned int) * numElems));
         getCounts<<<gridSize, blockSize>>>(d_ipVals, d_globalPrefixSum, d_blockwise0s,
                                            d_blockwise1s, d_blockwise2s, d_blockwise3s,
                                            numElems, i);
         //cudaMemcpy(h_check, d_check, sizeof(unsigned int) * numElems, cudaMemcpyDeviceToHost);
         performParallelPrefixSum<<<1, gridSize, 2 * sizeof(unsigned int) * gridSize.x>>>
                                  (d_blockwise0s, 10000);
         performParallelPrefixSum<<<1, gridSize, 2 * sizeof(unsigned int) * gridSize.x>>>
                                  (d_blockwise1s, 10000);
         performParallelPrefixSum<<<1, gridSize, 2 * sizeof(unsigned int) * gridSize.x>>>
                                  (d_blockwise2s, 10000);
         performParallelPrefixSum<<<1, gridSize, 2 * sizeof(unsigned int) * gridSize.x>>>
                                  (d_blockwise3s, 10000);

         cudaMemcpy(h_globalPrefixSum, d_globalPrefixSum, sizeof(unsigned int) * numBins,
                    cudaMemcpyDeviceToHost);
         getCountsSerial(h_inputVals, h_globalPrefixSum, 216, numElems, i);
         unsigned int *h_temp = (unsigned int *)malloc(sizeof(unsigned int) * numBins);
         memset(h_temp, 0, sizeof(unsigned int) * numBins);
         for (int j = 1; j < numBins; j++)
            h_temp[j] = h_temp[j-1] + h_globalPrefixSum[j-1];

         cudaMemcpy(d_globalPrefixSum, h_temp, sizeof(unsigned int) * numBins, cudaMemcpyHostToDevice);

         /*scatter<<<gridSize, blockSize>>>(d_ipVals, d_ipPos, d_globalPrefixSum, d_blockwise0s,
                                          d_blockwise1s, d_blockwise2s, d_blockwise3s,
                                          d_opVals, d_opPos,
                                          numElems, i);*/
        /*unsigned int *d_temp;
        d_temp = d_opVals;
        d_opVals = d_ipVals;
        d_ipVals = d_temp;

        d_temp = d_opPos;
        d_opPos = d_ipPos;
        d_ipPos = d_temp;*/
   //     mycopy<<<gridSize, blockSize>>>(d_opVals, d_opPos, d_ipVals, d_ipPos, numElems);
     }
     mycopy<<<gridSize, blockSize>>>(d_ipVals, d_ipPos, d_outputVals, d_outputPos, numElems);
     checkCudaErrors(cudaFree(d_globalPrefixSum));
     checkCudaErrors(cudaFree(d_blockwise0s));
     checkCudaErrors(cudaFree(d_blockwise1s));
     checkCudaErrors(cudaFree(d_blockwise2s));
     checkCudaErrors(cudaFree(d_blockwise3s));
     checkCudaErrors(cudaFree(d_ipVals));
     checkCudaErrors(cudaFree(d_ipPos));
     //free(h_binsPrefixSum);
}

