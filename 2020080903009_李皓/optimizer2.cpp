#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define ODD_TO_INDEX(odd) (((odd)-3)/2)     //Odd-to-index
#define INDEX_TO_ODD(index) (2*(index)+3)   //对序号进行转换
#define BLOCK_LOW(id,p,n) (id * (long long)(n) / p)
#define BLOCK_HIGH(id,p,n) ((id+1) * (long long)(n) / p -1)


int main(int argc, char* argv[])
{
    int    count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    long long    first;        /* Index of first multiple */
    int    global_count=0; /* Global prime count */
    long long    high_value;   /* Highest value on this proc */
    long long    i;
    int    id=0;           /* Process ID number */
    long long    index;        /* Index of current prime */
    long long    low_value;    /* Lowest value on this proc */
    char*   marked;       /* Portion of 2,...,'n' */
    char*   marked_prime; // 用于记录串行计算的素数结果
    long long    n;            /* Sieving from 2, ..., 'n' */
    long long    m;            // 将实际的编号 ODD 转化为序号 INDEX 
    int    p;            /* Number of processes */
    int    proc0_size;   /* Size of proc 0's subarray */
    long long    prime;        /* Current prime */
    long long    size;         /* Elements in 'marked' */
    long long    prime_size;   // to record prime size.
    
    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    if (argc != 2) {
        printf("argc is %d", argc);
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    //{//Just for debugger.if not,comment it.
    //    int temp;
    //    if (id == 0)
    //    {
    //        std::cin >> temp;
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD); // All threads will wait here until you
    //    //give thread 0 an input
    //}

    n = std::stoll(argv[1]);
    m = ODD_TO_INDEX(n) +1;// +1 确保了不出问题
    // printf("m: %lld ,n: %lld \n",m,n);


    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */
    // 对数组的长度进行修改
    low_value = INDEX_TO_ODD(BLOCK_LOW(id,p,m));
    high_value = INDEX_TO_ODD(BLOCK_HIGH(id,p,m));
    size = (high_value - low_value)/2 + 1;
    prime_size = ODD_TO_INDEX(sqrt(n))+1;
    // printf("low_val:%lld, high:%lld, size:%lld \n",low_value,high_value,size);


    /* Bail out if all the primes used for sieving are
       not all held by process 0 */

    proc0_size = m / p;

    if (INDEX_TO_ODD(proc0_size-1) < (int)sqrt((double)n)) {
        if (!id) printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */
    
    marked = (char*)malloc(size);
    marked_prime = (char*)malloc(prime_size);
    // printf(" primesize: %lld  \n",prime_size);
    if (marked == NULL||marked_prime== NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    for (i = 0; i < prime_size;i++) marked_prime[i] = 0;
    for (i = 0; i < size; i++) marked[i] = 0;

    index = 0;
    prime = 3;
    do{
        for (i = ODD_TO_INDEX(prime * prime); i < prime_size; i += prime)
            marked_prime[i] = 1;
        while (marked_prime[++index]);
        prime = INDEX_TO_ODD(index);
    }while(prime*prime<=sqrt(n));

    index = 0;
    prime = 3;
    do {
        if (prime * prime > low_value)
            first = ODD_TO_INDEX( prime * prime) -ODD_TO_INDEX( low_value);
            // 从prime*prime 开始找，因为（prime-1）之前的因数已经找过了
            // 对于分段后的数据而言，是先前一段的内容都被找到了
        else {
            if (!(low_value % prime)) 
                first = 0;
            // 如果开始的数就是因子的整数倍，那么就从这个开始找
            else {
                first = prime - (low_value % prime);
                if (!((low_value + first)%2)) 
                    first += prime;           
                    // 跳过偶数
                first /= 2;
                // 映射回index
            }
            // 如果不是，则为prime - (low_value % prime)
            // 假设 low_value 是 10 ，prime 为 3
            // first 就是 3-1= 2
            // 对应的就是mark[2]的位置，也就是 10，11，[12] 处
        }
        for (i = first; i < size; i += prime) marked[i] = 1;
        
            while (marked_prime[++index]);
            // 找到下一个素数因子
            prime = INDEX_TO_ODD(index);
            // 从映射回序数
            // malloc 中的内容是从 0 开始索引，计算素数从2 开始。
        
    } while (prime * prime <= n);
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i]) count++;

    MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
        0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();


    /* Print the results */

    if (!id) {
        printf("There are %d primes less than or equal to %lld\n",
            global_count+1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}
