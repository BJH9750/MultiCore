#include <cstdio>
#include <cstdlib>
#include <omp.h>

int c[5];

void x(int k){
    c[k%5] = k;
    printf("%d\n", omp_get_thread_num());
}

int main(){
    omp_set_num_threads(5);
    #pragma omp parallel
    {
        x(omp_get_thread_num());
    }
}