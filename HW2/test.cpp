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
    #pragma omp parallel sections
    {
            #pragma omp section
            {
                x(1);
            }
            #pragma omp section
            {
                x(2);
            }
            #pragma omp section
            {
                x(3);
            }
            #pragma omp section
            {
                x(4);
            }
            #pragma omp section
            {
                x(5);
            }
            
        }
        #pragma omp barrier
    printf("asd\n");
}