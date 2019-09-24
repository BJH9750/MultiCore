#include "taskqueue.h"
#include <time.h>
#include <string>

char *fn;
int thread_num;
bool done = false;
pthread_mutex_t mtx;
pthread_cond_t cond;
pthread_rwlock_t rwlock = PTHREAD_RWLOCK_INITIALIZER;
pthread_barrier_t barrier;

void *thread_work(void * args){
    uint64_t task;
    TaskQueue *tasks = (TaskQueue *)args;
    printf("%llu created\n", pthread_self());
    //pthread_barrier_wait(&barrier);
    while(!done){
        task &= 0;
        pthread_mutex_lock(&mtx);
        pthread_cond_wait(&cond, &mtx);
        if(!(tasks->empty())){
            task = tasks->pop();
            printf("%llu : %u %u\n", pthread_self(), task >> 62, (uint32_t) task);
        }
        pthread_mutex_unlock(&mtx);
        
    }
    printf("%llu finished\n", pthread_self());
    pthread_exit(NULL);
}

void *thread_main(void *args){
    TaskQueue *tasks = (TaskQueue *)args;
    
    char action;
    uint64_t num;
    uint64_t task;

    FILE* fin = fopen(fn, "r");
    pthread_t *workers = (pthread_t *)malloc(sizeof(pthread_t) * thread_num);

    printf("main thread : %llu\n", pthread_self());

    for(int i = 0; i < thread_num; ++i){
        pthread_create(&workers[i], NULL, thread_work, tasks);
    }

    while (fscanf(fin, "%c %lu\n", &action, &num) == 2) {
        task &= 0;
        printf("read : %c %ld\n", action, num);
        if (action == 'i') {            // insert
            task = static_cast<uint64_t>(Flag::insert) | num;
        }else if (action == 'q') {      // qeury
            task = static_cast<uint64_t>(Flag::query) | num;
        }else if (action == 'w') {     // wait
            task = static_cast<uint64_t>(Flag::wait) | num; 
        }else {
            printf("ERROR: Unrecognized action: '%c'\n", action);
        }
        tasks->push(task);
        pthread_cond_signal(&cond);     
    }
    
    while(!(tasks->empty())) pthread_cond_signal(&cond);

    done = true;
    pthread_cond_broadcast(&cond);

    for(int i = 0; i < thread_num; ++i){
        
        pthread_join(workers[i], NULL);
    }

    fclose(fin);
    tasks->print();
    pthread_mutex_destroy(&mtx);
    pthread_cond_destroy(&cond);
    pthread_barrier_destroy(&barrier);
    pthread_exit(NULL);
}



int main(int argc, char** argv){
    thread_num = atoi(argv[2]);
    fn = argv[1];
    int v = atoi(argv[3]);
    struct timespec start, stop;
    pthread_t tmain;
    int status;
    TaskQueue tasks(10, v);

    //clock_gettime( CLOCK_REALTIME, &start);

    pthread_barrier_init(&barrier, NULL, thread_num);
    pthread_mutex_init(&mtx, NULL);
    pthread_cond_init(&cond, NULL);

    pthread_create(&tmain, NULL, thread_main, &tasks);
    pthread_join(tmain, NULL);
    
    //clock_gettime( CLOCK_REALTIME, &stop);

    //cout << "Elapsed time: " << (stop.tv_sec - start.tv_sec) + ((double) (stop.tv_nsec - start.tv_nsec))/1000000000 << " sec" << endl;
    
}