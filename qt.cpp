#include "taskqueue.h"
#include <time.h>
#include <string>

char *fn;
int thread_num;
int thread_cnt;
bool done = false;
pthread_mutex_t mtx;
pthread_cond_t cond;
pthread_barrier_t barrier;

void *thread_work(void * args){
    uint64_t task;
    TaskQueue *tasks = (TaskQueue *)args;
    printf("%lu created\n", pthread_self());

    pthread_barrier_wait(&barrier);

    while(!done){
        task &= 0;
        pthread_mutex_lock(&mtx);
        pthread_cond_wait(&cond, &mtx);
        if(!(tasks->empty())){
            task = tasks->pop();
            printf("%lu : %u %u\n", pthread_self(), (uint32_t)(task >> 62), (uint32_t) task);
        }
        pthread_mutex_unlock(&mtx);
    }
    ++thread_cnt;
    printf("%lu finished\n", pthread_self());
    pthread_exit(NULL);
}

void *thread_main(void *args){
    TaskQueue *tasks = (TaskQueue *)args;
    
    char action;
    uint64_t num;
    uint64_t task;

    FILE* fin = fopen(fn, "r");
    pthread_t *workers = (pthread_t *)malloc(sizeof(pthread_t) * thread_num);

    printf("main thread : %lu\n", pthread_self());

    for(int i = 0; i < thread_num; ++i){
        pthread_create(&workers[i], NULL, thread_work, tasks);
    }

    pthread_barrier_wait(&barrier);

    while (fscanf(fin, "%c %lu\n", &action, &num) == 2) {
        task &= 0;
        switch (action){
            case 'i':
                task = static_cast<uint64_t>(Flag::insert) | num;
                break;
            case 'q':
                task = static_cast<uint64_t>(Flag::query) | num;
                break;
            case 'w':
                task = static_cast<uint64_t>(Flag::wait) | num;
                break;
            default:
                printf("ERROR: Unrecognized action: '%c'\n", action);
                break;
        }
        tasks->push(task);
        pthread_cond_signal(&cond);     
    }
    
    while(!(tasks->empty())) pthread_cond_signal(&cond);

    done = true;
    while(thread_cnt != thread_num) pthread_cond_signal(&cond);

    for(int i = 0; i < thread_num; ++i){
        pthread_join(workers[i], NULL);
    }

    tasks->print();

    fclose(fin);
    pthread_mutex_destroy(&mtx);
    pthread_cond_destroy(&cond);
    pthread_barrier_destroy(&barrier);
    pthread_exit(NULL);
}

int main(int argc, char** argv){
    thread_num = atoi(argv[2]);
    fn = argv[1];

    struct timespec start, stop;
    pthread_t tmain;
    TaskQueue tasks(10);

    pthread_barrier_init(&barrier, NULL, thread_num + 1);
    pthread_mutex_init(&mtx, NULL);
    pthread_cond_init(&cond, NULL);

    pthread_create(&tmain, NULL, thread_main, &tasks);
    pthread_join(tmain, NULL);
    
}