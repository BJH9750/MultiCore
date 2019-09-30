#include "taskqueue.h"
#include "skiplist2.h"
#include <time.h>
#include <string>

char *fn;
int thread_num;
int thread_cnt;
bool done = false;
pthread_mutex_t mtx;
pthread_cond_t cond;
pthread_cond_t full_cv;
pthread_barrier_t barrier;
pthread_rwlock_t rwlock = PTHREAD_RWLOCK_INITIALIZER;
skiplist<uint32_t, uint32_t> list(0, 1000000);
int verbose = 0;
long sum, odd;

void *thread_work(void * args){
    uint32_t task;
    TaskQueue *tasks = (TaskQueue *)args;
    if(verbose) printf("%lu created\n", pthread_self());

    pthread_barrier_wait(&barrier);

    while(!done){
        task &= 0;
        pthread_mutex_lock(&mtx);
        while(tasks->empty() && !done)
            pthread_cond_wait(&cond, &mtx);
        
        task = tasks->pop();
        pthread_cond_signal(&full_cv);
        pthread_mutex_unlock(&mtx);

        uint32_t flag = task & mask_t;
        task &= mask_v;
        if(verbose) printf("%lu %u \n",pthread_self(), task);
        switch (flag){
            case INSERT:
                pthread_rwlock_wrlock(&rwlock);
                if(verbose) printf("%lu insert %u \n",pthread_self(), task);
                list.insert(task, task);
                pthread_rwlock_unlock(&rwlock);
                break;
            case QUERY:
                pthread_rwlock_rdlock(&rwlock);
                if(verbose) printf("%lu query %u \n",pthread_self(), task);
                if(list.find(task) != task)
                    if(verbose)
		                printf("ERROR: Not Found: %u\n", task);
                pthread_rwlock_unlock(&rwlock);
                break;
            case WAIT:
                if(verbose) printf("%lu wait %u \n",pthread_self(), task);
                struct timeval tv;
                tv.tv_sec = task / 1000;
                tv.tv_usec = task % 1000;
                select(0, NULL, NULL, NULL, &tv);
                break;
            default:
                break;
        }
    }
    pthread_mutex_lock(&mtx);
    ++thread_cnt;
    pthread_mutex_unlock(&mtx);
    
    if(verbose) printf("%lu finished\n", pthread_self());
    pthread_exit(NULL);
}

void *thread_main(void *args){
    TaskQueue *tasks = (TaskQueue *)args;
    
    char action;
    uint32_t num;
    uint32_t task;

    FILE* fin = fopen(fn, "r");
    pthread_t *workers = (pthread_t *)malloc(sizeof(pthread_t) * thread_num);

    if(verbose) printf("main thread : %lu\n", pthread_self());

    for(int i = 0; i < thread_num; ++i){
        pthread_create(&workers[i], NULL, thread_work, tasks);
    }

    pthread_barrier_wait(&barrier);

    while (fscanf(fin, "%c %u\n", &action, &num) == 2) {
        sum += num;
        if (num % 2 == 1) odd++;
        task &= 0;
        switch (action){
            case 'i':
                task = INSERT | num;
                break;
            case 'q':
                task = QUERY | num;
                break;
            case 'w':
                task = WAIT | num;
                break;
            default:
                if(verbose) printf("ERROR: Unrecognized action: '%c'\n", action);
                break;
        }

        pthread_mutex_lock(&mtx);
        while(tasks->full())
            pthread_cond_wait(&full_cv, &mtx);
        tasks->push(task);
        pthread_cond_signal(&cond);   
        pthread_mutex_unlock(&mtx); 
    }
    
    while(!(tasks->empty())) pthread_cond_signal(&cond);

    done = true;

    while(thread_cnt != thread_num) pthread_cond_signal(&cond);

    for(int i = 0; i < thread_num; ++i){
        pthread_join(workers[i], NULL);
    }

    fclose(fin);
    pthread_exit(NULL);
}

int main(int argc, char** argv){
    thread_num = atoi(argv[2]);
    fn = argv[1];
    if(argv[3] != NULL) verbose = atoi(argv[3]);

    struct timespec start, stop;
    pthread_t tmain;
    TaskQueue tasks(thread_num);

    pthread_barrier_init(&barrier, NULL, thread_num + 1);
    pthread_mutex_init(&mtx, NULL);
    pthread_cond_init(&cond, NULL);
    pthread_cond_init(&full_cv, NULL);

    clock_gettime( CLOCK_REALTIME, &start);

    pthread_create(&tmain, NULL, thread_main, &tasks);
    pthread_join(tmain, NULL);

    clock_gettime( CLOCK_REALTIME, &stop);

    cout << list.printList() << endl;
    cout << sum << " " << odd << endl;
    cout << "Elapsed time: " << (stop.tv_sec - start.tv_sec) + ((double) (stop.tv_nsec - start.tv_nsec))/BILLION << " sec" << endl;

    if(verbose) tasks.print();
    pthread_barrier_destroy(&barrier);
    pthread_mutex_destroy(&mtx);
    pthread_cond_destroy(&cond);
    pthread_cond_destroy(&full_cv);

    return (EXIT_SUCCESS);
}