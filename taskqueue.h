#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cstdint>
#include <unistd.h>
#include <pthread.h>

//enum class Flag: uint32_t{ insert = (1U << 30), query = (1U << 31), wait = (3U << 30)}; // 01XXXXX, 10XXXXX, 11XXXXX
const uint32_t INSERT = (1U << 30), QUERY = (1U << 31), WAIT = (3U << 30);
const uint32_t mask_v = (1U << 30) - 1;
const uint32_t mask_t = (3U << 30);

class TaskQueue
{
public: 

    const static uint32_t ERROR = (3U << 30);

    TaskQueue(int _size){
        size =  _size  * 2;
        head = 0;
        tail = 0;
        arr = new uint32_t[size];
        pthread_mutexattr_init(&tqlockattr);
        pthread_mutexattr_settype(&tqlockattr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&tqlock, &tqlockattr);
    }

    ~TaskQueue(){
        delete[] arr;
        pthread_mutex_destroy(&tqlock);
        pthread_mutexattr_destroy(&tqlockattr);
    }

    void push(uint32_t _value){
        //pthread_mutex_lock(&tqlock);
        if((tail + 1)  %  size == head){
            size *= 2;
            uint32_t *tmp = new uint32_t[size];
            memmove(tmp, arr, (size / 2) * sizeof(uint32_t));
            delete[] arr;
            arr = tmp;
        }
        tail = (tail + 1) % size;
        arr[tail] = _value;
        //pthread_mutex_unlock(&tqlock);
    }

    uint32_t pop(){
        //pthread_mutex_lock(&tqlock);
        uint32_t hval = ERROR;
        if(head != tail){
            head = (head + 1) % size;
            hval = arr[head];
        }
        //pthread_mutex_unlock(&tqlock);
        return hval;
    }

    bool empty(){
        //pthread_mutex_lock(&tqlock);
        bool isEmpty = (head == tail);
        //pthread_mutex_unlock(&tqlock);
        return isEmpty;
    }

    void print(){
        printf("head : %d tail : %d size : %d\n", head, tail, size);
        for(int i = 0; i < size; ++i){
            printf("(%u %u) ", (arr[i] >> 30), (arr[i] << 2) >> 2);
        }
        printf("\n");
    }

protected:
    uint32_t head;
    uint32_t tail;
    uint32_t size;
    uint32_t * arr;
    pthread_mutex_t tqlock;
    pthread_mutexattr_t tqlockattr;
};