#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cstdint>
#include <unistd.h>
#include <pthread.h>

enum class Flag: uint32_t{ insert = (1U << 30), query = (1U << 31), wait = (3U << 30)}; // 01XXXXX, 10XXXXX, 11XXXXX
const uint32_t mask = (1U << 30) - 1;

class TaskQueue
{
public: 

    const uint32_t ERROR = (3U << 30);

    TaskQueue(int _size){
        size =  _size  * 2 + 1;
        head = 0;
        tail = 0;
        arr = new uint32_t[size];
        pthread_rwlock_init(&rwlock, NULL);
    }

    ~TaskQueue(){
        delete[] arr;
        pthread_rwlock_destroy(&rwlock);
    }

    void push(uint32_t _value){
        pthread_rwlock_wrlock(&rwlock);
        if((tail + 1)  %  size == head){
            size *= 2;
            uint32_t *tmp = new uint32_t[size];
            memcpy(tmp, arr, (size / 2) * sizeof(uint32_t));
            delete[] arr;
            arr = tmp;
        }
        tail = (tail + 1) % size;
        arr[tail] = _value;
        pthread_rwlock_unlock(&rwlock);
    }

    uint32_t pop(){
        uint32_t hval = ERROR;
        pthread_rwlock_wrlock(&rwlock);
        if(head != tail){
            head = (head + 1) % size;
            hval = arr[head];
        }
        pthread_rwlock_unlock(&rwlock);
        return hval;
    }

    uint32_t front(){
        return head;
    }

    uint32_t back(){
        return tail;
    }

    bool empty(){
        bool isEmpty;
        pthread_rwlock_rdlock(&rwlock);
        isEmpty = (head == tail);
        pthread_rwlock_unlock(&rwlock);
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
    pthread_rwlock_t rwlock;
};