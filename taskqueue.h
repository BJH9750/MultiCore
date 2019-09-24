#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cstdint>
#include <unistd.h>
#include <pthread.h>

enum class Flag: uint64_t{ insert = (1ULL << 62), query = (1ULL << 63), wait = (3ULL << 62)};

class TaskQueue
{
public: 

    const uint64_t ERROR = (3ULL << 62);

    TaskQueue(int _size){
        size =  _size  * 2 + 1;
        head = 0;
        tail = 0;
        arr = new uint64_t[size];
        pthread_rwlock_init(&rwlock, NULL);
    }

    ~TaskQueue(){
        delete[] arr;
        pthread_rwlock_destroy(&rwlock);
    }

    void push(uint64_t _value){
        pthread_rwlock_wrlock(&rwlock);
        if((tail + 1)  %  size == head){
            size *= 2;
            uint64_t *tmp = new uint64_t[size];
            memcpy(tmp, arr, (size / 2) * sizeof(uint64_t));
            delete[] arr;
            arr = tmp;
        }
        tail = (tail + 1) % size;
        arr[tail] = _value;
        pthread_rwlock_unlock(&rwlock);
    }

    uint64_t pop(){
        uint64_t hval = ERROR;
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
            printf("(%u %u) ", (uint32_t)(arr[i] >> 62), (uint32_t)arr[i]);
        }
        printf("\n");
    }

protected:
    uint32_t head;
    uint32_t tail;
    uint32_t size;
    uint64_t * arr;
    pthread_rwlock_t rwlock;
};