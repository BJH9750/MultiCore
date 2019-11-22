#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <pthread.h>
#include "typedef.h"
//enum class Flag: uint32_t{ insert = (1U << 30), query = (1U << 31), wait = (3U << 30)}; // 01XXXXX, 10XXXXX, 11XXXXX
//typedef unsigned int uint32_t;
#define QUEUE_SIZE 100000

const uint32_t INSERT = (1U << 30), QUERY = (1U << 31), WAIT = (3U << 30);
const uint32_t mask_v = (1U << 30) - 1;
const uint32_t mask_t = (3U << 30);

class TaskQueue
{
public: 

    const static uint32_t ERROR = (3U << 30);

    TaskQueue(int _size){
        head = 0;
        tail = 0;
        active = 0;
        arr = new uint32_t[QUEUE_SIZE];
    }

    ~TaskQueue(){
        delete[] arr;
    }

    void push(uint32_t _value){
        tail = (tail + 1) % QUEUE_SIZE;
        arr[tail] = _value;
    }

    uint32_t pop(){
        uint32_t hval = ERROR;
        head = (head + 1) % QUEUE_SIZE;
        hval = arr[head];
        return hval;
    }

    bool empty(){
        return (head == tail);
    }

    bool full(){
        return (tail + 1)  %  QUEUE_SIZE == head;
    }

    void print(){
        printf("head : %d tail : %d size : %d\n", head, tail, QUEUE_SIZE);
        for(int i = 0; i <= tail; ++i){
            printf("(%u %u) ", (arr[i] >> 30), (arr[i] << 2) >> 2);
        }
        printf("\n");
    }

protected:
    uint32_t head;
    uint32_t tail;
    uint32_t * arr;
    uint32_t active;
};