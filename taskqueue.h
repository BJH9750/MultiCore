#include <cstdlib>
#include <cstdio>
#include <climits>
#include <unistd.h>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <iostream>
#include <pthread.h>


using namespace std;

enum class Flag: uint64_t{ insert = (1ULL << 62), query = (1ULL << 63), wait = (3ULL << 62)};

class TaskQueue
{
public: 

    uint64_t ERROR = (3ULL << 62);

    TaskQueue(int _size, int _v){
        size =  _size  * 2 + 1;
        head = 0;
        tail = 0;
        arr = new uint64_t[size];
        pthread_mutex_init(&mtx, NULL);

        v = _v;
    }

    void push(uint64_t _value){
        //pthread_mutex_lock(&mtx);
        if((tail + 1)  %  size == head){
            size *= 2;
            uint64_t *tmp = new uint64_t[size];
            memcpy(tmp, arr, (size / 2) * sizeof(uint64_t));
            delete[] arr;
            arr = tmp;
        }
        tail = (tail + 1) % size;
        arr[tail] = _value;
        if(v) printf("%lu %d inserted\n", _value >> 62, (int)_value);
        //pthread_mutex_unlock(&mtx);
    }

    uint64_t pop(){
        uint64_t hval = ERROR;
        //pthread_mutex_lock(&mtx);
        if(head != tail){
            head = (head + 1) % size;
            hval = arr[head];
        }
        if(v){
            if(hval == ERROR) printf("pop error\n");
            else printf("%lu %d popped\n", hval >> 62, (int)hval);
        }
        
        //pthread_mutex_unlock(&mtx);
        return hval;
    }

    uint32_t front(){
        return head;
    }

    uint32_t back(){
        return tail;
    }

    bool empty(){
        return head == tail;
    }

    void print(){
        printf("head : %d tail : %d size : %d\n", head, tail, size);
        for(int i = 0; i < size; ++i){
            cout << "(" << (arr[i] >> 62) << " " << (uint32_t)arr[i] << ") ";
        }
        cout << endl;
    }

    ~TaskQueue(){
        delete[] arr;
        pthread_mutex_destroy(&mtx);
    }

protected:
    uint32_t head;
    uint32_t tail;
    uint32_t size;
    uint64_t * arr;
    pthread_mutex_t mtx;
    bool v;
};