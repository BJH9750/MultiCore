#include <iostream>
#include <sstream>
#include <pthread.h>

#define BILLION  1000000000L
#define NPAIRS  4

using namespace std;

pthread_rwlock_t skplock = PTHREAD_RWLOCK_INITIALIZER;

template<class K,class V,int MAXLEVEL>
class skiplist_node
{
public:

    skiplist_node():forwards{}, cnt(0), cur(0){}
 
    skiplist_node(K searchKey):forwards{}, cnt(1), cur(0){
        key[0] = searchKey;
    }
 
    skiplist_node(K searchKey, V val):forwards{}, cnt(1), cur(0){
	    key[0] = searchKey;
	    value[0] = val;
    }
 
    skiplist_node(K *_key, V* _value, int off, int _cnt){
        memcpy(key, _key + off, _cnt * sizeof(K));
        memcpy(value, _value + off, _cnt * sizeof(V));
        cnt = _cnt;
    }
    
    virtual ~skiplist_node(){}

    void insert(K k, V v){
	    for(int i = 0; i < cnt; ++i){
	        if(key[i] < k) 
		        continue;

            memmove(key + i + 1, key + i, (cnt - i) * sizeof(K));
            memmove(value  + i + 1, value + i, (cnt - i) * sizeof(V));

	        key[i] = k;
	        value[i] = v;
		    cnt++;
		    return;
	    }
        key[cnt] = k;
        value[cnt] = v;
        cnt++;
        return;
    }
    
    void print(){
        for(int i = 0; i < cnt; ++i){
            printf("(%d %d) ", key[i], value[i]);
        }
        printf("\n");
    }

    int cnt;
    int cur;
    K key[NPAIRS];
    V value[NPAIRS]; 
    skiplist_node<K,V,MAXLEVEL>* forwards[MAXLEVEL+1];
};
 
///////////////////////////////////////////////////////////////////////////////
 
template<class K, class V, int MAXLEVEL = 16>
class skiplist
{
public:
    typedef K KeyType;
    typedef V ValueType;
    typedef skiplist_node<K,V,MAXLEVEL> NodeType;
 
    skiplist(K minKey, K maxKey):m_pHeader(NULL),m_pTail(NULL),
                                max_curr_level(1),max_level(MAXLEVEL),
                                m_minKey(minKey),m_maxKey(maxKey)
    {
        m_pHeader = new NodeType(m_minKey);
        m_pTail = new NodeType(m_maxKey);
        for (int i = 1 ; i <= MAXLEVEL; ++i){
            m_pHeader->forwards[i] = m_pTail;
        }
    }
 
    virtual ~skiplist(){
        NodeType* currNode = m_pHeader->forwards[1];
        while (currNode != m_pTail) {
            NodeType* tempNode = currNode;
            currNode = currNode->forwards[1];
            delete tempNode;
        }
        delete m_pHeader;
        delete m_pTail;
    }
 
    void insert(K searchKey, V newValue){
        skiplist_node<K,V,MAXLEVEL>* update[MAXLEVEL];
        NodeType* currNode = m_pHeader;

        //pthread_rwlock_rdlock(&skplock);
        for(int level = max_curr_level; level >= 1; --level){
            while (currNode->forwards[level]->key[0] <= searchKey){
                currNode = currNode->forwards[level];
            }
            update[level] = currNode;
        }
        //pthread_rwlock_unlock(&skplock);
        //pthread_rwlock_wrlock(&skplock);

	    if(currNode->cnt < NPAIRS){
	        currNode->insert(searchKey, newValue);
	    }else{ // split
            int newlevel = randomLevel();
            if(newlevel > max_curr_level){
                for(int level = max_curr_level + 1; level <= newlevel; ++level){
                    update[level] = m_pHeader;
                }
                max_curr_level = newlevel;
            }
            int cnt = currNode->cnt;
            int mid = cnt / 2; 
            NodeType* newNode = new NodeType(currNode->key, currNode->value, mid, cnt - mid);
            
            currNode->cnt = mid;

            if(newNode->key[0] < searchKey){
                newNode->insert(searchKey, newValue);
            }else{
                currNode->insert(searchKey, newValue);
            }

            for(int lv = 1; lv <= max_curr_level; ++lv){
                newNode->forwards[lv] = update[lv]->forwards[lv];
                update[lv]->forwards[lv] = newNode; // make previous node point to new node
            }
        }
    }
 
    void erase(K searchKey){
	    /*
        skiplist_node<K,V,MAXLEVEL>* update[MAXLEVEL];
        NodeType* currNode = m_pHeader;
        for(int level=max_curr_level; level >=1; level--) {
            while ( currNode->forwards[level]->key < searchKey ) {
                currNode = currNode->forwards[level];
            }
            update[level] = currNode;
        }
        currNode = currNode->forwards[1];
        if ( currNode->key == searchKey ) {
            for ( int lv = 1; lv <= max_curr_level; lv++ ) {
                if ( update[lv]->forwards[lv] != currNode ) {
                    break;
                }
                update[lv]->forwards[lv] = currNode->forwards[lv];
            }
            delete currNode;
            // update the max level
            while ( max_curr_level > 1 && m_pHeader->forwards[max_curr_level] == NULL ) {
                max_curr_level--;
            }
        }
	*/
    }

    V find(K searchKey){
        NodeType* currNode = m_pHeader;
        for(int level = max_curr_level; level >= 1; --level){
            while(currNode->forwards[level]->key[0] <= searchKey){
                currNode = currNode->forwards[level];
            }
        }

        for(int i = 0; i < currNode->cnt; ++i){
            if(currNode->key[i] == searchKey){
                return currNode->value[i];
            }
        }
        return -1;
    }
 
    bool empty() const{
        return ( m_pHeader->forwards[1] == m_pTail );
    }
 
    std::string printList(){
	    int i=0;
        cout << "list" << endl;
        std::stringstream sstr;
        NodeType* currNode = m_pHeader; //->forwards[1];
        while ( currNode != m_pTail ) {
            sstr << "(" ;
            for(int i=0;i<currNode->cnt;i++){
                sstr << currNode->key[i] << "," ;
                }
            sstr << ")";
            currNode = currNode->forwards[1];
            i++;
            if(i>200) break;
        }
        return sstr.str();
    }
 
    const int max_level;
 
protected:
    double uniformRandom(){
        return rand() / double(RAND_MAX);
    }
 
    int randomLevel(){
        int level = 1;
        double p = 0.5;
        while(uniformRandom() < p && level < MAXLEVEL){
            level++;
        }
        return level;
    }

    K m_minKey;
    K m_maxKey;
    int max_curr_level;
    skiplist_node<K,V,MAXLEVEL>* m_pHeader;
    skiplist_node<K,V,MAXLEVEL>* m_pTail;
};
 
///////////////////////////////////////////////////////////////////////////////
 
