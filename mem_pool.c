/*
 * Created by Ivo Georgiev on 2/9/16.
 */

#include <stdlib.h>
#include <assert.h>
#include <stdio.h> // for perror()

#include "mem_pool.h"

/*************/
/*           */
/* Constants */
/*           */
/*************/
#define     MEM_FILL_FACTOR                  0.75;
#define     MEM_EXPAND_FACTOR                2;

static const unsigned   MEM_POOL_STORE_INIT_CAPACITY    = 20;
static const float      MEM_POOL_STORE_FILL_FACTOR      = MEM_FILL_FACTOR;
static const unsigned   MEM_POOL_STORE_EXPAND_FACTOR    = MEM_EXPAND_FACTOR;

static const unsigned   MEM_NODE_HEAP_INIT_CAPACITY     = 40;
static const float      MEM_NODE_HEAP_FILL_FACTOR       = MEM_FILL_FACTOR;
static const unsigned   MEM_NODE_HEAP_EXPAND_FACTOR     = MEM_EXPAND_FACTOR;

static const unsigned   MEM_GAP_IX_INIT_CAPACITY        = 40;
static const float      MEM_GAP_IX_FILL_FACTOR          = MEM_FILL_FACTOR;
static const unsigned   MEM_GAP_IX_EXPAND_FACTOR        = MEM_EXPAND_FACTOR;



/*********************/
/*                   */
/* Type declarations */
/*                   */
/*********************/
typedef struct _node {
    alloc_t alloc_record;
    unsigned used;
    unsigned allocated;
    struct _node *next, *prev; // doubly-linked list for gap deletion
} node_t, *node_pt;

typedef struct _gap {
    size_t size;
    node_pt node;
} gap_t, *gap_pt;

typedef struct _pool_mgr {
    pool_t pool;
    node_pt node_heap;
    unsigned total_nodes;
    unsigned used_nodes;
    gap_pt gap_ix;
    unsigned gap_ix_capacity;
} pool_mgr_t, *pool_mgr_pt;



/***************************/
/*                         */
/* Static global variables */
/*                         */
/***************************/
static pool_mgr_pt *pool_store = NULL; // an array of pointers, only expand
static unsigned pool_store_size = 0;
static unsigned pool_store_capacity = 0;



/********************************************/
/*                                          */
/* Forward declarations of static functions */
/*                                          */
/********************************************/
static alloc_status _mem_resize_pool_store();
static alloc_status _mem_resize_node_heap(pool_mgr_pt pool_mgr);
static alloc_status _mem_resize_gap_ix(pool_mgr_pt pool_mgr);
static alloc_status
        _mem_add_to_gap_ix(pool_mgr_pt pool_mgr,
                           size_t size,
                           node_pt node);
static alloc_status
        _mem_remove_from_gap_ix(pool_mgr_pt pool_mgr,
                                size_t size,
                                node_pt node);
static alloc_status _mem_sort_gap_ix(pool_mgr_pt pool_mgr);



/****************************************/
/*                                      */
/* Definitions of user-facing functions */
/*                                      */
/****************************************/
alloc_status mem_init() {
    // ensure that it's called only once until mem_free
    // allocate the pool store with initial capacity
    // note: holds pointers only, other functions to allocate/deallocate
    if(pool_store == NULL) {
        pool_store = (pool_mgr_pt *) malloc(sizeof(pool_t) * MEM_POOL_STORE_INIT_CAPACITY);
        for(size_t t = 0 ;  t < MEM_POOL_STORE_INIT_CAPACITY; t++){
            pool_store[t] = NULL;
        }
        pool_store_capacity =  MEM_POOL_STORE_INIT_CAPACITY;
        return ALLOC_OK;
    }
    return ALLOC_CALLED_AGAIN;

}

alloc_status mem_free() {
    // ensure that it's called only once for each mem_init
    // make sure all pool managers have been deallocated
    // can free the pool store array
    // update static variables

    if(pool_store != NULL) {
        for(size_t t = 0 ; t < pool_store_size; t++){
            if(&pool_store[t] != NULL) {
                mem_pool_close(&pool_store[t]->pool);
            }
        }
        pool_store = NULL;
        pool_store_size = 0;

        return ALLOC_OK;
    }
    return ALLOC_CALLED_AGAIN;

}

pool_pt mem_pool_open(size_t size, alloc_policy policy) {
    // make sure there the pool store is allocated
    // expand the pool store, if necessary
    // allocate a new mem pool mgr
    // check success, on error return null
    // allocate a new memory pool
    // check success, on error deallocate mgr and return null
    // allocate a new node heap
    // check success, on error deallocate mgr/pool and return null
    // allocate a new gap index
    // check success, on error deallocate mgr/pool/heap and return null
    // assign all the pointers and update meta data:
    //   initialize top node of node heap
    //   initialize top node of gap index
    //   initialize pool mgr
    //   link pool mgr to pool store
    // return the address of the mgr, cast to (pool_pt)
    if(!pool_store){
        mem_init();
    }

    if(pool_store_size + 1 >= pool_store_capacity){
        _mem_resize_pool_store();
    }

    if(pool_store){
        _mem_resize_pool_store();
        pool_mgr_pt   pmt;
        pmt = (pool_mgr_pt) pool_store;

        if(pmt) {
            pool_t pl;
            pl.alloc_size = 0;
            pl.policy = policy;
            pl.num_gaps = 0;
            pl.total_size = size;
            pl.num_allocs = 0;

            pmt->pool = pl;

            if(!&pmt->pool){
                free(pmt);
                return NULL;
            }

            pmt->node_heap = (node_pt) malloc(sizeof(node_t)*MEM_NODE_HEAP_INIT_CAPACITY);
            for(size_t i=0;i  < MEM_NODE_HEAP_INIT_CAPACITY;i++){
                node_t nd;
                pmt->node_heap[i] = nd;
                pmt->node_heap[i].allocated = 0;
                pmt->node_heap[i].used = 0;
                pmt->node_heap[i].prev = NULL;
                pmt->node_heap[i].next = NULL;
                alloc_t all;
                all.size = 0;
                all.mem = NULL;
                pmt->node_heap[i].alloc_record = all;
            }

            if(!pmt->node_heap){
                free(&pmt->pool);
                free(pmt);
                return NULL;
            }
            pmt->total_nodes = MEM_NODE_HEAP_INIT_CAPACITY;

            pmt->gap_ix = (gap_pt) malloc(sizeof(gap_t) * MEM_GAP_IX_INIT_CAPACITY);
            for(size_t i =0; i < MEM_GAP_IX_INIT_CAPACITY; i++){
                gap_t ng;
                pmt->gap_ix[i] = ng;
                pmt->gap_ix[i].node = NULL;
                pmt->gap_ix[i].size = 0;
            }

            if(_mem_add_to_gap_ix(pmt,pmt->pool.total_size,&pmt->node_heap[0])==ALLOC_FAIL){
                printf("FAILED TO ALLOCATE");
                return NULL;
            }else{
                pmt->node_heap[0].used = 1;
                pmt->used_nodes += 1;
            }

            if(!pmt->gap_ix) {
                if(pmt->node_heap != NULL){
                    for(size_t t = pmt->total_nodes; t>0 ; t--){
                        if(&pmt->node_heap[t] != NULL){
                            free(&pmt->node_heap[t]);
                        }
                    }
                }
                mem_pool_close(&pmt->pool);
                free(&pmt->node_heap);
                free(pmt);
            }

            (*pool_store) = pmt;
        }
    }
    return (pool_pt) pool_store;
}

alloc_status mem_pool_close(pool_pt pool) {
    // get mgr from pool by casting the pointer to (pool_mgr_pt)
    // check if this pool is allocated
    // check if pool has only one gap
    // check if it has zero allocations
    // free memory pool
    // free node heap
    // free gap index
    // find mgr in pool store and set to null
    // note: don't decrement pool_store_size, because it only grows
    // free mgr

    pool_mgr_pt pmt = *pool_store;

    if(pool != NULL){

        if (pool->num_gaps > 0 && pool->alloc_size == 0) {
            free(pmt->gap_ix);


            if (pool->mem) {
                free(pool->mem);
            }

            for (size_t t = pmt->total_nodes - 1; t > -1; t--) {
                if (pmt->node_heap[t].alloc_record.mem != NULL) {
                    free(pmt->node_heap[t].alloc_record.mem);
                }
                free(&pmt->node_heap[t]);
            }

            pmt->node_heap = NULL;

            gap_pt gpt = pmt->gap_ix;

            for (size_t t = pmt->gap_ix_capacity - 1; t > -1; t--) {
                free(&pmt->gap_ix[t]);
            }
            pmt->gap_ix = NULL;

            pool_pt ppt = NULL;
            int poolIdx = 0;
            int found = 0;
            for (size_t t = 0; t < pool_store_size; t++) {
                if (pool_store == &pmt) {
                    free(&pool_store[t]->pool);
                    pool_store[t] = NULL;
                    break;
                }
            }
            return ALLOC_OK;
        }
    }

    return ALLOC_NOT_FREED;
}

alloc_pt mem_new_alloc(pool_pt pool, size_t size) {
    // get mgr from pool by casting the pointer to (pool_mgr_pt)
    // check if any gaps, return null if none
    // expand heap node, if necessary, quit on error
    // check used nodes fewer than total nodes, quit on error
    // get a node for allocation:
    // if FIRST_FIT, then find the first sufficient node in the node heap
    // if BEST_FIT, then find the first sufficient node in the gap index
    // check if node found
    // update metadata (num_allocs, alloc_size)
    // calculate the size of the remaining gap, if any11111
    // remove node from gap index
    // convert gap_node to an allocation node of given size
    // adjust node heap:
    //   if remaining gap, need a new node
    //   find an unused one in the node heap
    //   make sure one was found
    //   initialize it to a gap node
    //   update metadata (used_nodes)
    //   update linked list (new node right after the node for allocation)
    //   add to gap index
    //   check if successful
    // return allocation record by casting the node to (alloc_pt)
    pool_mgr_pt pmt = (pool_mgr_pt) pool;


    if(pool->num_gaps > 0){

        //expand node heap if necessary
        _mem_resize_node_heap(pmt);

        //check if used is less than total nodes
        if(pmt->used_nodes >= pmt->total_nodes){
            return NULL;
        }

        //iterate through gaps and look for a sufficient gap to allocate to
        gap_pt bestGap = NULL;
        size_t bestSize = pool->total_size;
        int idx = 0;


        for(size_t i = 0; i < pool->num_gaps;i++){
            if(pmt->gap_ix[i].size >= size){
                if(pool->policy == FIRST_FIT){
                    //this is the first fit break loop for assignment
                        idx = i;
                        bestSize = pmt->gap_ix[i].size;
                        break;
                }else {
                    //check for best fit and set variables
                    if (pmt->gap_ix[i].size < bestSize){
                        idx = 0;
                        bestSize = pmt->gap_ix[i].size;
                    }
                }
            }
        }

        if(idx < pool->num_gaps) {
            bestGap = &pmt->gap_ix[idx];
            node_pt bestNode = bestGap->node;

            if (bestNode->alloc_record.mem != NULL) {
                free(bestNode->alloc_record.mem);
            }

            //if gap found, allocate new memory
            if (bestNode->alloc_record.size >= size) {
                int remainingSize = bestSize - size;

                _mem_remove_from_gap_ix(pmt, bestGap->size, bestNode);

                //convert to allocation node
                bestNode->allocated = 1;
                bestNode->alloc_record.size = size;
                bestNode->used = 1;

                pool->num_allocs += 1;


                //if remaining size, add gap node record
                if (remainingSize > 0) {
                    //create a new gap
                    node_pt newGapNode = NULL;
                    int idx = 0;
                    while (idx < pmt->total_nodes) {
                        if (pmt->node_heap[idx].allocated == 0) {
                            break;
                        }
                        idx += 1;
                    }

                    newGapNode = &pmt->node_heap[idx];

                    if (newGapNode) {
                        if (_mem_add_to_gap_ix(pmt, remainingSize, newGapNode) == ALLOC_OK) {
                            newGapNode->used = 1;
                            pmt->used_nodes += 1;

                            //pmt->pool.num_allocs -= 1;
                            newGapNode->prev = bestNode;
                            newGapNode->next = bestNode->next;
                            bestNode->next = newGapNode;


                            if (newGapNode->next != NULL) {
                                newGapNode->next->prev = newGapNode;
                            }
                        } else {
                            return NULL;
                        }
                    }
                }


                bestNode->alloc_record.mem = (char *) malloc(sizeof(char) * size);
                pool->alloc_size += size;
                return (alloc_pt) bestNode;
            }
        }
    }

    return NULL;
}

alloc_status mem_del_alloc(pool_pt pool, alloc_pt alloc) {
    // get mgr from pool by casting the pointer to (pool_mgr_pt)
    // get node from alloc by casting the pointer to (node_pt)
    // find the node in the node heap
    // this is node-to-delete
    // make sure it's found
    // convert to gap node
    // update metadata (num_allocs, alloc_size)
    // if the next node in the list is also a gap, merge into node-to-delete
    //   remove the next node from gap index
    //   check success
    //   add the size to the node-to-delete
    //   update node as unused
    //   update metadata (used nodes)
    //   update linked list:
    /*
                    if (next->next) {
                        next->next->prev = node_to_del;
                        node_to_del->next = next->next;
                    } else {
                        node_to_del->next = NULL;
                    }
                    next->next = NULL;
                    next->prev = NULL;
     */

    // this merged node-to-delete might need to be added to the gap index
    // but one more thing to check...
    // if the previous node in the list is also a gap, merge into previous!
    //   remove the previous node from gap index
    //   check success
    //   add the size of node-to-delete to the previous
    //   update node-to-delete as unused
    //   update metadata (used_nodes)
    //   update linked list
    /*
                    if (node_to_del->next) {
                        prev->next = node_to_del->next;
                        node_to_del->next->prev = prev;
                    } else {
                        prev->next = NULL;
                    }
                    node_to_del->next = NULL;
                    node_to_del->prev = NULL;
     */
    //   change the node to add to the previous node!
    // add the resulting node to the gap index
    // check success

    pool_mgr_pt pmpt = (pool_mgr_pt) pool;


    node_pt node = NULL;
    int idx = 0;

    while( idx < pmpt->total_nodes){
        if(&pmpt->node_heap[idx].alloc_record == alloc){
            break;
        }
        idx += 1;
    }


    if(idx < pmpt->total_nodes){
        node = &pmpt->node_heap[idx];
    }else{
        printf("Failed to Find Allocation\n");
        return ALLOC_FAIL;
    }

    if(node != NULL && node->allocated == 1){
        //convert
        gap_t gap;

        _mem_add_to_gap_ix(pmpt,node->alloc_record.size,node);
        pool->num_allocs -= 1;
        pool->alloc_size -= alloc->size;
        node->allocated = 0;
        gap.size = node->alloc_record.size;

        //free allocated memory
        if(alloc->mem != NULL){
            free(alloc->mem);
            alloc->mem = NULL;
        }

        //check if next node is not allocated and then remove if not
        while(node->next != NULL && node->next->allocated == 0 && node->next->used == 1){
            //find in gap_idx
            int gapIdx =0;
            for(size_t i = 0 ; i < pool->num_gaps;i++){
                if(pmpt->gap_ix[i].node == node->next){
                    break;
                }
                gapIdx++;
            }

            int currGapIdx = 0;
            for(size_t i = 0 ; i < pool->num_gaps;i++){
                if(pmpt->gap_ix[i].node == node){
                    break;
                }
                currGapIdx++;
            }

            if(node->next != NULL){
                //resize
                node->alloc_record.size += node->next->alloc_record.size;
                pmpt->gap_ix[currGapIdx].size += node->next->alloc_record.size;


                //remove from gap index
                node_pt gpn = pmpt->gap_ix[gapIdx].node;
                _mem_remove_from_gap_ix(pmpt,pool->num_gaps,pmpt->gap_ix[gapIdx].node);

                //reset meta data
                gpn->used = 0;
                gpn->allocated = 0;
                gpn->alloc_record.size = 0;
                pmpt->used_nodes -= 1;


                //reset pointers
                node->next = gpn->next;
                if(gpn->next){
                    gpn->next->prev = node;
                }

                gpn->next = NULL;
                gpn->prev = NULL;

                if(&gpn->alloc_record != NULL){
                    if(gpn->alloc_record.mem != NULL){
                        free(gpn->alloc_record.mem);
                        gpn->alloc_record.mem = NULL;
                    }

                }
            }
        }



        //check if previous node is a gap
        while(node->prev != NULL && &node->prev->alloc_record != NULL) {
            if(node->prev->allocated != 0){
                break;
            }
            node_pt gpn = node->prev;

            gpn->alloc_record.size += node->alloc_record.size;
            gpn->next = node->next;
            if(gpn->next != NULL){
                gpn->next->prev  = gpn;
            }

            _mem_remove_from_gap_ix(pmpt,node->alloc_record.size,node);
            node->next = NULL;
            node->prev = NULL;
            node->alloc_record.size = 0;
            node->used = 0;
            pmpt->used_nodes -= 1;

            node = gpn;
        }
        return ALLOC_OK;
    }
    return ALLOC_FAIL;
}

void mem_inspect_pool(pool_pt pool,pool_segment_pt *segments, unsigned *num_segments) {
    // get the mgr from the pool
    // allocate the segments array with size == used_nodes
    // check successful
    // loop through the node heap and the segments array
    //    for each node, write the size and allocated in the segment
    // "return" the values:
    /*
                    *segments = segs;
                    *num_segments = pool_mgr->used_nodes;
     */

    pool_mgr_pt pmpt = (pool_mgr_pt) pool;

    pool_segment_pt segs;

    segs = (pool_segment_pt) malloc(sizeof(pool_segment_t) * pmpt->used_nodes);
    int segIdx = 0;

    for(int i = 0;i< pmpt->total_nodes;i++){
        if(&pmpt->node_heap[i] != NULL && pmpt->node_heap[i].used){
            pool_segment_t seg;
            segs[segIdx] = seg;
            segs[segIdx].size = pmpt->node_heap[i].alloc_record.size;
            segs[segIdx].allocated = pmpt->node_heap[i].allocated;
            segIdx += 1;
        }

    }

    *segments = segs;
    *num_segments = pmpt->used_nodes;
}



/***********************************/
/*                                 */
/* Definitions of static functions */
/*                                 */
/***********************************/
static alloc_status _mem_resize_pool_store() {
    // check if necessary
    /*
                if (((float) pool_store_size / pool_store_capacity)
                    > MEM_POOL_STORE_FILL_FACTOR) {...}
     */

    if(pool_store_size == pool_store_capacity){
        return ALLOC_FAIL;
    }

    // don't forget to update capacity variables
    if(((float)pool_store_size / pool_store_capacity) > MEM_POOL_STORE_FILL_FACTOR){
        pool_store = realloc(pool_store,sizeof(pool_store) *MEM_POOL_STORE_EXPAND_FACTOR);
        for(size_t t = pool_store_size; t < pool_store_size * MEM_POOL_STORE_EXPAND_FACTOR;t++){
                pool_store[t] = NULL;
        }
    }


    pool_store_size = pool_store_size * MEM_POOL_STORE_EXPAND_FACTOR;
    if(pool_store){
        return ALLOC_OK;
    }
    return ALLOC_FAIL;
}

static alloc_status _mem_resize_node_heap(pool_mgr_pt pool_mgr) {
    // see above
    if(pool_mgr->used_nodes / pool_mgr->total_nodes > MEM_NODE_HEAP_FILL_FACTOR){
        pool_mgr->node_heap = realloc(pool_mgr->node_heap, sizeof(node_t)* pool_mgr->total_nodes * MEM_NODE_HEAP_EXPAND_FACTOR);

        pool_mgr->total_nodes *= MEM_NODE_HEAP_EXPAND_FACTOR;

        //reserve initial space for our nodes
        for(int i = pool_mgr->used_nodes; i < pool_mgr->total_nodes;i++){
            if(&pool_mgr->node_heap[i]) {
                node_t nd;
                pool_mgr->node_heap[i] = nd;
                pool_mgr->node_heap[i].allocated = 0;
                pool_mgr->node_heap[i].used = 0;
                alloc_t all;
                all.size = 0;
                pool_mgr->node_heap[i].alloc_record = all;
            }
        }

    }
    return ALLOC_OK;
}

static alloc_status _mem_resize_gap_ix(pool_mgr_pt pool_mgr) {
    // see above
    if((*pool_mgr).pool.num_gaps >= MEM_GAP_IX_FILL_FACTOR) {
        (*pool_store)->gap_ix = realloc((*pool_store)->gap_ix,
                                        sizeof((*(*pool_store)->gap_ix)) * MEM_GAP_IX_EXPAND_FACTOR);
        pool_mgr->gap_ix_capacity = MEM_GAP_IX_EXPAND_FACTOR * pool_mgr->pool.num_gaps;
        for (size_t i = 0; i < pool_mgr->pool.num_gaps; i++) {
            gap_t ng;
            pool_mgr->gap_ix[i] = ng;
            pool_mgr->gap_ix[i].node = NULL;
            pool_mgr->gap_ix[i].size = 0;
        }
    }

    return ALLOC_FAIL;
}

static alloc_status _mem_add_to_gap_ix(pool_mgr_pt pool_mgr,
                                       size_t size,
                                       node_pt node) {

    // expand the gap index, if necessary (call the function)
    // add the entry at the end
    // update metadata (num_gaps)
    // sort the gap index (call the function)
    // check success
    if(pool_mgr->pool.num_gaps + 1 == pool_mgr->gap_ix_capacity) {
        _mem_resize_gap_ix(pool_mgr);
    }


    if(pool_mgr == NULL || pool_mgr->gap_ix == NULL){
        return ALLOC_FAIL;
    }

    pool_mgr->pool.num_gaps += 1;
    node->alloc_record.size = size;
    gap_t  gp;
    gp.node = node;
    gp.size = node->alloc_record.size;
    int gapIdx = pool_mgr->pool.num_gaps-1;
    pool_mgr->gap_ix[gapIdx]= gp;

    if(node->alloc_record.mem != NULL){
        free(node->alloc_record.mem);
        node->alloc_record.mem = NULL;
    }

    node->allocated = 0;
    node->used = 1;

    if(pool_mgr->gap_ix[pool_mgr->pool.num_gaps-1].node != node){
        return ALLOC_FAIL;
    }

    _mem_sort_gap_ix(pool_mgr);
    return ALLOC_OK;
}

static alloc_status _mem_remove_from_gap_ix(pool_mgr_pt pool_mgr,
                                            size_t size,
                                            node_pt node) {
    // find the position of the node in the gap index
    // loop from there to the end of the array:
    //    pull the entries (i.e. copy over) one position up
    //    this effectively deletes the chosen node
    // update metadata (num_gaps)
    // zero out the element at position num_gaps!
    gap_pt gpt = NULL;

    int idx = 0;
    while(idx < pool_mgr->pool.num_gaps){
        if(pool_mgr->gap_ix[idx].node ==  node){
            break;
        }

        idx +=1;
    }

    if(idx < pool_mgr->pool.num_gaps) {
        gpt = &pool_mgr->gap_ix[idx];


        for (size_t i = idx; i < pool_mgr->pool.num_gaps; i++) {
            pool_mgr->gap_ix[i] = pool_mgr->gap_ix[i+1];
        }

        //free(gpt);
    }else{
        printf("REMOVING FROM GAP IX FAILED\n");
        return ALLOC_FAIL;
    }


    if(pool_mgr->gap_ix == NULL){
        return ALLOC_FAIL;
    }
    pool_mgr->pool.num_gaps -= 1;
    return ALLOC_OK;
}

// note: only called by _mem_add_to_gap_ix, which appends a single entry
static alloc_status _mem_sort_gap_ix(pool_mgr_pt pool_mgr) {
    // the new entry is at the end, so "bubble it up"
    // loop from num_gaps - 1 until but not including 0:
    //    if the size of the current entry is less than the previous (u - 1)
    //    or if the sizes are the same but the current entry points to a
    //    node with a lower address of pool allocation address (mem)
    //       swap them (by copying) (remember to use a temporary variable)


    gap_pt gidx = pool_mgr->gap_ix;

    if (gidx != NULL) {


       for(int i = pool_mgr->pool.num_gaps-1 ; i > 0 ; i--){
           if(gidx[i].size < gidx[i-1].size){
               gap_t temp = gidx[i];
               gidx[i] = gidx[i-1];
               gidx[i-1] = temp;
           }
        }

        return ALLOC_OK;
    }

    return ALLOC_FAIL;
}


