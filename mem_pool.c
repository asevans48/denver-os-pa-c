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

    if(!pool_store){
        pool_store_capacity = MEM_POOL_STORE_INIT_CAPACITY;
        pool_store = (pool_mgr_pt *) malloc(sizeof(pool_pt) * MEM_POOL_STORE_INIT_CAPACITY);
        for(size_t t =0; t< MEM_POOL_STORE_INIT_CAPACITY;t++){
            pool_store[t] = NULL;
        }
        return ALLOC_OK;
    }

    return ALLOC_CALLED_AGAIN;

}

alloc_status mem_free() {
    // ensure that it's called only once for each mem_init
    // make sure all pool managers have been deallocated
    // can free the pool store array
    // update static variables
    if(pool_store != NULL){
        for(size_t i = 0; i< pool_store_size;i++){
            if(pool_store[i] != NULL) {

                for(size_t i = pool_store[i]->pool.num_gaps ; i >= 0; i--){
                    free(&pool_store[i]->gap_ix[i]);
                }

                mem_pool_close(&pool_store[i]->pool);
                node_pt pt = pool_store[i]->node_heap;
                for(size_t i = pool_store[i]->total_nodes; i >= 0; i--){
                    free(&pool_store[i]->node_heap[i]);
                }

                free(&pool_store[i]);
                pool_store[i] = NULL;
            }
        }
        free(pool_store);
        pool_store = NULL;
        pool_store_size = 0;
        pool_store_capacity = MEM_POOL_STORE_INIT_CAPACITY;
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
    int idx = 0;
    if(!pool_store){
        if(mem_init() == ALLOC_CALLED_AGAIN){
            return NULL;
        }
    }else{
        _mem_resize_pool_store();
        while(pool_store[idx] != NULL){
            idx +=1;
        }

        if(idx >= pool_store_capacity){
            printf("NEW POOL NOT FOUND\n");
            return NULL;
        }

        pool_t newPool;
        newPool.alloc_size = 0;
        newPool.mem = (char *) malloc(sizeof(char *) *size);
        newPool.num_allocs = 0;
        newPool.num_gaps = 1;
        newPool.policy = policy;
        newPool.total_size = size;
        pool_store[idx] = (pool_mgr_pt) malloc(sizeof(pool_mgr_t));

        pool_store[idx]->pool = newPool;

        pool_store[idx]->gap_ix = (gap_pt) malloc(sizeof(gap_t)*MEM_GAP_IX_INIT_CAPACITY);
        pool_store[idx]->node_heap = (node_pt) malloc(sizeof(node_t) * MEM_NODE_HEAP_INIT_CAPACITY);
        pool_store[idx]->total_nodes = MEM_NODE_HEAP_INIT_CAPACITY;
        pool_store[idx]->gap_ix_capacity = MEM_GAP_IX_INIT_CAPACITY;
        pool_store[idx]->used_nodes = 1;


        if(&pool_store[idx]->pool == NULL){
            printf("FAILED TO ALLOCATE POOL\n");
            free(&pool_store[idx]->pool);
            return NULL;
        }

        for(size_t i =0 ; i< MEM_NODE_HEAP_INIT_CAPACITY; i++){
            node_t newNode;
            newNode.allocated = 0;
            alloc_t alloc;
            alloc.mem = NULL;
            alloc.size = 0;
            newNode.alloc_record = alloc;
            if(i == 0){
                newNode.used = 1;
            }else{
                newNode.used = 0;
            }
            newNode.next = NULL;
            newNode.prev = NULL;
            pool_store[idx]->node_heap[i] = newNode;
        }

        if(pool_store[idx]->node_heap == NULL){
            printf("NODE HEAP NOT SET\n");
            free(&pool_store[idx]->pool);
            return NULL;
        }


        for(size_t i = 0; i < MEM_GAP_IX_INIT_CAPACITY; i++){
            if(i ==0){
                gap_t gp;
                gp.node = &(pool_store[idx]->node_heap[0]);
                gp.size = pool_store[idx]->pool.total_size;
                pool_store[idx]->gap_ix[0] = gp;
            }else {
                gap_t gp;
                gp.size = 0;
                gp.node = NULL;
                pool_store[idx]->gap_ix[i] = gp;
            }
        }


        if(pool_store[idx]->gap_ix == NULL){
            node_pt  nd = pool_store[idx]->node_heap;

            for(size_t i = 0; i < MEM_NODE_HEAP_INIT_CAPACITY; i++){
                free(&nd[i]);
            }

            pool_store[idx]->node_heap = NULL;
            pool_t pl = pool_store[idx]->pool;
            free(pl.mem);
            return NULL;
        }
    }
    pool_pt pt =&(pool_store[idx]->pool);
    return pt;
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

    pool_mgr_pt pmt = (pool_mgr_pt) pool_store;

    if(pmt != NULL){
        unsigned int allocs = 0;
        unsigned int gaps = 1;
        if(pool->num_allocs != allocs && pool->num_gaps == gaps){
            return ALLOC_NOT_FREED;
        }

        free(pool->mem);

        for(size_t i = 0; i < pmt->total_nodes && &pmt->node_heap[i] != NULL; i++){
            if(&pmt->node_heap[i].alloc_record != NULL){
                if(pmt->node_heap[i].alloc_record.mem != NULL){
                    free(pmt->node_heap[i].alloc_record.mem);
                }
                free(&pmt->node_heap[i].alloc_record);
                alloc_pt alloc = NULL;
                pmt->node_heap[i].alloc_record = *alloc;
            }
            pmt->node_heap[i].next = NULL;
            pmt->node_heap[i].prev = NULL;
            free(&pmt->node_heap[i]);
        }
        free(pmt->node_heap);

        //free gap IDX
        for(size_t i = 0; i < pmt->pool.num_gaps && &pmt[i] != NULL;i++){
            gap_pt gpt = &pmt->gap_ix[i];
            gpt->node = NULL;
            free(gpt);
        }


        free(pmt->gap_ix);

        if(&pool != NULL) {
            free(pool);
            pool = NULL;
        }

        pmt= NULL;
        free(pmt);


        return ALLOC_OK;
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
        _mem_resize_node_heap(pmt);

        if(pmt->used_nodes < pmt->total_nodes){
            size_t idx =0;
            size_t bestIdx = 0;
            unsigned int bestSize = pmt->pool.total_size;
            while(idx < pmt->pool.num_gaps){
                if(pmt->gap_ix[idx].size >= size){
                    if(pool->policy != BEST_FIT){
                        bestIdx = idx;
                         break;
                    }else{
                        if(pmt->gap_ix[idx].size < bestSize){
                            bestIdx = idx;
                            bestSize = pmt->gap_ix[idx].size;
                        }
                    }
                }
                idx += 1;
            }

            idx = bestIdx;
            if(idx < pmt->pool.num_gaps){
                //assign new allocation
                node_pt bestNode = pmt->gap_ix[idx].node;
                size_t remaining = pmt->gap_ix[idx].size - size;

                if(remaining > pool->total_size - pool->alloc_size){
                    printf("Remaining Allocation Size Too Large\n");
                    return NULL;
                }

                _mem_remove_from_gap_ix(pmt,size,bestNode);

                bestNode->allocated = 1;
                bestNode->alloc_record.size = size;
                bestNode->alloc_record.mem = pmt->pool.mem;
                pmt->pool.alloc_size += size;



                //take remaining and assign new gap if so
                if(remaining > 0){

                    //get new node
                    size_t ridx = 0;
                    while(idx < pmt->total_nodes && pmt->node_heap[ridx].used == 1){
                        ridx += 1;
                    }

                    if(ridx < pmt->total_nodes && pmt->node_heap[ridx].used == 0) {
                        if(&pmt->node_heap[ridx] == NULL){
                           return NULL;
                        }else{
                            pmt->node_heap[ridx].alloc_record.size = 0;
                            pmt->node_heap[ridx].alloc_record.mem = pool->mem;
                            pmt->node_heap[ridx].used = 1;
                        }
                        pmt->used_nodes += 1;

                        //add to gap_ix
                        node_pt ng = &pmt->node_heap[ridx];
                        ng->used = 1;
                        ng->allocated = 0;
                        ng->next = NULL;
                        ng->prev = NULL;
                        alloc_t nalloc;
                        nalloc.size =0;
                        nalloc.mem = pmt->pool.mem;
                        ng->alloc_record = nalloc;

                        if(_mem_add_to_gap_ix(pmt, remaining,ng) == ALLOC_FAIL){
                            printf("FAILED TO ADD TO GAPS\n");
                            return NULL;
                        }

                        //reset pointers
                        //set next of gap to actual next
                        ng->next = bestNode->next;
                        if(bestNode->next != NULL){
                            bestNode->next->prev = ng;
                        }

                        ng->prev = bestNode;
                        bestNode->next = ng;

                        //check next to see if gap and fold
                        if(ng->next != NULL && ng->next->allocated == 0 && ng->next->used ==1){
                            size_t gapIdx = 0;
                            while(gapIdx < pool->num_gaps && &pmt->gap_ix[gapIdx] != NULL && pmt->gap_ix[gapIdx].node != pmt->node_heap[idx].next){
                                gapIdx += 1;
                            }

                            //roll gap
                            if(gapIdx < pool->num_gaps) {
                                pmt->gap_ix[pmt->pool.num_gaps - 1].size += pmt->gap_ix[gapIdx].size;
                                pmt->gap_ix[gapIdx].node->used = 0;
                                _mem_remove_from_gap_ix(pmt,pmt->gap_ix[gapIdx].size,pmt->gap_ix[gapIdx].node);
                                pmt->used_nodes -= 1;
                            }else{
                                printf("Failed to Find Next Gap\n");
                                return NULL;
                            }
                        }

                    }else{
                        printf("IDX > Available Nodes\n");
                        return NULL;
                    }
                }
                bestNode->alloc_record.mem = (char *) malloc(sizeof(char)*size);
                pmt->pool.num_allocs += 1;
                return (alloc_pt) bestNode;
            }
        }
    }
    printf("WARNING: No Allocation\n");
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

    pool_mgr_pt pmt = (pool_mgr_pt) pool;

    //get allocation node and heap index
    node_pt  node = (node_pt) alloc;

    if(node != NULL){
        //setup gap vars
        size_t gap_size = node->alloc_record.size;
        _mem_add_to_gap_ix(pmt,gap_size,node);
        pmt->pool.num_allocs -= 1;
        pmt->pool.alloc_size -= node->alloc_record.size;


        //zero out node vars
        node->allocated = 0;
        //look for next gaps and roll
        while(node->next != NULL && node->next->allocated == 0 && node->next->used == 1) {
            //get gap from gap index
            size_t gapIdx = pool->num_gaps;
            size_t nextGapIdx = pool->num_gaps;
            for (size_t i = 0; i < pmt->gap_ix_capacity; i++) {
                if (pmt->gap_ix[i].node == node) {
                    gapIdx = i;
                }else if(pmt->gap_ix[i].node == node->next){
                   nextGapIdx = i;
                }
            }

            if (gapIdx < pmt->pool.num_gaps && gapIdx < pool->num_gaps) {
                //remove gap and reset size
                pmt->gap_ix[gapIdx].size += pmt->gap_ix[nextGapIdx].size;


                node_pt next = node->next;
                if(next->next != NULL){
                   next->next->prev = node;
                }
                node->next = next->next;
                next->next= NULL;
                next->prev = NULL;
                next->used = 0;
                _mem_remove_from_gap_ix(pmt,pmt->gap_ix[nextGapIdx].size,next);
                pmt->used_nodes -= 1;
            } else {
                printf("Failed to Roll Into Next Node\n");
                return ALLOC_FAIL;
            }

        }

        //look for previous gaps and roll or add to gap index
        while(node->prev != NULL && node->prev->used == 1 && node->prev->allocated ==0){
            // /get gap from gap index
            size_t gapIdx = pool->num_gaps;
            size_t prevGapIdx = pool->num_gaps;
            for (size_t i = 0; i < pmt->gap_ix_capacity; i++) {
                if (pmt->gap_ix[i].node == node) {
                    gapIdx = i;
                }else if(pmt->gap_ix[i].node == node->prev){
                    prevGapIdx = i;
                }
            }

            if(gapIdx < pmt->pool.num_gaps && prevGapIdx < pmt->pool.num_gaps){
                pmt->gap_ix[prevGapIdx].size += pmt->gap_ix[gapIdx].size;

                node_pt oldNode = node;
                node = oldNode->prev;

                //reset pointers
                node->next = oldNode->next;

                if(oldNode->next != NULL){
                    oldNode->next->prev = node;
                }
                oldNode->next = NULL;
                oldNode->prev = NULL;

                //reset values
                oldNode->used = 0;
                _mem_remove_from_gap_ix(pmt,pmt->gap_ix[gapIdx].size,oldNode);
                pmt->used_nodes -= 1;

            }else{
                return ALLOC_FAIL;
            }

        }


        return ALLOC_OK;

    }


    printf("Allocation FAILED\n");
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

    pool_mgr_pt  pmt = (pool_mgr_pt) pool;
    pool_segment_pt  segPtrs =  calloc(pmt->used_nodes,sizeof(pool_segment_t));
    *num_segments = pmt->used_nodes;

    if(pmt->used_nodes > 0 && segPtrs != NULL){

        if(segPtrs) {
            size_t idx = 0;
            node_pt node = &pmt->node_heap[0];//head node
            while (node->used == 0) {
                node = node->next;
            }

            pool_segment_pt  allocSeg = segPtrs;

            for (size_t i = 0; i < pmt->used_nodes; i++) {
                if (node->allocated == 1) {
                    allocSeg->size = node->alloc_record.size;
                } else {
                    gap_pt gp = &pmt->gap_ix[0];
                    while (gp->node != node) {
                        gp++;
                    }
                    allocSeg->size = gp->size;
                }
                allocSeg->allocated = node->allocated;
                allocSeg++;
                node = node->next;

            }

            *segments = segPtrs;
        }
    }

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

    if(((float)pool_store_size / pool_store_capacity) > MEM_POOL_STORE_FILL_FACTOR){
        unsigned int oldSize = pool_store_size;
        pool_store_size = pool_store_size * MEM_POOL_STORE_EXPAND_FACTOR;
        pool_store = (pool_mgr_pt  *)realloc(pool_store,sizeof(pool_mgr_pt) * pool_store_size);
        for(size_t i = oldSize; i< pool_store_size; i++){
            pool_store[i] = NULL;
        }

        if(pool_store != NULL){
            return ALLOC_OK;
        }
    }

    return ALLOC_FAIL;
}

static alloc_status _mem_resize_node_heap(pool_mgr_pt pool_mgr) {
    // see above
    if(((float)pool_mgr->used_nodes / pool_mgr->total_nodes) > MEM_NODE_HEAP_FILL_FACTOR){
        pool_mgr->total_nodes *= MEM_NODE_HEAP_EXPAND_FACTOR;
        pool_mgr->node_heap = (node_pt) realloc(pool_mgr->node_heap, sizeof(node_pt) * pool_mgr->total_nodes);

        for(size_t i = pool_mgr->used_nodes ; i < pool_mgr->total_nodes; i++){
            node_t nd;
            nd.used = 0;
            nd.allocated = 0;
            pool_mgr->node_heap[i] = nd;
        }

        if(pool_mgr == NULL){
            return ALLOC_FAIL;
        }
    }
    return ALLOC_OK;
}

static alloc_status _mem_resize_gap_ix(pool_mgr_pt pool_mgr) {
    // see above
    if(((float)pool_mgr->pool.num_gaps /pool_mgr->gap_ix_capacity) > MEM_GAP_IX_FILL_FACTOR){
        unsigned int oldSize = pool_mgr->gap_ix_capacity;
        pool_mgr->gap_ix_capacity *= MEM_GAP_IX_EXPAND_FACTOR;
        pool_mgr->gap_ix = (gap_pt) realloc(pool_mgr->gap_ix,sizeof(gap_pt)*pool_mgr->gap_ix_capacity);

        for(size_t i = oldSize; i < pool_mgr->gap_ix_capacity;i++){
            gap_t gp;
            gp.node = NULL;
            gp.size = 0;
            pool_mgr->gap_ix[i] = gp;
        }
        if(pool_mgr->gap_ix != NULL) {
            return ALLOC_OK;
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
    _mem_resize_gap_ix(pool_mgr);

    gap_t gap;
    gap.node = node;
    gap.size = size;


    pool_mgr->gap_ix[pool_mgr->pool.num_gaps] = gap;

    gap_pt cg = &pool_mgr->gap_ix[pool_mgr->pool.num_gaps];

    if(cg->node != node){
        printf("Failed to Allocate to Gaps\n");
        return ALLOC_FAIL;
    }

    pool_mgr->pool.num_gaps += 1;
    return _mem_sort_gap_ix(pool_mgr);
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
    size_t idx =0;
    while(idx < pool_mgr->pool.num_gaps && pool_mgr->gap_ix[idx].node != node){
        idx += 1;
    }

    if(idx == pool_mgr->pool.num_gaps){
        printf("Failed to Remove Gap\n");
        return ALLOC_FAIL;
    }

    for(size_t i = idx;i < pool_mgr->pool.num_gaps;i++){
        pool_mgr->gap_ix[i] = pool_mgr->gap_ix[i+1];
    }

    gap_t gap;
    gap.size =0;
    gap.node = NULL;
    pool_mgr->gap_ix[pool_mgr->pool.num_gaps] = gap;


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

    size_t start = pool_mgr->pool.num_gaps - 1;
    while(start > 0 && pool_mgr->gap_ix[start].size < pool_mgr->gap_ix[start -1].size){
        gap_t temp = pool_mgr->gap_ix[start];
        pool_mgr->gap_ix[start] = pool_mgr->gap_ix[start -1];
        start -=1;
        pool_mgr->gap_ix[start] = temp;
    }


    //using start loop until lowest address
    while(start > 0 && pool_mgr->gap_ix[start].size == pool_mgr->gap_ix[start - 1].size && pool_mgr->gap_ix[start].node->alloc_record.mem < pool_mgr->gap_ix[start -1].node->alloc_record.mem){
        gap_t temp = pool_mgr->gap_ix[start];
        pool_mgr->gap_ix[start] = pool_mgr->gap_ix[start -1];
        start -=1;
        pool_mgr->gap_ix[start] = temp;
    }

    return ALLOC_OK;
}


