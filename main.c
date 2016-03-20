#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "test_suite.h"
#include "mem_pool.h"

/* forward declarations */
void print_pool(pool_pt pool);

/* main */
int main(int argc, char *argv[]) {
    run_test_suite();
    return 0;
}

/* function definitions */
void print_pool(pool_pt pool) {
    pool_segment_pt segs = NULL;
    unsigned size = 0;

    assert(pool);

    mem_inspect_pool(pool, &segs, &size);
    printf("Number of Nodes %d\n",size);

    printf("Asserting Segs\n\n");
    assert(segs);

    printf("Assert Size\n");
    assert(&size);


    for (unsigned u = 0; u < size; u ++) {
        printf("Used Status %lu\n",segs[u].allocated);
        printf("%10lu - %s\n", (unsigned long) segs[u].size, (segs[u].allocated) ? "alloc" : "gap");
    }
    free(segs);

    printf("\n");

}
