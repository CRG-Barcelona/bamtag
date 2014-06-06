#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* The basic bam reading procedure. */
#include <htslib/sam.h>

int main(int argc, char *argv[])
/* Process command line. */
{
    samFile *s = sam_open(argv[1], "r");
    samFile *out = sam_open(argv[2], "wb");
    bam_hdr_t *h = sam_hdr_read(s);
    bam1_t *b = bam_init1();
    sam_hdr_write(out, h);
    while (sam_read1(s, h, b) >= 0)
    {
	bam1_core_t *core = &b->core;
	if (core->flag & BAM_FMUNMAP)
	    sam_write1(out, h, b);
    }
    bam_hdr_destroy(h);
    sam_close(s);
    sam_close(out);
}
