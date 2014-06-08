#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* bamtag - tag BAMS our typical way */

#include <jkweb/common.h>
#include <jkweb/dystring.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/rangeTree.h>
#include <jkweb/genomeRangeTree.h>
#include <jkweb/bigWig.h>

#include <htslib/sam.h>

#include <limits.h>

#define ZERO_QUALITY_TAG "ZL"
#define BAD_REGION_TAG "ZR"
#define DUPLICATE_TAG "ZD"
#define ENZYME_TAG "ZE"

#define REV_BAM_FQCFAIL 1535
#define REV_BAM_FDUP 1023

void usage()
/* Explain usage and exit. */
{
errAbort(
  "bamtag - adds specific nonstandard Z-tags to the BAM file relevant to the\n"
  "   Beato Lab. All of the Z-tags will exclude the read (or read-pair) from\n"
  "   analysis.  The tags include:\n\n"
  "      ZL:Z   (string)        excluded because it is a B-read beyond specified\n"
  "                             cutoff value. (Also marks flag with 0x200 bit)\n"
  "      ZR:Z   (string)        excluded because it overlaps an HDR or other blacklist\n"
  "                             region (Also marks with 0x200 bit)\n"
  "      ZD:Z   (string)        exluded because it is an extraneous duplicate read.\n"
  "                             (Also marks flag with 0x400 bit)\n"
  "      ZE:Z   (string)        excluded because the read or its mate is not within the\n"
  "                             specified distance downstream of the read\'s start.\n"
  "                             (Also marks with 0x200 bit)\n"
  "usage:\n"
  "   bamtag TAGS input.bam output.bam\n"
  "where TAGS one or more of:\n"
  "   -clear                     clears 0x200 and 0x400 bits and ZL/ZD/ZR/ZE tags\n"
  "   -dupes                     *memory-intensive* and it's better to run this after\n"
  "                              merging whatever.\n"
  "   -add-MC                    add cigar string of the mate to the read (MC:Z tag)\n"
  "   -no-dupe-mark              when adding duplicates tags, don't mark\n"
  "   -B-read=percent            remove zero-quality or nearly such when a given\n"
  "                              percentage of the read (0.0-1.0] is specified\n"
  "                              (if percent unspecified, default is 0.4)\n"
  "                              (for paired-end, should be sorted by name)\n"
  "   -mapping-quality=q         minimum mapping quality\n"
  "   -same-chrom                mark different-chromosome pairings as QC-fail\n"
  "   -bad-bed=bed,name,percent  remove blacklist regions from BAM with an optional\n"
  "                              percent of the read overlapping one of the bad\n"
  "                              regions in question\n"
  "   -enzymes=distance,tag,bed1,bed2,...\n"
  "                              mark paired-end reads where one or both of the read\n"
  "                              starts is not located within the distance specified.\n"
  "In the case of paired-end reads, a read is tagged if its mate is also tagged.\n"
  "In the case that a tag already exists, it is appended or replaced with the new\n"
  "information.\n"
  );
}

static struct optionSpec options[] = 
{
    {"clear", OPTION_BOOLEAN},
    {"B-read", OPTION_DOUBLE},
    {"add-MC", OPTION_BOOLEAN},
    {"dupes", OPTION_BOOLEAN},
    {"no-dupe-mark", OPTION_BOOLEAN},
    {"bad-bed", OPTION_STRING},
    {"mapping-quality", OPTION_INT},
    {"good-sizes", OPTION_STRING},
    {"same-chrom", OPTION_BOOLEAN},
    {"report-only", OPTION_BOOLEAN},
    {"regions", OPTION_STRING},
    {"enzymes", OPTION_STRING},
    {NULL, 0},
};

void clear_tags(bam1_t *read1, bam1_t *read2, const char tag[2])
/* Just clear out the specific tag prior to writing a new one */
{
    uint8_t *tag1 = bam_aux_get(read1, tag);
    uint8_t *tag2 = NULL;
    if (tag1)
	bam_aux_del(read1, tag1);
    if (read2)
    {
	tag2 = bam_aux_get(read2, tag);
	if (tag2)
	    bam_aux_del(read2, tag2);
    }
}

boolean isa_B_read(bam1_t *bam, float perc)
/* check to see if a read has a high percent of phred quality 2 (##########) as the */
/* quality level */
{
    const bam1_core_t *core = &bam->core;
    int qlen = core->l_qseq;
    uint8_t *quals = bam_get_qual(bam);
    int i;
    int B = 0;
    for (i = 0; i < qlen; i++)
	if (quals[i] == 2)
	    B++;
    if (((float)B / (float)qlen) >= perc)
	return TRUE;
    return FALSE;
}

int mark_B_read(bam1_t *read1, bam1_t *read2, float perc)
/* pretty straight-forward, although I'm not in love with the design here. */
/* maybe could be a bit more elegant. */
{
    int ret = 0;
    boolean bad_read = FALSE;
    bam1_core_t *core1 = &read1->core;
    bam1_core_t *core2 = NULL;
    if (read2)
	core2 = &read2->core;
    char s[16];
    char *copy;
    safef(s, sizeof(s), "%0.2f", perc);
    copy = cloneString(s);
    if (!read1 && read2)
	uglyf("strangeness\n");
    clear_tags(read1, read2, ZERO_QUALITY_TAG);
    if (!read2)
    /* singleton case */
    {
	if (isa_B_read(read1, perc))
	    bad_read = TRUE;
	if (bad_read)
	{
	    bam_aux_append(read1, ZERO_QUALITY_TAG, 'Z', strlen(copy)+1, (uint8_t*)copy);	    
	    core1->flag |= BAM_FQCFAIL;
	    ret = 1;
	}
    }
    else
    /* paired end case */
    {
	if (isa_B_read(read1, perc) && isa_B_read(read2, perc))
	    bad_read = TRUE;
	if (bad_read)
	{
	    bam_aux_append(read1, ZERO_QUALITY_TAG, 'Z', strlen(copy)+1, (uint8_t*)copy);
	    bam_aux_append(read2, ZERO_QUALITY_TAG, 'Z', strlen(copy)+1, (uint8_t*)copy);
	    core1->flag |= BAM_FQCFAIL;
	    core2->flag |= BAM_FQCFAIL;
	    ret = 2;
	}
    }
    freeMem(copy);
    return ret;
}

boolean read_in_blacklist(const bam_hdr_t *header, bam1_t *read, struct genomeRangeTree *black_grt, float percent)
/* return true or false if a given % of bases of a read overlaps a bad region */
{
    if (read == NULL)
	return FALSE;
    const bam1_core_t *core = &read->core;
    if (core->tid < 0)
	return FALSE;
    int start, end;
    char *chrom = header->target_name[core->tid];
    int overlap_num = 0;
    start = core->pos;
    end = bam_endpos(read);
    if (!chrom)
	return FALSE;
    overlap_num = genomeRangeTreeOverlapSize(black_grt, chrom, start, end);
    if ((overlap_num > 0) && (((float)overlap_num / (float)(end - start)) >= percent))
	return TRUE;
    return FALSE;    
}

boolean read_upstream_enz(const bam_hdr_t *header, bam1_t *read, struct genomeRangeTree *enz_grt, unsigned enz_dist)
/* maybe this is too much like read_in_blacklist .. :( not great code reuse here */
{
    const bam1_core_t *core = &read->core;
    if (core->tid < 0)
	return FALSE;
    int start, end;
    char *chrom = header->target_name[core->tid];
    boolean rev_strand = FALSE;
    if (core->flag & BAM_FREVERSE)
	rev_strand = TRUE;
    start = core->pos;
    end = bam_endpos(read);
    if (rev_strand)
    {
	start -= enz_dist;
	if (start < 0)
	    start = 0;
    }
    else
    {
	end += enz_dist;
	if (end > header->target_len[core->tid])
	    end = header->target_len[core->tid];
    }

    if (!genomeRangeTreeOverlaps(enz_grt, chrom, start, end))
	return FALSE;
    return TRUE;    
}

int mark_bad_regions(const bam_hdr_t *header, bam1_t *read1, bam1_t *read2, 
		     struct genomeRangeTree *black_grt, float percent, char *bl_name, boolean just_one_bad)
/* check insert size and mark accordingly */
{
    int ret = 0;
    boolean bad = FALSE;
    bam1_core_t *core1 = &read1->core;
    bam1_core_t *core2 = NULL;
    if (read2)
	core2 = &read2->core;
    /* delete the old tags, if they exist */
    clear_tags(read1, read2, BAD_REGION_TAG);
    boolean bad_read1 = read_in_blacklist(header, read1, black_grt, percent);
    boolean bad_read2 = read_in_blacklist(header, read2, black_grt, percent);
    if ((!read2 && bad_read1) || (bad_read1 && bad_read2) || (just_one_bad && (bad_read1 || bad_read2)))
	bad = TRUE;
    if (read2 && bad)
    {
	bam_aux_append(read1, BAD_REGION_TAG, 'Z', strlen(bl_name)+1, (uint8_t*)bl_name);
	bam_aux_append(read2, BAD_REGION_TAG, 'Z', strlen(bl_name)+1, (uint8_t*)bl_name);
	core1->flag |= BAM_FQCFAIL;
	core2->flag |= BAM_FQCFAIL;
	ret = 2;
    }
    else if (bad)
    {
	bam_aux_append(read1, BAD_REGION_TAG, 'Z', strlen(bl_name)+1, (uint8_t*)bl_name);
	core1->flag |= BAM_FQCFAIL;
	ret = 1;
    }
    return ret;
}

int mark_no_enzymes(const bam_hdr_t *header, bam1_t *read1, bam1_t *read2, struct genomeRangeTree *enz_grt, 
		    unsigned enz_dist, char *enz_tag)
/* pretty straight-forward.  given a read, extend it and see if it overlaps in the GRT. */
{
    int ret = 0;
    boolean bad_read = FALSE;
    bam1_core_t *core1 = &read1->core;
    bam1_core_t *core2 = NULL;
    if (read2)
	core2 = &read2->core;
    /* delete the old tags, if they exist */
    clear_tags(read1, read2, ENZYME_TAG);
    if (!read2 || !read_upstream_enz(header, read1, enz_grt, enz_dist) || !read_upstream_enz(header, read2, enz_grt, enz_dist))
	bad_read = TRUE;
    if (read2 && bad_read)
    {
	bam_aux_append(read1, ENZYME_TAG, 'Z', strlen(enz_tag)+1, (uint8_t*)enz_tag);
	bam_aux_append(read2, ENZYME_TAG, 'Z', strlen(enz_tag)+1, (uint8_t*)enz_tag);
	core1->flag |= BAM_FQCFAIL;
	core2->flag |= BAM_FQCFAIL;
	ret = 2;
    }
    else if (bad_read)
    {
	bam_aux_append(read1, ENZYME_TAG, 'Z', strlen(enz_tag)+1, (uint8_t*)enz_tag);
	core1->flag |= BAM_FQCFAIL;
	ret = 1;
    }
    return ret;
}

void mark_reads(const bam_hdr_t *header, bam1_t *read1, bam1_t *read2, float B_read, 
		unsigned low_size, unsigned high_size, struct genomeRangeTree *black_grt, float black_percent, char *bl_name,
		boolean just_one_bad, struct genomeRangeTree *enz_grt, unsigned enz_dist, char *enz_tag, boolean same_chrom)
/* do all the different markings */
{
    bam1_core_t *core1 = &read1->core;
    bam1_core_t *core2 = (read2) ? &read2->core : NULL;
    if (B_read != -1)
	mark_B_read(read1, read2, B_read);
    /* only do this tagging on mapped reads */
    if (enz_grt && read2 && !((core1->flag & BAM_FUNMAP) || (core2->flag & BAM_FUNMAP)))
	mark_no_enzymes(header, read1, read2, enz_grt, enz_dist, enz_tag);
    else if (enz_grt && !read2 && !(core1->flag & BAM_FUNMAP))
	mark_no_enzymes(header, read1, NULL, enz_grt, enz_dist, enz_tag);
    if (black_grt && read2 && !((core1->flag & BAM_FUNMAP) || (core2->flag & BAM_FUNMAP)))
	mark_bad_regions(header, read1, read2, black_grt, black_percent, bl_name, just_one_bad);
    else if (black_grt && !read2 && !(core1->flag & BAM_FUNMAP))
	mark_bad_regions(header, read1, NULL, black_grt, black_percent, bl_name, just_one_bad);
    if (same_chrom)
    {
	if ((core1->tid != core1->mtid) && (core1->mtid > 0))
	    core1->flag |= BAM_FQCFAIL;
	if (core2 && ((core2->tid != core2->mtid) && (core2->mtid > 0)))
	    core2->flag |= BAM_FQCFAIL;
    }
}

void mark_quality(bam1_t *read, bam1_t *mate, unsigned mapping_quality)
/* reads and their mates under a certain mapping quality */
{
    if (read && mate)
    {
	bam1_core_t *read_core = &read->core;
	bam1_core_t *mate_core = &mate->core;
	if ((read_core->qual < mapping_quality) || (mate_core->qual < mapping_quality)) 
	{
	    read_core->flag |= BAM_FQCFAIL;
	    mate_core->flag |= BAM_FQCFAIL;
	}
    }
    if (read)
    {
	bam1_core_t *read_core = &read->core;
	if (read_core->qual < mapping_quality)
	    read_core->flag |= BAM_FQCFAIL;
    }
}

void fix_hic(bam1_t* read1, bam1_t *read2)
{
    bam1_core_t *core1 = &read1->core;
    bam1_core_t *core2 = &read2->core;
    core1->flag |= (BAM_FPAIRED | BAM_FREAD1);
    core2->flag |= (BAM_FPAIRED | BAM_FREAD2);
    if (!(core1->flag & BAM_FUNMAP) && !(core2->flag & BAM_FUNMAP))
    {
	core1->flag |= BAM_FPROPER_PAIR;
	core2->flag |= BAM_FPROPER_PAIR;
    }
}

char *cigar_string(bam1_t *b)
/* return the string representation of the cigar, with the mate's n_cigar size */
/* embedded, in case it's needed */
{
    int i;
    const bam1_core_t *c = &b->core;
    struct dyString *ds = dyStringNew(128);
    if (c->n_cigar) 
    {
	uint32_t *cigar = bam_get_cigar(b);
	for (i = 0; i < c->n_cigar; ++i) 
	{
	    dyStringPrintf(ds, "%d", bam_cigar_oplen(cigar[i]));
	    dyStringAppendC(ds, bam_cigar_opchr(cigar[i]));
	}
    }
    else
	dyStringAppendC(ds, '0');
    return dyStringCannibalize(&ds);
}

void add_mate_cigar(bam1_t *read1, bam1_t *read2)
/* add the mate's cigar string and vice versa */
{
    char *cigar1 = cigar_string(read1);
    char *cigar2 = cigar_string(read2);
    clear_tags(read1, read2, "MC");
    bam_aux_append(read1, "MC", 'Z', strlen(cigar2)+1, (uint8_t*)cigar2);
    bam_aux_append(read2, "MC", 'Z', strlen(cigar1)+1, (uint8_t*)cigar1);
    freeMem(cigar1);
    freeMem(cigar2);
}

void tag_the_bam(samFile *in, samFile *out, float B_read, unsigned low_size, 
		 unsigned high_size, struct genomeRangeTree *black_grt, float black_percent, char *bl_name,
                 boolean just_one_bad, struct genomeRangeTree *enz_grt, unsigned enz_dist, char *enz_tag, 
		 unsigned mapping_quality, boolean same_chrom, boolean add_mc)
/* Do the loop to do the checks */
{
    bam1_t *b[2];
    bam1_t *cur = NULL;
    bam1_t *pre = NULL;
    int curr = 0; 
    boolean has_prev = FALSE;
    bam_hdr_t *header = sam_hdr_read(in);
    sam_hdr_write(out, header);
    b[0] = bam_init1();
    b[1] = bam_init1();
    while (sam_read1(in, header, b[curr]) >= 0)
    {
	cur = b[curr];
	pre = b[1-curr];
	if (has_prev)
	{
	    char *name1 = bam_get_qname(cur);
	    char *name2 = bam_get_qname(pre);
	    if (name1 && name2 && sameString(name1, name2))
	    /* we have both in a pair of paired-end reads */
	    {
		if (enz_grt)
		    fix_hic(pre, cur);
		mark_reads(header, pre, cur, B_read, low_size, high_size, black_grt, black_percent, bl_name, just_one_bad, enz_grt, enz_dist, enz_tag, same_chrom);
		mark_quality(pre, cur, mapping_quality);
		if (add_mc)
		    add_mate_cigar(pre, cur);
		sam_write1(out, header, pre);
		sam_write1(out, header, cur);
		has_prev = FALSE;
	    }
	    else 
	    /* unpaired or single-end */
	    {
		mark_reads(header, pre, NULL, B_read, low_size, high_size, black_grt, black_percent, bl_name, just_one_bad, enz_grt, enz_dist, enz_tag, same_chrom);
		mark_quality(pre, NULL, mapping_quality);
		sam_write1(out, header, pre);
	    }
	}
	else 
	    has_prev = TRUE;
	curr = 1 - curr;
    }
    if (has_prev)
    /* lingering single-end read in odd number of reads */
    {
	mark_reads(header, pre, NULL, B_read, low_size, high_size, black_grt, black_percent, bl_name, just_one_bad, enz_grt, enz_dist, enz_tag, same_chrom);
	mark_quality(pre, NULL, mapping_quality);
	sam_write1(out, header, pre);
    }
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    bam_hdr_destroy(header);
}

void check_good_sizes(char *sizes_s, unsigned *low, unsigned *high)
/* check sizes option and return parsed values inside pointers */
{
    char *s = strchr(sizes_s, ':');
    unsigned l, h;
    if (s == NULL)
	errAbort("bad format for the sizes: should be low:high");
    *s++ = '\0';
    l = sqlUnsigned(sizes_s);
    h = sqlUnsigned(s);
    if (l >= h)
	errAbort("low should be lower than high size");
    *low = l;
    *high = h;
}

struct genomeRangeTree *load_black_grt(char *black_bed, char **p_bl_name, float *p_black_perc)
/* load bed file, then make a genomeRangeTree */
/* string put at **p_bl_name should be freed */
{
    struct genomeRangeTree *grt = genomeRangeTreeNew();
    char *words[3];
    int size = chopCommas(black_bed, words);
    char *bed_file = words[0];
    float black_perc = 0.4;
    char *name = NULL;
    if (size == 1)
	name = cloneString(bed_file);
    else if (size == 2)
	name = cloneString(words[1]);
    else if (size == 3)
    {
	name = cloneString(words[1]);
	black_perc = sqlFloat(words[2]);
    }
    struct bed *beds = bedLoadNAll(bed_file, 3);
    struct bed *bed;
    for (bed = beds; bed != NULL; bed = bed->next)
	genomeRangeTreeAdd(grt, bed->chrom, bed->chromStart, bed->chromEnd);
    bedFreeList(&beds);
    *p_bl_name = name;
    *p_black_perc = black_perc;
    return grt;
}

struct genomeRangeTree *load_enzyme_grt(char *params, unsigned *pDist, char **pEnz_tag)
/* load the bed of */
{
    struct slName *files = slNameListFromComma(params);
    struct slName *file;
    if (slCount(files) < 3)
	errAbort("expecting more parameters with -enzymes");
    struct slName *dist_name = slPopHead(&files);
    struct slName *enz_tag_sl = slPopHead(&files);
    *pEnz_tag = cloneString(enz_tag_sl->name);
    unsigned dist = sqlUnsigned(dist_name->name);
    struct genomeRangeTree *grt = genomeRangeTreeNew();;
    for (file = files; file != NULL; file = file->next)
    {
	struct bed *bedList = bedLoadNAll(file->name, 6);
	struct bed *bed;
	for (bed = bedList; bed != NULL; bed = bed->next)
	    genomeRangeTreeAdd(grt, bed->chrom, bed->chromStart, bed->chromEnd);
	bedFreeList(&bedList);
    }
    slNameFreeList(&files);
    slNameFree(&dist_name);
    *pDist = dist;
    return grt;
}

void flush_dupes_list(samFile *out, bam_hdr_t *header, struct hash *ch, struct slRef **pDupList, boolean nomarkdupe)
{
    struct slRef *ref;
    int dupCount = slCount(*pDupList);
    if (dupCount == 1)
    {
	/* it's not a duplicate, therefore it won't get tagged */
	ref = *pDupList;
	bam1_t *bam = (bam1_t *)ref->val;
	sam_write1(out, header, bam);
	bam_destroy1(bam);
	freez(pDupList);
    }
    else
    {
	/* it's a duplicate, so the hash is involved in the tagging */
	slReverse(pDupList);
	int i = 1;
	while ((ref = slPopHead(pDupList)) != NULL)
	{
	    bam1_t *bam = (bam1_t *)ref->val;
	    bam1_core_t *core = &bam->core;
	    char *name = bam_get_qname(bam);
	    int ix = hashIntValDefault(ch, name, 0);
	    char s[16];
	    char *cpy;
	    if (ix == 0)
	    {
		ix = i;
		hashAddInt(ch, name, ix);
	    }
	    else
		hashRemove(ch, name);
	    safef(s, sizeof(s), "%d.%d", dupCount, ix);
	    cpy = cloneString(s);
	    /* write out the bam */
	    bam_aux_append(bam, DUPLICATE_TAG, 'Z', strlen(cpy)+1, (uint8_t*)cpy);
	    if ((ix > 1) && !nomarkdupe)
		core->flag |= BAM_FDUP;
	    sam_write1(out, header, bam);
	    i++;
	    bam_destroy1(bam);
	    freez(&ref);
	}
    }
}

int bamstrandCmp(const void *a, const void *b)
/* for sorting after clustering */
{
    const struct slRef *ref_a = *((struct slRef **)a);
    const struct slRef *ref_b = *((struct slRef **)b);    
    const bam1_t *ba = ref_a->val;
    const bam1_t *bb = ref_b->val;
    const bam1_core_t *ca = &ba->core;
    const bam1_core_t *cb = &bb->core;
    int strand_a = (ca->flag & BAM_FREVERSE);
    int strand_b = (cb->flag & BAM_FREVERSE);
    return strand_a - strand_b;
}

void flush_dupes(samFile *out, bam_hdr_t *header, struct hash *ch, struct slRef **pDupList, boolean nomarkdupe)
/* separate into +/- strand first in case this isn't sorted already */
{
    struct slRef *ref;
    int dupCount = slCount(*pDupList);
    if (dupCount > 1)
    {
	struct slRef *pos_list = NULL;
	struct slRef *neg_list = NULL;
	struct slRef *cur = NULL;
	while ((cur = slPopHead(pDupList)) != NULL)
	{
	    bam1_t *cur_b = cur->val;
	    bam1_core_t *cur_core = &cur_b->core;
	    if (cur_core->flag & BAM_FREVERSE)
		slAddHead(&neg_list, cur);
	    else
		slAddHead(&pos_list, cur);
	}
	slReverse(&neg_list);
	slReverse(&pos_list);
	flush_dupes_list(out, header, ch, &pos_list, nomarkdupe);
	flush_dupes_list(out, header, ch, &neg_list, nomarkdupe);
	freez(&pos_list);
	freez(&neg_list);
	freez(pDupList);
    }
    else
	flush_dupes_list(out, header, ch, pDupList, nomarkdupe);
}

void dupetag(char *bigfile, samFile *out)
{
    struct slRef *dupList = NULL;
    struct hash *ch = newHash(28);
    samFile *in = sam_open(bigfile, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *b = bam_init1();
    boolean nomarkdupe = optionExists("no-dupe-mark");
    bam_hdr_write(out->fp.bgzf, header);
    while (bam_read1(in->fp.bgzf, b) >= 0)
    {
	bam1_core_t *cur_core = &b->core;
	cur_core->flag &= REV_BAM_FDUP;	
	clear_tags(b, NULL, DUPLICATE_TAG);
	if ((cur_core->flag & BAM_FUNMAP) || (cur_core->flag & BAM_FMUNMAP) || (cur_core->flag & BAM_FQCFAIL))
	    bam_write1(out->fp.bgzf, b);
	else if (dupList == NULL)
	{
	    /* very first thing read */
	    refAdd(&dupList, bam_dup1(b));
	}
	else
	{
	    bam1_t *list_bam = (bam1_t *)dupList->val;
	    bam1_core_t *list_core = &list_bam->core;
	    if ((list_core->tid != cur_core->tid) || (list_core->pos != cur_core->pos) || 
		(list_core->mtid != cur_core->mtid) || (list_core->mpos != cur_core->mpos))
	    {
		/* flush list */
		flush_dupes(out, header, ch, &dupList, nomarkdupe);
		refAdd(&dupList, bam_dup1(b));
	    }
	    else
	    {
		/* lengthen list */
		refAdd(&dupList, bam_dup1(b));
	    }
	}
    }
    if (dupList != NULL)
    {
	/* flush one more time */
	flush_dupes(out, header, ch, &dupList, nomarkdupe);
    }
    freeHash(&ch);
    hts_close(in);
    bam_destroy1(b);
    bam_hdr_destroy(header);
}

void cleartags(char *bigfile, samFile *out)
{
    samFile *in = sam_open(bigfile, "r");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *b = bam_init1();
    sam_hdr_write(out, header);
    while (sam_read1(in, header, b) >= 0)
    {
	bam1_core_t *cur_core = &b->core;
	clear_tags(b, NULL, DUPLICATE_TAG);
	clear_tags(b, NULL, BAD_REGION_TAG);
	clear_tags(b, NULL, DUPLICATE_TAG);
	clear_tags(b, NULL, ENZYME_TAG);
	cur_core->flag &= REV_BAM_FQCFAIL;
	cur_core->flag &= REV_BAM_FDUP;	
	sam_write1(out, header, b);
    }
    bam_hdr_destroy(header);
    sam_close(in);
    bam_destroy1(b);
}

void bamtag(char *bigfile, char *outputfile)
/* bamtag - main */
{
    samFile *out = NULL;
    samFile *in = NULL;
    if (optionExists("dupes"))
    {
	if (outputfile)
	    out = sam_open(outputfile, "wb");
	dupetag(bigfile, out);
	if (out)
	    sam_close(out);
    }
    else if (optionExists("clear"))
    {
	if (outputfile)
	    out = sam_open(outputfile, "wb");
	cleartags(bigfile, out);
	if (out)
	    sam_close(out);
    }
    else
    {
	float B_read = (float)optionDouble("B-read", 0.40);
	char *low_high = optionVal("good-sizes", NULL);
	unsigned low_size = 0;
	unsigned high_size = 0;
	bam_hdr_t* header;
	char *black_bed = optionVal("bad-bed", NULL);
	char *bl_name = NULL;
	boolean do_B_read = optionExists("B-read");
	boolean same_chrom = optionExists("same-chrom");
	boolean add_mc = optionExists("add-MC");
	float black_perc = 0;
	boolean just_one_bad = optionExists("bad-either");
	if (outputfile)
	    out = sam_open(outputfile, "wb");
	struct genomeRangeTree *black_grt = NULL;
	struct genomeRangeTree *enz_grt = NULL;
	char *enz_tag = NULL;
	unsigned enz_dist = 0;
	int mapping_quality = sqlUnsigned(optionVal("mapping-quality", "20"));
	if (optionExists("enzymes"))
	    enz_grt = load_enzyme_grt(optionVal("enzymes", NULL), &enz_dist, &enz_tag);
	if (black_bed)
	{
	    black_grt = load_black_grt(black_bed, &bl_name, &black_perc);
	}
	in = sam_open(bigfile, "r");
	/* first check to see at least one tagging is being done */
	if (!do_B_read)
	    B_read = -1;
	if (low_high)
	    check_good_sizes(low_high, &low_size, &high_size);
	tag_the_bam(in, out, B_read, low_size, high_size, black_grt, black_perc, bl_name, just_one_bad, enz_grt, enz_dist, enz_tag, mapping_quality, same_chrom, add_mc);
	if (black_grt)
	    genomeRangeTreeFree(&black_grt);
	if (enz_tag)
	    freeMem(enz_tag);
	if (enz_grt)
	    genomeRangeTreeFree(&enz_grt);
	sam_close(in);
	if (bl_name)
	    freeMem(bl_name);
	if (out)
	    sam_close(out);
    }
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if ((argc != 3) && (argc != 2))
    usage();
if (argc == 2)
    bamtag(argv[1], NULL);
else 
    bamtag(argv[1], argv[2]);
return 0;
}
