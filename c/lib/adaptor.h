#ifndef __ADAPTOR_H__
#define __ADAPTOR_H__

#include "seq.h"
#include "aln.h"

Seq *adaptor_read_seq(const char *filename);

void adaptor_mask_clear(Seq *seq, unsigned char *mask);

void adaptor_mask_left(Seq *seq, unsigned char *mask,
		       AlnNode *aln_end);

void adaptor_mask_right(Seq *seq, unsigned char *mask,
			AlnNode *aln_end);
#endif
