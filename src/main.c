// Copyright (c) 2013 Genome Research Limited.
//
// This file is part of samtools addreplacereadgroups.
//
// samtools addreplacereadgroups  is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see L<http://www.gnu.org/licenses/>.

#define _GNU_SOURCE

#include <htslib/sam.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>

// BEGIN: this bit should be in htslib
#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

typedef khash_t(s2i) sdict_t;

//#define CLONE(name) \
//extern kh_##name##_t *kh_clone_##name(kh_##name##_t *h);
//CLONE(s2i)

#define kh_clone(name,h) kh_clone_##name(h)


#define __ac_fsize(m) ((m) < 16? 1 : (m)>>4)


#define KHASH_INIT_CLONE(name, SCOPE, khkey_t, khval_t) \
SCOPE kh_##name##_t* kh_clone_##name(kh_##name##_t *h)						\
{ \
    if (h) { \
        kh_##name##_t* r = calloc(1, sizeof(kh_##name##_t));        \
        r->n_buckets = h->n_buckets;                                \
        r->size = h->size;                                          \
        r->n_occupied = h->n_occupied;                              \
        r->upper_bound = h->upper_bound;                            \
        r->keys = calloc(h->n_buckets, sizeof(khkey_t));			\
        memcpy(r->keys, h->keys, h->n_buckets * sizeof(khkey_t) );  \
        h->flags = calloc(h->n_buckets, sizeof(khkey_t));			\
        memcpy(r->flags, h->flags, __ac_fsize(h->n_buckets) * sizeof(khkey_t)); \
        h->vals = calloc(h->n_buckets, sizeof(khval_t));			\
        memcpy(r->vals, h->vals, h->n_buckets * sizeof(khval_t) );  \
        return r;                                                   \
    }                                                               \
    return NULL;                                                    \
}

KHASH_INIT_CLONE(s2i, static, kh_cstr_t, int64_t)
// END: this bit should be in htslib
// BEGIN: this bit should be in htslib sam.c

bam_hdr_t* dup_header( const bam_hdr_t* input_header ) {
    
    if (input_header == NULL) return NULL;
    bam_hdr_t* retval = (bam_hdr_t*)malloc(sizeof(bam_hdr_t));
    // TODO: Need to clone interior of this?
    bzero((void*)retval, sizeof(bam_hdr_t));
    retval->n_targets = input_header->n_targets;
    retval->l_text = input_header->l_text;
    retval->ignore_sam_err = input_header->ignore_sam_err;
    
    int32_t i;
	if (input_header->target_name) {
        retval->target_len = calloc(input_header->n_targets, sizeof(uint32_t));
        memcpy(retval->target_len, input_header->target_len, sizeof(uint32_t)*input_header->n_targets);
        retval->target_name = calloc(input_header->n_targets, sizeof(char*));
		for (i = 0; i < input_header->n_targets; ++i)
			retval->target_name[i] = strdup(input_header->target_name[i]);
	}
	retval->text = strdup(input_header->text);
    if (input_header-> cigar_tab) {
        retval->cigar_tab = (uint8_t*)calloc(128, sizeof(uint8_t));
        memcpy(retval->cigar_tab, input_header->cigar_tab, 128 * sizeof(uint8_t));
    }
	if (input_header->sdict)
    {
        retval->sdict = (void*)kh_clone(s2i, (sdict_t*)input_header->sdict);
    }
    
    return retval;
}
// END: this bit should be in htslib sam.c

struct parsed_opts {
    char* input_name;
    char* output_name;
};

struct state {
    samFile* input_file;
    bam_hdr_t* input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
};

typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 2) {
        dprintf(STDERR_FILENO, "Arguments should be: readgroupise <input.bam> <output.bam>\r\n");
        return NULL;
    }
    
    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) return NULL;

    retval->input_name = strdup(argv[1]);

    retval->output_name = strdup(argv[2]);

    return retval;
}

bool init(parsed_opts_t* opts, state_t** state_out) {
    state_t* retval = (state_t*) malloc(sizeof(state_t));
    if (retval == NULL) {
        dprintf(STDERR_FILENO, "Out of memory");
        return false;
    }
    *state_out = retval;
    
    // Open files
    retval->input_file = sam_open(opts->input_name, "rb", 0);
    if (retval->input_file == NULL) {
        dprintf(STDERR_FILENO, "Could not open input file: %s\r\n", opts->input_name);
        return false;
    }
    retval->input_header = sam_hdr_read(retval->input_file);
    
    retval->output_header = dup_header(retval->input_header);
    retval->output_file = sam_open(opts->output_name, "wb", 0);
    
    if (retval->output_file == NULL) {
        dprintf(STDERR_FILENO, "Could not open output file: %s\r\n", opts->output_name);
        return false;
    }
    
    return true;
}

char* getRGLine(const char* text) {
    char* rg = strstr(text,"\n@RG");
    if (rg == NULL) {
        return NULL;
    }
    rg++;//skip initial \n
    char* end = strchr(rg, '\n');
    char* line;
    if (end)
    {
        line = strndup(rg,(end-rg));
    } else {
        line = strdup(rg);
    }
    
    return line;
}
char* getRGID(const char* text) {
    char *line, *next;
    line = getRGLine(text);
 
    assert(line!=NULL);
    
    next = line;
    char* token = strsep(&next, "\t");
    token = strsep(&next,"\t"); // skip first token it should always be "@RG"
    while (next != NULL) {
        char* key = strsep(&token,":");
        if (!strcmp(key,"ID")) {
            char* retval = strdup(token);
            free(line);
            return retval;
        }
        token = strsep(&next,"\t");
    }
    free (line);
}

bool readgroupise(state_t* state) {
    char* id = getRGID(state->output_header->text);

    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        free(id);
        return false;
    }
    
    bam1_t* file_read = bam_init1();
    if (sam_read1(state->input_file, state->input_header, file_read) < 0) {
        bam_destroy1(file_read);
        file_read = NULL;
    }
    while (file_read != NULL) {
        uint8_t* data = (uint8_t*)strdup(id);
        int len = strlen(id)+1;
        // If the old exists delete it
        uint8_t* old = bam_aux_get(file_read, "RG");
        if (old != NULL) {
            bam_aux_del(file_read, old);
        }
        bam_aux_append(file_read, "RG",'Z',len,data);
        sam_write1(state->output_file, state->output_header, file_read);
        if (sam_read1(state->input_file, state->input_header, file_read) < 0) {
            bam_destroy1(file_read);
            file_read = NULL;
        }
    }
    free(id);

    return true;
}

void cleanup_opts(parsed_opts_t* opts) {
    free(opts->output_name);
    free(opts->input_name);
}

void cleanup_state(state_t* state) {
    sam_close(state->output_file);
    bam_hdr_destroy(state->output_header);
    sam_close(state->input_file);
    bam_hdr_destroy(state->input_header);
}

int main(int argc, char** argv) {

    parsed_opts_t* opts = parse_args(argc, argv);
    state_t* state = NULL;
    if (!opts || !init(opts, &state)) return -1;
    
    if (!readgroupise(state)) return -1;
    
    cleanup_opts(opts);
    cleanup_state(state);
    
    return 0;
}
