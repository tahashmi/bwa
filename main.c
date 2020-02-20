#include <arrow-glib/arrow-glib.h>
#include <plasma-glib/plasma-glib.h>

#include <stdint-gcc.h>

#include <stdio.h>
#include <string.h>
#include "kstring.h"
#include "utils.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.17-r1198-dirty"
#endif

#define ASCII_START 32
#define ASCII_END 126

typedef char   gchar;
extern int64_t counts1,counts2,counts3,counts4,counts5,counts5,counts6,counts7,counts8,counts9,counts10,counts11,
        counts12,counts13,counts14,counts15,counts15,counts16,counts17,counts18,counts19,counts20,counts21,counts22,
        countsX,countsY,countsM;

int bwa_fa2pac(int argc, char *argv[]);
int bwa_pac2bwt(int argc, char *argv[]);
int bwa_bwtupdate(int argc, char *argv[]);
int bwa_bwt2sa(int argc, char *argv[]);
int bwa_index(int argc, char *argv[]);
int bwt_bwtgen_main(int argc, char *argv[]);

int bwa_aln(int argc, char *argv[]);
int bwa_sai2sam_se(int argc, char *argv[]);
int bwa_sai2sam_pe(int argc, char *argv[]);

int bwa_bwtsw2(int argc, char *argv[]);

int main_fastmap(int argc, char *argv[]);
int main_mem(int argc, char *argv[]);
int main_shm(int argc, char *argv[]);

int main_pemerge(int argc, char *argv[]);
int main_maxk(int argc, char *argv[]);
	
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bwa (alignment via Burrows-Wheeler transformation)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         mem           BWA-MEM algorithm\n");
	fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
	fprintf(stderr, "         pemerge       merge overlapping paired ends (EXPERIMENTAL)\n");
	fprintf(stderr, "         aln           gapped/ungapped alignment\n");
	fprintf(stderr, "         samse         generate alignment (single ended)\n");
	fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
	fprintf(stderr, "         bwasw         BWA-SW for long queries\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         shm           manage indices in shared memory\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
	fprintf(stderr, "\n");
	fprintf(stderr,
"Note: To use BWA, you need to first index the genome with `bwa index'.\n"
"      There are three alignment algorithms in BWA: `mem', `bwasw', and\n"
"      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'\n"
"      first. Please `man ./bwa.1' for the manual.\n\n");
	return 1;
}

/////////////////////////////////////////////////////
//		Arrow-Functionality                //
/////////////////////////////////////////////////////
char* generateRandomString(int size) {
    int i;
    char *res = malloc(size + 1);
    for(i = 0; i < size; i++) {
        res[i] = (char) (rand()%(ASCII_END-ASCII_START))+ASCII_START;
    }
    res[size + 1] = '\0';
    return res;
}

void create_plasma_object(GArrowRecordBatch *batch_genomics)
{
    guint8 id_arr[20];
    char objID_file[] = "objID.txt";
    memcpy(id_arr, generateRandomString(20),20);

    FILE *fOut;
    fOut = fopen(objID_file, "a");
    fputs(id_arr, fOut);
    fputc('\n', fOut);
    fclose(fOut);

    gboolean success = TRUE;
    GError *error = NULL;

    GPlasmaClient *gPlasmaClient;
    GPlasmaObjectID *object_id;
    GPlasmaClientCreateOptions *create_options;
    GPlasmaClientOptions *gplasmaClient_options;
    GPlasmaCreatedObject *Object;
    GArrowBuffer *arrowBuffer;

    arrowBuffer = GSerializeRecordBatch(batch_genomics);
    gint64 size = garrow_buffer_get_size(arrowBuffer);

    //g_print("obj_id: %s\n", id_arr);
    fprintf(stderr, "[%s] obj_id: %s , size: %ld \n", __func__, id_arr, size);

    create_options = gplasma_client_create_options_new();
    gplasmaClient_options = gplasma_client_options_new();
    gPlasmaClient = gplasma_client_new("/tmp/store0",gplasmaClient_options, &error);
    object_id = gplasma_object_id_new(id_arr, 20, &error);

    {
        /* It should be guint8 instead of gchar. We use gchar here
           just for convenient. */
        guint8 metadata[] = "metadata";
        gplasma_client_create_options_set_metadata(create_options, (const guint8 *)metadata, sizeof(metadata));
    }
    Object = gplasma_client_create(gPlasmaClient, object_id, size, create_options, &error);

    g_object_unref(create_options);
    {
        GArrowBuffer *data;
        g_object_get(Object, "data", &data, NULL);
        //garrow_mutable_buffer_set_data(GARROW_MUTABLE_BUFFER(data),0,arrowBuffer,&error);
        garrow_mutable_buffer_set_data(GARROW_MUTABLE_BUFFER(data),0, garrow_buffer_get_databytes(arrowBuffer),size,&error);
        g_object_unref(data);
    }

    gplasma_created_object_seal(Object, &error);
    g_object_unref(Object);
    gplasma_client_disconnect(gPlasmaClient, &error);
    g_object_unref(gPlasmaClient);
    g_object_unref(object_id);
}

void arrow_plasma_create(int id, long size) {
    GArrowRecordBatch * rb_genomics;
    rb_genomics=arrow_builders_finish(id, size);
    {
        //g_print("%s", garrow_record_batch_to_string(rb_genomics, NULL));
        //fprintf(stderr, "[%s] Arrow: %s\n", __func__, garrow_record_batch_to_string(rb_genomics, NULL));
    }
    create_plasma_object(rb_genomics);
    g_object_unref(rb_genomics);
}
int exist(const char *name)
{
    struct stat   buffer;
    return (stat (name, &buffer) == 0);
}
////////////////////////END////////////////////////////

int main(int argc, char *argv[])
{
	extern char *bwa_pg;
	int i, ret;
	double t_real;
	kstring_t pg = {0,0,0};
	
	char file[] = "objID.txt";
        if(exist(file))
          remove(file);
        arrow_builders_start();
	
	t_real = realtime();
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "fa2pac") == 0) ret = bwa_fa2pac(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwt") == 0) ret = bwa_pac2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwtgen") == 0) ret = bwt_bwtgen_main(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtupdate") == 0) ret = bwa_bwtupdate(argc-1, argv+1);
	else if (strcmp(argv[1], "bwt2sa") == 0) ret = bwa_bwt2sa(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) ret = bwa_index(argc-1, argv+1);
	else if (strcmp(argv[1], "aln") == 0) ret = bwa_aln(argc-1, argv+1);
	else if (strcmp(argv[1], "samse") == 0) ret = bwa_sai2sam_se(argc-1, argv+1);
	else if (strcmp(argv[1], "sampe") == 0) ret = bwa_sai2sam_pe(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtsw2") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "dbwtsw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "bwasw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "fastmap") == 0) ret = main_fastmap(argc-1, argv+1);
	else if (strcmp(argv[1], "mem") == 0) ret = main_mem(argc-1, argv+1);
	else if (strcmp(argv[1], "shm") == 0) ret = main_shm(argc-1, argv+1);
	else if (strcmp(argv[1], "pemerge") == 0) ret = main_pemerge(argc-1, argv+1);
	else if (strcmp(argv[1], "maxk") == 0) ret = main_maxk(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}
	free(bwa_pg);
	
	arrow_plasma_create(1,counts1); //Chrms
	arrow_plasma_create(2,counts2);
	arrow_plasma_create(3,counts3);
	arrow_plasma_create(4,counts4);
	arrow_plasma_create(5,counts5);
	arrow_plasma_create(6,counts6);
	arrow_plasma_create(7,counts7);
	arrow_plasma_create(8,counts8);
	arrow_plasma_create(9,counts9);
	arrow_plasma_create(10,counts10);
	arrow_plasma_create(11,counts11);
	arrow_plasma_create(12,counts12);
	arrow_plasma_create(13,counts13);
	arrow_plasma_create(14,counts14);
	arrow_plasma_create(15,counts15);
	arrow_plasma_create(16,counts16);
	arrow_plasma_create(17,counts17);
	arrow_plasma_create(18,counts18);
	arrow_plasma_create(19,counts19);
	arrow_plasma_create(20,counts20);
	arrow_plasma_create(21,counts21);
	arrow_plasma_create(22,counts22);
	arrow_plasma_create(23,countsX);
	arrow_plasma_create(24,countsY);
	arrow_plasma_create(25,countsM);
	
	return ret;
}
