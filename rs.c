/*
	This program handles CODAR Range Series (RS) files, as ncdump and ncgen do for NetCDF data.
	The same executable can be called rsdump or rsgen to act as follows:
	- rsdump reads a binary RS file and generates an ASCII text representation of the data that can then be edited.
	- rsgen reads an ascii file produced by rsdump and converts it into a binary RS file.

	(c) 2021 Marcel Losekoot, Bodega Marine Laboratory, UC Davis.
	Based on ts.c, added Debug, added fprintf for error messages, added hexdump for undocumented blocks.
	Bugs: only deals with fbin settings cviq and flt4, data type for iqdata is hardcoded to float
	Doc Bugs: block sign.nOwner is really sitecode, undefined block hasi

	Notes: the binary RS file is bigendian by definition, so the program tests itself and corrects accordingly.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>		// ctime()
#include <ctype.h>		// islower(), isprint()
#include <unistd.h>		// fseek constant SEEK_END
#include <string.h>		// memset()
#include <stdint.h>		// uint32_t
#include <libgen.h>		// basename()
#include <math.h>		// round()


char Version[] = "rs.c version 1.0a 2021-02-15" ;
int Debug = 0 ;			// Global debug flag to produce verbose output to stderr

#define SIZE_LINE 256		// Maximum length of a line of ascii text


typedef uint32_t fourcc ;	// four bytes that are subject to byte swapping

// declare a node for the linked list
struct node
{
	fourcc key ;		// the key name
	uint32_t size ;		// the size of the node's data block
	unsigned char *data ;	// the data block
	struct node *next ;	// the next node
} ;

// define key value codes, used to recognize and label block types
#define KEY_AQFT	(fourcc )0x41514654	// "AQFT"
#define KEY_HEAD	(fourcc )0x48454144	// "HEAD"
#define KEY_sign	(fourcc )0x7369676e	// "sign"
#define KEY_mcda	(fourcc )0x6d636461	// "mcda"
#define KEY_dbrf	(fourcc )0x64627266	// "dbrf"
#define KEY_cnst	(fourcc )0x636e7374	// "cnst"
#define KEY_hasi	(fourcc )0x68617369	// "hasi"
#define KEY_swep	(fourcc )0x73776570	// "swep"
#define KEY_fbin	(fourcc )0x6662696e	// "fbin"
#define KEY_BODY	(fourcc )0x424f4459	// "BODY"
#define KEY_rtag	(fourcc )0x72746167	// "rtag"
#define KEY_gps1	(fourcc )0x67707331	// "gps1"
#define KEY_indx	(fourcc )0x696e6478	// "indx"
#define KEY_scal	(fourcc )0x7363616c	// "scal"
#define KEY_afft	(fourcc )0x61666674	// "afft"
#define KEY_ifft	(fourcc )0x69666674	// "ifft"
#define KEY_END 	(fourcc )0x454e4420	// "END "

// define size constants for character strings
#define SIZE_DESCRIPTION	64
#define SIZE_OWNERNAME		64
#define SIZE_COMMENT		64

// define codes for recognizing binary format and type options
#define BINFORMAT_CVIQ	(fourcc )0x63766971	// "cviq"
#define BINFORMAT_DBRA	(fourcc )0x64627261	// "dbra"
#define BINTYPE_FLT8	(fourcc )0x666c7438	// "flt8"
#define BINTYPE_FLT4	(fourcc )0x666c7434	// "flt4"
#define BINTYPE_FIX2	(fourcc )0x66697832	// "fix2"
#define BINTYPE_FIX3	(fourcc )0x66697833	// "fix3"
#define BINTYPE_FIX4	(fourcc )0x66697834	// "fix4"

// define a struct to hold values that are referenced across blocks
struct config
{
	fourcc bin_format ;			// the binary format: 'cviq' or 'dbra'
	fourcc bin_type ;			// the type of binary data: one of 'flt8', 'flt4', 'fix2', 'fix3' or 'fix4'
	int32_t index ; 			// index (0 to dopplercells-1)
	double scalar_one ;			// scaling value for I
	double scalar_two ;			// scaling value for Q
} ;

struct block_header
{
	fourcc key ;
	uint32_t size ;
} __attribute__((packed)) ;	// disable padding to make the header line up with the file data

struct block_functions						// this struct is used to relate a key name with a set of functions
{
	fourcc key ;						// a 4 byte block key
	int (*fixup)(struct node *) ;				// a pointer to a function that is called to perform endian fixup on the data block
	int (*make)(struct node *, struct config *, FILE *) ;	// a pointer to a function that is called to create a data block from text
	int (*dump)(struct node *, struct config *, FILE *) ;	// a pointer to a function that is called to produce text output from a data block
	int (*gen)(struct node *, FILE *) ;			// a pointer to a function that is called to write out a binary version of the block
} ;


int check_little_endian(void) ;
void usage_rsdump(char *) ;
void usage_rsgen(char *) ;
int rsdump(FILE *, FILE *, int) ;
int rsgen(FILE *, FILE *) ;
int read_binary_file(FILE *, unsigned long, unsigned char *) ;
int check_header(unsigned char *) ;
struct node *parse_file(unsigned char *, unsigned long) ;
int parse_block(struct node *, unsigned char *, unsigned long) ;
int superblock(fourcc) ;
void show_list(struct node *) ;
void free_all_nodes(struct node *) ;
void free_all_nodes_and_data(struct node *) ;
void debugdump(unsigned char *, int) ;
void hexdump(unsigned char *, unsigned int, FILE *) ;
void endian_fixup(void *, int) ;
void swapcopy(unsigned char *, unsigned char *, int) ;
void swapcopy2(unsigned char *, unsigned char *) ;
void swapcopy4(unsigned char *, unsigned char *) ;
void swapcopy8(unsigned char *, unsigned char *) ;
char *strkey(fourcc) ;
int fixup_sizes(struct node *) ;
uint32_t calculate_body_size(struct node *) ;
uint32_t calculate_head_size(struct node *) ;
int set_block_size(struct node *, fourcc , uint32_t) ;
int dump_list(struct node *, FILE *, int) ;
void chomp(char *, int) ;
int read_parameter(FILE *, char [], void *) ;
int fixup_block(struct node *) ;
struct block_functions *find_block_functions(fourcc) ;
int rs_write(struct node *, FILE *) ;
int count_iqdata_lines(FILE *) ;

// a set of functions that dump the contents of a specific type of block
int dump_block_aqft(struct node *, struct config *, FILE *) ;
int dump_block_head(struct node *, struct config *, FILE *) ;
int dump_block_sign(struct node *, struct config *, FILE *) ;
int dump_block_mcda(struct node *, struct config *, FILE *) ;
int dump_block_dbrf(struct node *, struct config *, FILE *) ;
int dump_block_cnst(struct node *, struct config *, FILE *) ;
int dump_block_hasi(struct node *, struct config *, FILE *) ;
int dump_block_swep(struct node *, struct config *, FILE *) ;
int dump_block_fbin(struct node *, struct config *, FILE *) ;
int dump_block_body(struct node *, struct config *, FILE *) ;
int dump_block_rtag(struct node *, struct config *, FILE *) ;
int dump_block_gps1(struct node *, struct config *, FILE *) ;
int dump_block_indx(struct node *, struct config *, FILE *) ;
int dump_block_scal(struct node *, struct config *, FILE *) ;
int dump_block_afft(struct node *, struct config *, FILE *) ;
int dump_block_ifft(struct node *, struct config *, FILE *) ;
int dump_block_end(struct node *, struct config *, FILE *) ;

// a set of functions that fixup the data portion of a specific type of block
int fixup_data_aqft(struct node *) ;
int fixup_data_head(struct node *) ;
int fixup_data_sign(struct node *) ;
int fixup_data_mcda(struct node *) ;
int fixup_data_dbrf(struct node *) ;
int fixup_data_cnst(struct node *) ;
int fixup_data_hasi(struct node *) ;
int fixup_data_swep(struct node *) ;
int fixup_data_fbin(struct node *) ;
int fixup_data_body(struct node *) ;
int fixup_data_rtag(struct node *) ;
int fixup_data_gps1(struct node *) ;
int fixup_data_indx(struct node *) ;
int fixup_data_scal(struct node *) ;
int fixup_data_afft(struct node *) ;
int fixup_data_ifft(struct node *) ;
int fixup_data_end(struct node *) ;

// a set of functions that create a node for a specific type of block
int make_node_aqft(struct node *, struct config *config, FILE *) ;
int make_node_head(struct node *, struct config *config, FILE *) ;
int make_node_sign(struct node *, struct config *config, FILE *) ;
int make_node_mcda(struct node *, struct config *config, FILE *) ;
int make_node_dbrf(struct node *, struct config *config, FILE *) ;
int make_node_cnst(struct node *, struct config *config, FILE *) ;
int make_node_hasi(struct node *, struct config *config, FILE *) ;
int make_node_swep(struct node *, struct config *config, FILE *) ;
int make_node_fbin(struct node *, struct config *config, FILE *) ;
int make_node_body(struct node *, struct config *config, FILE *) ;
int make_node_rtag(struct node *, struct config *config, FILE *) ;
int make_node_gps1(struct node *, struct config *config, FILE *) ;
int make_node_indx(struct node *, struct config *config, FILE *) ;
int make_node_scal(struct node *, struct config *config, FILE *) ;
int make_node_afft(struct node *, struct config *config, FILE *) ;
int make_node_ifft(struct node *, struct config *config, FILE *) ;
int make_node_end(struct node *, struct config *config, FILE *) ;

// a set of functions that generate binary file data for a specific type of block
int gen_block_aqft(struct node *, FILE *) ;
int gen_block_head(struct node *, FILE *) ;
int gen_block_sign(struct node *, FILE *) ;
int gen_block_mcda(struct node *, FILE *) ;
int gen_block_dbrf(struct node *, FILE *) ;
int gen_block_cnst(struct node *, FILE *) ;
int gen_block_hasi(struct node *, FILE *) ;
int gen_block_swep(struct node *, FILE *) ;
int gen_block_fbin(struct node *, FILE *) ;
int gen_block_body(struct node *, FILE *) ;
int gen_block_rtag(struct node *, FILE *) ;
int gen_block_gps1(struct node *, FILE *) ;
int gen_block_indx(struct node *, FILE *) ;
int gen_block_scal(struct node *, FILE *) ;
int gen_block_afft(struct node *, FILE *) ;
int gen_block_ifft(struct node *, FILE *) ;
int gen_block_end(struct node *, FILE *) ;


int Global_flag_little_endian = 1 ;	// 1 indicates this code is little endian, 0 means it's big endian. The binary file is big endian.


int main(int argc, char *argv[])
{
	char *program_name = basename(argv[0]) ;
	Global_flag_little_endian = check_little_endian() ;
	int err = 0 ;
	FILE *fdin ;
	FILE *fdout ;
	if( strcmp(program_name,"rsdump") == 0 )		// the program name must be rsdump or rsgen
	{
		// do rsdump
		if( argc < 2 )
		{
			usage_rsdump(program_name) ;
			return 0 ;
		}
		int just_header = 0 ;
		if( strcmp(argv[1],"-h") == 0 )
		{
			just_header = 1 ;
			argv++ ;
			argc-- ;
		}
		char *infilename = argv[1] ;
		if( (fdin = fopen(infilename,"rb")) == NULL )
		{
			fprintf(stderr,"Cannot open input file '%s'\n",infilename) ;
			return 1 ;
		}
		if( argc > 2 )	// an outfile is supplied
		{
			char *outfilename = argv[2] ;
			if( (fdout = fopen(outfilename,"wt")) == NULL )
			{
				fprintf(stderr,"Cannot open output file '%s'\n",outfilename) ;
				fclose(fdin) ;
				return 1 ;
			}
		}
		else
		{
			fdout = stdout ;
		}
		err = rsdump(fdin,fdout,just_header) ;
	}
	if( strcmp(program_name,"rsgen") == 0 )
	{
		// do rsgen
		if( argc < 3 )
		{
			usage_rsgen(program_name) ;
			return 0 ;
		}
		char *infilename = argv[1] ;
		if( (fdin = fopen(infilename,"rt")) == NULL )
		{
			fprintf(stderr,"Cannot open input file '%s'\n",infilename) ;
			return 1 ;
		}
		char *outfilename = argv[2] ;
		if( (fdout = fopen(outfilename,"wb")) == NULL )
		{
			fprintf(stderr,"Cannot open output file '%s'\n",outfilename) ;
			if( fdout != stdout )
				fclose(fdout) ;
			return 1 ;
		}
		err = rsgen(fdin,fdout) ;
	}
	fclose(fdin) ;
	if( fdout != stdout )
		fclose(fdout) ;
	return err ;
}

int check_little_endian(void)	// check whether this code is big or little endian, returns 1 for little.
{
	unsigned int i = 1 ;
	char *c = (char *)&i ;
	if( *c )		
		return 1 ;	// little endian
	return 0 ;
}

void usage_rsdump(char *name)
{
	fprintf(stderr,"Usage: %s [-h] infile [outfile]\n",name) ;
	fprintf(stderr,"Processes CODAR SeaSonde RangeSeries data files.\n") ;
	fprintf(stderr,"%s\n",Version) ;
}

void usage_rsgen(char *name)
{
	fprintf(stderr,"Usage: %s infile outfile\n",name) ;
	fprintf(stderr,"Processes CODAR SeaSonde RangeSeries data files.\n") ;
	fprintf(stderr,"Reads an ascii text infile and writes a binary version to outfile.\n") ;
	fprintf(stderr,"%s\n",Version) ;
}

int rsdump(FILE *infile, FILE *outfile, int just_header)	// top level function in rsdump mode
// read the binary file into a buffer
// parse the buffer for RIFF blocks
// make a linked list of nodes
// write a description for each node to a text file
{
	fseek(infile,0L,SEEK_END) ;	// seek to the end of the file
	unsigned long filesize = ftell(infile) ;
	rewind(infile) ;
	unsigned char *filedata = malloc(filesize) ;	// try to make a buffer sized to read the entire file
	if( filedata == NULL )
	{
		fprintf(stderr,"Cannot get memory for file with %ld bytes\n",filesize) ;
		return 1 ;
	}
	int err = 0 ;
	if( read_binary_file(infile,filesize,filedata) == 0 )
	{
		struct node *list = parse_file(filedata,filesize) ;
		if( list != NULL )
		{
			err = dump_list(list,outfile,just_header) ;
			free_all_nodes(list) ;
		}
	}
	free(filedata) ;
	return err ;
}

int rsgen(FILE *infile, FILE *outfile)	// top level function in rsgen mode
// read lines of text from a text file
// parse the block key names
// call the relevant make function to read related data from the text file and make a linked list node for the block
// write the linked list to a binary RS file
{
	char line[SIZE_LINE] ;
	long line_count = 0 ;
	struct config config ;
	memset(&config,0,sizeof(struct config)) ;
	struct node root ;
	struct node *list = &root ;
	memset(list,0,sizeof(struct node)) ;
	while( fgets(line,SIZE_LINE,infile) )
	{
		chomp(line,SIZE_LINE) ;				// remove newline
		if( Debug ) { fprintf(stderr,"debug: chomped line is '%s'\n",line) ; }
		line_count++ ;
		if( strlen(line) <= 1 ) continue ;		// skip empty lines
		if( index(line,':') != NULL ) continue ;	// skip parameter lines
		fourcc key = *(fourcc *)line ;			// extract the block type
		endian_fixup(&key,sizeof(key)) ;
		struct block_functions *block_functions = find_block_functions(key) ;	// returns a set of functions from the Global_functions_dictionary for this block type
		if( block_functions == NULL )
		{
			fprintf(stderr,"Cannot gen block '%s'\n",strkey(key)) ;
			return 1 ;
		}
		int (*make_function)(struct node *, struct config *, FILE *) = block_functions->make ;
		int err = (*make_function)(list,&config,infile) ;	// calls the 'make' function from Function_dictionary corresponding to the block type
		if( err )
		{
			fprintf(stderr,"Error in '%s' block starting at line %ld\n",strkey(key),line_count) ;
			return 1 ;
		}
		if( list->next != NULL ) list = list->next ;	// advance the list pointer to the newly created node
	}
	printf("Read %ld lines\n",line_count) ;
	fixup_sizes(&root) ;	// calculate body, head and aqft block sizes, update nodes
	// write to outfile
	int err = rs_write(root.next,outfile) ;
	free_all_nodes_and_data(list) ;
	return err ;
}

int rs_write(struct node *list, FILE *outfile)		// writes a binary RS file using the data from the linked list
{
	if( Debug ) { fprintf(stderr,"debug: rs_write: start\n") ; }
	while( list != NULL )
	{
		fourcc key = list->key ;
		struct block_functions *block_functions = find_block_functions(key) ;	// gets a set of functions from Global_functions_dictionary for this block type
		if( block_functions == NULL )
		{
			fprintf(stderr,"Cannot write block '%s'\n",strkey(list->key)) ;
			return 1 ;
		}
		int (*gen_function)(struct node *, FILE *) = block_functions->gen ;
		int err = (*gen_function)(list,outfile) ;	// calls the 'gen' function corresponding to the block type
		if( err )
		{
			fprintf(stderr,"Error in '%s' block\n",strkey(key)) ;
			return 1 ;
		}
		list = list->next ;
	}
	if( Debug ) { fprintf(stderr,"debug: rs_write: finish\n") ; }
	return 0 ;
}

void chomp(char *line, int max)		// delete trailing newline character
{
	line[max-1] = '\0' ;
	char *nl = index(line,0x0a) ;
	if( nl ) *nl = '\0' ;
}

int read_parameter(FILE *fd, char format[], void *buffer)	// looks for a line with text to match the given format, copies the value to buffer
{
	size_t block_start = ftell(fd) ;	// put a finger in the file at the first line of the block
	int err = 0 ;
	char line[SIZE_LINE] ;
	if( Debug ) { fprintf(stderr,"debug: read_parameters: format='%s'\n",format) ; }
	while( fgets(line,SIZE_LINE,fd) )
	{
		chomp(line,SIZE_LINE) ;
		if( strlen(line) == 0 )
		{
			fprintf(stderr,"Cannot find parameter '%s'\n",format) ;
			err = 1 ;
			break ;			// stop on blank line
		}
		if( Debug ) { fprintf(stderr,"debug: read_parameters: line='%s'\n",line) ; }
		int count = sscanf(line,format,buffer) ;
		if( Debug ) { fprintf(stderr,"debug: read_parameters: sscanf returned %d\n",count) ; }
		if( count != 0 )
		{
			if( Debug ) { fprintf(stderr,"debug: read_parameters: using line %s\n",line) ; }
			err = 0 ;
			break ;
		}
	}
	fseek(fd,block_start,SEEK_SET) ;
	return err ;
}

int read_binary_file(FILE * rsfile, unsigned long filesize, unsigned char *buffer)		// reads the entire binary RS file into memory
{
	unsigned long count = fread(buffer,1,filesize,rsfile) ;
	if( count != filesize )
	{
		fprintf(stderr,"Error reading rs file, only read %lu bytes out of %lu\n",count,filesize) ;
		return 1 ;
	}
	if( count > sizeof(struct block_header) )
	{
		if( check_header(buffer) )
			return 1 ;
	}
	return 0 ;
}

int check_header(unsigned char *buffer)			// make sure the data in buffer is from an RS file
{
	struct block_header *header = (struct block_header *)buffer ;
	fourcc key = header->key ;
	endian_fixup(&key,sizeof(key)) ;
	if( Debug ) { fprintf(stderr,"debug: check_header: read %x as %x\n",header->key,key) ; }
	if( key != KEY_AQFT )
	{
		fprintf(stderr,"Bad header key: %x\n",key) ;
		return 1 ;
	}
	// TODO check file size (meh?)
	return 0 ;
}

struct node *parse_file(unsigned char *buffer, unsigned long length)	// convert the binary RIFF blocks into a linked list
{
	struct node dummy_root ;	// use a dummy root node to prime the parser
	if( parse_block(&dummy_root,buffer,length) )	// parses the whole file
	{
		fprintf(stderr,"Parser error\n") ;
		return NULL ;
	}
	return dummy_root.next ;	// return the next node as the true root
}

int parse_block(struct node *root, unsigned char *buffer, unsigned long length)		// parser for RIFF blocks, builds the linked list
{
	while( length > 0 )
	{
		// make a new node to describe this block
		struct node *newnode = malloc(sizeof(struct node)) ;
		if( newnode == NULL )
		{
			fprintf(stderr,"Malloc error\n") ;
			return 1 ;
		}
		memset(newnode,0,sizeof(struct node)) ;
		root->next = newnode ;					// link the previous node to the new node
		struct block_header *header = (struct block_header *)buffer ;	// make buffer contents accessible through struct block_header
		endian_fixup(&(header->key),sizeof(newnode->key)) ;	// fixup the endian order
		newnode->key = header->key ;				// copy the block type, aka key, from the header
		endian_fixup(&(header->size),sizeof(newnode->size)) ;	// fixup the endian order
		newnode->size = header->size ;				// copy the size of the data block (excluding header)
		length -= sizeof(struct block_header) ;			// reduce the block length by the size of the header
		buffer += sizeof(struct block_header) ;			// advance the buffer pointer by the size of the header
		if( newnode->size > length )
		{
			fprintf(stderr,"Block '%s' size truncted from %u to %lu bytes\n",strkey(newnode->key),newnode->size,length) ;
			newnode->size = length ;
		}
		newnode->data = buffer ;				// point at the data portion of the block
		if( superblock(newnode->key) )				// if the block is a superblock, recursively parse its data block
		{
			if( parse_block(newnode,newnode->data,newnode->size) )
				return 1 ;
		}
		else
		{
			if( fixup_block(newnode) )			// otherwise, the node's data portion just needs an endian fixup
				;
				//return 1 ;				// continue on error
		}
		// move on to the next block in the buffer
		length -= newnode->size ;				// reduce the block length by the size of the data block
		buffer += newnode->size ;				// advance the buffer pointer by the size of the data block
		do{							// iterate to deal with recursive parsing
			root = root->next ;				// advance the root pointer to the new node, ready to iterate
		}while( root->next != NULL ) ;
	}
	return 0 ;
}

struct block_functions Global_function_dictionary[] =		// a list of RIFF keys and associated functions, used to lookup which function to call
{
	{ KEY_AQFT, fixup_data_aqft, make_node_aqft, dump_block_aqft, gen_block_aqft  },
	{ KEY_HEAD, fixup_data_head, make_node_head, dump_block_head, gen_block_head  },
	{ KEY_sign, fixup_data_sign, make_node_sign, dump_block_sign, gen_block_sign  },
	{ KEY_mcda, fixup_data_mcda, make_node_mcda, dump_block_mcda, gen_block_mcda  },
	{ KEY_dbrf, fixup_data_dbrf, make_node_dbrf, dump_block_dbrf, gen_block_dbrf  },
	{ KEY_cnst, fixup_data_cnst, make_node_cnst, dump_block_cnst, gen_block_cnst  },
	{ KEY_hasi, fixup_data_hasi, make_node_hasi, dump_block_hasi, gen_block_hasi  },
	{ KEY_swep, fixup_data_swep, make_node_swep, dump_block_swep, gen_block_swep  },
	{ KEY_fbin, fixup_data_fbin, make_node_fbin, dump_block_fbin, gen_block_fbin  },
	{ KEY_BODY, fixup_data_body, make_node_body, dump_block_body, gen_block_body  },
	{ KEY_rtag, fixup_data_rtag, make_node_rtag, dump_block_rtag, gen_block_rtag  },
	{ KEY_gps1, fixup_data_gps1, make_node_gps1, dump_block_gps1, gen_block_gps1  },
	{ KEY_indx, fixup_data_indx, make_node_indx, dump_block_indx, gen_block_indx  },
	{ KEY_scal, fixup_data_scal, make_node_scal, dump_block_scal, gen_block_scal  },
	{ KEY_afft, fixup_data_afft, make_node_afft, dump_block_afft, gen_block_afft  },
	{ KEY_ifft, fixup_data_ifft, make_node_ifft, dump_block_ifft, gen_block_ifft  },
	{ KEY_END, fixup_data_end, make_node_end, dump_block_end, gen_block_end  },
	{ 0, NULL, NULL, NULL, NULL }
} ;

struct block_functions *find_block_functions(fourcc key)	// search for the given RIFF key in the list of keys/functions, return a struct of functions
{
	if( key == 0 )
	{
		fprintf(stderr,"Bad key (zero!)\n") ;
		return NULL ;
	}
	for( struct block_functions *block_functions = Global_function_dictionary ; block_functions->key != 0 ; block_functions++ )
	{
		if( key == block_functions->key )
			return block_functions ;
	}
	fprintf(stderr,"Cannot locate functions for key '%s'\n",strkey(key)) ;
	return NULL ;
}

int fixup_block(struct node *node)	// calls the fixup function that is associated with the given node's RIFF key
{
	struct block_functions *block_functions = find_block_functions(node->key) ;	// returns a set of functions from the Global_functions_dictionary for this block type
	if( block_functions == NULL )
	{
		fprintf(stderr,"Cannot fixup block '%s'\n",strkey(node->key)) ;
		return 1 ;
	}
	int (*fixup_function)(struct node *) = block_functions->fixup ;
	int err = (*fixup_function)(node) ;	// calls the fixup function corresponding to the block key
	if( err )
	{
		fprintf(stderr,"Error fixing block %s\n",strkey(node->key)) ;
		return 1 ;
	}
	return 0 ;
}

void endian_fixup(void *original, int size)	// performs byte swapping if needed, as determined by the global endian flag
{
	if( Global_flag_little_endian )
	{
		unsigned char tmp[8] ;
		memcpy(tmp,original,size) ;	// uses a temporary copy
		swapcopy(original,tmp,size) ;
	}
}

void swapcopy(unsigned char *dest, unsigned char *source, int size)	// swaps bytes according to the size of the object
{
	switch( size )
	{
		case 2:
			swapcopy2(dest,source) ;
		break ;
		case 4:
			swapcopy4(dest,source) ;
		break ;
		case 8:
			swapcopy8(dest,source) ;
		break ;
	}
}

void swapcopy2(unsigned char *dest, unsigned char *source)
{
	dest[0] = source[1] ;
	dest[1] = source[0] ;
}

void swapcopy4(unsigned char *dest, unsigned char *source)
{
	dest[0] = source[3] ;
	dest[1] = source[2] ;
	dest[2] = source[1] ;
	dest[3] = source[0] ;
}

void swapcopy8(unsigned char *dest, unsigned char *source)
{
	dest[0] = source[7] ;
	dest[1] = source[6] ;
	dest[2] = source[5] ;
	dest[3] = source[4] ;
	dest[4] = source[3] ;
	dest[5] = source[2] ;
	dest[6] = source[1] ;
	dest[7] = source[0] ;
}

int superblock(fourcc key)	// returns 1 if the RIFF key denotes a 'superblock', i.e. one that is composed of sub-blocks
{
	// check for one of the known superkeys
	if( key == KEY_AQFT ) return 1 ;
	if( key == KEY_HEAD ) return 1 ;
	if( key == KEY_BODY ) return 1 ;
	if( key == KEY_END ) return 1 ;
	return 0 ;		// not a superblock key
}

void debugdump(unsigned char *buffer, int n)	// prettyprint hexascii for debug purposes
{
	for( int row = 0 ; row <= n/8 ; row++ )
	{
		for( int col = 0 ; col < 8 ; col++ )
		{
			int offset = row*8+col ;
			if( offset < n )
			{
				printf("%02x ",buffer[offset]) ;
			}
			else
			{
				printf("   ") ;
			}
		}
		printf("\t") ;
		for( int col = 0 ; col < 8 ; col++ )
		{
			int offset = row*8+col ;
			if( offset < n )
			{
				if( isprint(buffer[offset]) )
					printf("%c",buffer[offset]) ;
				else
					printf(".") ;
			}
			else
			{
				printf("   ") ;
			}
		}
		printf("\n") ;
	}
}

#define MAX_NAME 10
char Global_name[MAX_NAME] ;
char *strkey(fourcc name)		// converts a RIFF key to a string, with byte-swapping. Relies on a global variable.
{
	int size = sizeof(name) ;
	if( size >= MAX_NAME )
		size = MAX_NAME-1 ;
	memcpy(Global_name,(void *)&name,size) ;
	endian_fixup(Global_name,size) ;
	Global_name[size] = '\0' ;
	return Global_name ;
}

int fixup_sizes(struct node *list)	// updates the various block size counts from the data in the linked list
{
	uint32_t head_size = calculate_head_size(list) ;
	if( head_size == 0 ) return 1 ;
	uint32_t body_size = calculate_body_size(list) ;
	if( body_size == 0 ) return 1 ;
	uint32_t aqft_size = head_size + sizeof(struct block_header) + body_size + sizeof(struct block_header) ;
	if( set_block_size(list,KEY_AQFT,aqft_size) ) return 1 ;
	if( set_block_size(list,KEY_HEAD,head_size) ) return 1 ;
	if( set_block_size(list,KEY_BODY,body_size) ) return 1 ;
	return 0 ;
}

uint32_t calculate_body_size(struct node *list)			// returns the size of the blocks in the BODY superblock
{
	int in_body = 0 ;
	uint32_t size = 0 ;
	while( list != NULL )
	{
		if( list->key == KEY_END )
			in_body = 0 ;
		if( in_body )
			size += ( list->size + sizeof(struct block_header) ) ;	// remember to count the block header
		if( list->key == KEY_BODY )
			in_body = 1 ;
		list = list->next ;
	} ;
	return size ;
}

uint32_t calculate_head_size(struct node *list)			// returns the size of the blocks in the HEAD superblock
{
	int in_head = 0 ;
	uint32_t size = 0 ;
	while( list != NULL )
	{
		if( list->key == KEY_END )
			in_head = 0 ;
		if( list->key == KEY_BODY )
			in_head = 0 ;
		if( in_head )
			size += (list->size + sizeof(struct block_header) ) ;	// remember to count the block header
		if( list->key == KEY_HEAD )
			in_head = 1 ;
		list = list->next ;
	} ;
	return size ;
}

int set_block_size(struct node *list, fourcc key, uint32_t size)	// searches the linked list for a given RIFF key to update the size data
{
	while( list != NULL )
	{
		if( list->key == key )
		{
			list->size = size ;
			return 0 ;
		}
		list = list->next ;
	}
	return 1 ;	// returns 1 for failure, 0 for success
}

void show_list(struct node *list)	// show the linked list, for debugging
{
	unsigned int count = 0 ;
	while( list != NULL )
	{
		printf("Node %u: key %.4s size %u\n",count,(char *)&(list->key),list->size) ;
		list = list->next ;
		count++ ;
	}
}

int dump_list(struct node *list, FILE *outfile, int just_header) // goes through the list of nodes, writing an ascii text description of each node to outfile
{
	unsigned int count = 0 ;
	struct config config ;
	memset(&config,0,sizeof(struct config)) ;
	while( list != NULL )
	{
		if( Debug ) { fprintf(stderr,"debug: dump_list: node has key '%s'\n",strkey(list->key)) ; }
		struct block_functions *block_functions = find_block_functions(list->key) ;	// returns a set of functions from the Global_functions_dictionary for this block key
		if( block_functions == NULL )
		{
			fprintf(stderr,"Cannot dump block '%s'\n",strkey(list->key)) ;
			return 1 ;
		}
		int (*dump_function)(struct node *, struct config *, FILE *) = block_functions->dump ;	// extract the dump function
		if( just_header && dump_function == dump_block_body ) return 0 ;
		int err = (*dump_function)(list,&config,outfile) ;	// calls the dump function
		if( err )
		{
			fprintf(stderr,"Error dumping block '%s'\n",strkey(list->key)) ;
			return 1 ;
		}
		list = list->next ;
	}
	return 0 ;
}


// Start of the block-specific functions.
// For each block type there are four functions: fixup_data_xxxx, dump_block_xxxx, make_node_xxxx and gen_block_xxxx.
// The function fixup_data_xxxx performs endian fixup on data read from a binary RIFF. Used in rsdump mode.
// The function dump_block_xxxx writes a text version of a block from a linked list node. Used in rsdump mode.
// The function make_node_xxxx reads a text version of the block and makes a linked list node. Used in rsgen mode.
// The function gen_block_xxxx writes a binary RIFF block from a linked list node. Used in rsgen mode.


int fixup_data_aqft(struct node *node)
{
	return 0 ;
}

int dump_block_aqft(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_AQFT)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_aqft(struct node *list, struct config *config, FILE *fd)		// creates a new node for AQFT block, fd not used
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_AQFT ;
	// fixup size at the very end
	// no explicit data block, it's composed of sub blocks
	return 0 ;
}

int gen_block_aqft(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


int fixup_data_head(struct node *node)
{
	return 0 ;
}

int dump_block_head(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_HEAD)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_head(struct node *list, struct config *config, FILE *fd)
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_HEAD ;
	// fixup size at the very end
	// no explicit data block, it's composed of sub blocks
	return 0 ;
}

int gen_block_head(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_sign				// contains only the data portion for this type of block, i.e. no header data (key, size)
{
	fourcc version ;			// file version
	fourcc filetype ;			// file type
	fourcc sitecode ; 			// owner code
	uint32_t userflags ;			// user flags
	char description[SIZE_DESCRIPTION] ;	// file description
	char ownername[SIZE_OWNERNAME] ;	// owner name
	char comment[SIZE_COMMENT] ;		// comment
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_sign(struct node *node)
{
	if( node->size < sizeof(struct block_sign) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_sign)) ;
		return 1 ;
	}
	struct block_sign *sign = (struct block_sign *)(node->data) ;
	endian_fixup(&(sign->version),sizeof(sign->version)) ;
	endian_fixup(&(sign->filetype),sizeof(sign->filetype)) ;
	endian_fixup(&(sign->sitecode),sizeof(sign->sitecode)) ;
	endian_fixup(&(sign->userflags),sizeof(sign->userflags)) ;
	return 0 ;
}

int dump_block_sign(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_sign) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_sign)) ;
		return 1 ;
	}
	struct block_sign *sign = (struct block_sign *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_sign)) ;
	fprintf(outfile,"version:%s\n",strkey(sign->version)) ;
	fprintf(outfile,"filetype:%s\n",strkey(sign->filetype)) ;
	fprintf(outfile,"sitecode:%s\n",strkey(sign->sitecode)) ;
	fprintf(outfile,"userflags:%x\n",sign->userflags) ;
	fprintf(outfile,"description:%.*s\n",SIZE_DESCRIPTION,sign->description) ;
	fprintf(outfile,"ownername:%.*s\n",SIZE_OWNERNAME,sign->ownername) ;
	fprintf(outfile,"comment:%.*s\n",SIZE_OWNERNAME,sign->comment) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_sign(struct node *list, struct config *config, FILE *fd)
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_sign ;
	struct block_sign *sign = malloc(sizeof(struct block_sign)) ;
	if( sign == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(sign,0,sizeof(struct block_sign)) ;
	newnode->data = (unsigned char *)sign ;
	newnode->size = sizeof(struct block_sign) ;
	if( read_parameter(fd,"version:%4c",(void *)&(sign->version)) ) return 1 ;
	if( read_parameter(fd,"filetype:%4c",(void *)&(sign->filetype)) ) return 1 ;
	if( read_parameter(fd,"sitecode:%4c",(void *)&(sign->sitecode)) ) return 1 ;
	if( read_parameter(fd,"userflags:%x",(void *)&(sign->userflags)) ) return 1 ;
	char format[32] ;
	sprintf(format,"description:%%%dc",SIZE_DESCRIPTION) ;
	if( read_parameter(fd,format,(void *)(sign->description)) ) return 1 ;
	sprintf(format,"ownername:%%%dc",SIZE_OWNERNAME) ;
	if( read_parameter(fd,format,(void *)(sign->ownername)) ) return 1 ;
	sprintf(format,"comment:%%%dc",SIZE_COMMENT) ;
	if( read_parameter(fd,format,(void *)(sign->comment)) ) return 1 ;
	return 0 ;
}

int gen_block_sign(struct node *node, FILE *outfile)
{
	struct block_sign *sign = (struct block_sign *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->version),sizeof(sign->version)) ;
	if( fwrite(&(sign->version),sizeof(sign->version),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->filetype),sizeof(sign->filetype)) ;
	if( fwrite(&(sign->filetype),sizeof(sign->filetype),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->sitecode),sizeof(sign->sitecode)) ;
	if( fwrite(&(sign->sitecode),sizeof(sign->sitecode),1,outfile) != 1 ) return 1 ;
	//endian_fixup(&(sign->userflags),sizeof(sign->userflags)) ;
	if( fwrite(&(sign->userflags),sizeof(sign->userflags),1,outfile) != 1 ) return 1 ;
	if( fwrite(&(sign->description),SIZE_DESCRIPTION,1,outfile) != 1 ) return 1 ;
	if( fwrite(&(sign->ownername),SIZE_OWNERNAME,1,outfile) != 1 ) return 1 ;
	if( fwrite(&(sign->comment),SIZE_COMMENT,1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_mcda
{
	uint32_t filetimestamp ;		// mac timestamp of first sweep
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_mcda(struct node *node)
{
	if( node->size < sizeof(struct block_mcda) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_mcda)) ;
		return 1 ;
	}
	struct block_mcda *mcda = (struct block_mcda *)(node->data) ;
	endian_fixup(&(mcda->filetimestamp),sizeof(mcda->filetimestamp)) ;
	return 0 ;
}

int dump_block_mcda(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_mcda) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_mcda)) ;
		return 1 ;
	}
	struct block_mcda *mcda = (struct block_mcda *)(node->data) ;
	uint32_t mcda_time = mcda->filetimestamp ;
	time_t t = mcda_time ;
	fprintf(outfile,"%s\n",strkey(KEY_mcda)) ;
	if( t != 0 )
	{
		t -= 2082844800 ;	// move epoc from 1904-01-01 00:00:00 to 1970-01-01 00:00:00 
		fprintf(outfile,"filetimestamp:%lu (NB: seconds since 1970) (%.24s)\n",t,ctime(&t)) ;
	}
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_mcda(struct node *list, struct config *config, FILE *fd)		// creates a new node for mcda block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_mcda ;
	struct block_mcda *mcda = malloc(sizeof(struct block_mcda)) ;
	if( mcda == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(mcda,0,sizeof(struct block_mcda)) ;
	newnode->data = (unsigned char *)mcda ;
	newnode->size = sizeof(struct block_mcda) ;
	if( read_parameter(fd,"filetimestamp:%u",(void *)&(mcda->filetimestamp)) ) return 1 ;
	mcda->filetimestamp += 2082844800 ;	// move epoc from 1970-01-01 00:00:00 to 1904-01-01 00:00:00 
	return 0 ;
}

int gen_block_mcda(struct node *node, FILE *outfile)
{
	struct block_mcda *mcda = (struct block_mcda *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(mcda->filetimestamp),sizeof(mcda->filetimestamp)) ;
	if( fwrite(&(mcda->filetimestamp),sizeof(mcda->filetimestamp),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_dbrf
{
	double rxloss ;		// received power correction, in dB
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_dbrf(struct node *node)
{
	if( node->size < sizeof(struct block_dbrf) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_dbrf)) ;
		return 1 ;
	}
	struct block_dbrf *dbrf = (struct block_dbrf *)(node->data) ;
	endian_fixup(&(dbrf->rxloss),sizeof(dbrf->rxloss)) ;
	return 0 ;
}

int dump_block_dbrf(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_dbrf) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_dbrf)) ;
		return 1 ;
	}
	struct block_dbrf *dbrf = (struct block_dbrf *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_dbrf)) ;
	fprintf(outfile,"rxloss:%.4lf\n",dbrf->rxloss) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_dbrf(struct node *list, struct config *config, FILE *fd)		// creates a new node for dbrf block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_dbrf ;
	struct block_dbrf *dbrf = malloc(sizeof(struct block_dbrf)) ;
	if( dbrf == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(dbrf,0,sizeof(struct block_dbrf)) ;
	newnode->data = (unsigned char *)dbrf ;
	newnode->size = sizeof(struct block_dbrf) ;
	if( read_parameter(fd,"rxloss:%lf",(void *)&(dbrf->rxloss)) ) return 1 ;
	return 0 ;
}

int gen_block_dbrf(struct node *node, FILE *outfile)
{
	struct block_dbrf *dbrf = (struct block_dbrf *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(dbrf->rxloss),sizeof(dbrf->rxloss)) ;
	if( fwrite(&(dbrf->rxloss),sizeof(dbrf->rxloss),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_cnst
{
	int32_t nchannels ;		// number of antennas/channels (normally 3)
	int32_t nranges ;		// number of ranges
	int32_t nsweeps ;		// number of sweeps (normally 32)
	int32_t iqindicator ;		// iqindicator: 1=?, 2=IQ
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_cnst(struct node *node)
{
	if( node->size < sizeof(struct block_cnst) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_cnst)) ;
		return 1 ;
	}
	struct block_cnst *cnst = (struct block_cnst *)node->data ;
	endian_fixup(&(cnst->nchannels),sizeof(cnst->nchannels)) ;
	endian_fixup(&(cnst->nranges),sizeof(cnst->nranges)) ;
	endian_fixup(&(cnst->nsweeps),sizeof(cnst->nchannels)) ;
	endian_fixup(&(cnst->iqindicator),sizeof(cnst->iqindicator)) ;
	return 0 ;
}

int dump_block_cnst(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_cnst) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_cnst)) ;
		return 1 ;
	}
	struct block_cnst *cnst = (struct block_cnst *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_cnst)) ;
	fprintf(outfile,"nchannels:%d\n",cnst->nchannels) ;
	fprintf(outfile,"nranges:%d\n",cnst->nranges) ;
	fprintf(outfile,"nsweeps:%d\n",cnst->nsweeps) ;
	fprintf(outfile,"iqindicator:%d\n",cnst->iqindicator) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_cnst(struct node *list, struct config *config, FILE *fd)		// creates a new node for cnst block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_cnst ;
	struct block_cnst *cnst = malloc(sizeof(struct block_cnst)) ;
	if( cnst == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(cnst,0,sizeof(struct block_cnst)) ;
	newnode->data = (unsigned char *)cnst ;
	newnode->size = sizeof(struct block_cnst) ;
	if( read_parameter(fd,"nchannels:%d",(void *)&(cnst->nchannels)) ) return 1 ;
	if( read_parameter(fd,"nranges:%d",(void *)&(cnst->nranges)) ) return 1 ;
	if( read_parameter(fd,"nsweeps:%d",(void *)&(cnst->nsweeps)) ) return 1 ;
	if( read_parameter(fd,"iqindicator:%d",(void *)&(cnst->iqindicator)) ) return 1 ;
	return 0 ;
}

int gen_block_cnst(struct node *node, FILE *outfile)
{
	struct block_cnst *cnst = (struct block_cnst *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->nchannels),sizeof(cnst->nchannels)) ;
	if( fwrite(&(cnst->nchannels),sizeof(cnst->nchannels),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->nranges),sizeof(cnst->nranges)) ;
	if( fwrite(&(cnst->nranges),sizeof(cnst->nranges),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->nsweeps),sizeof(cnst->nsweeps)) ;
	if( fwrite(&(cnst->nsweeps),sizeof(cnst->nsweeps),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(cnst->iqindicator),sizeof(cnst->iqindicator)) ;
	if( fwrite(&(cnst->iqindicator),sizeof(cnst->iqindicator),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_hasi
{
	uint32_t hasi ;		// hasi (?)
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_hasi(struct node *node)
{
	if( node->size < sizeof(struct block_hasi) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_hasi)) ;
		return 1 ;
	}
	//struct block_hasi *hasi = (struct block_hasi *)node->data ;
	//endian_fixup(&(hasi->hasi),sizeof(hasi->hasi)) ;
	return 0 ;
}

int dump_block_hasi(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_hasi) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_hasi)) ;
		return 1 ;
	}
	fprintf(outfile,"%s\n",strkey(KEY_hasi)) ;
	// unknown structure, dump the raw data bytes
	hexdump(node->data,node->size,outfile) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

void hexdump(unsigned char *data, unsigned int size, FILE *outfile)
{
	int loop ;
	fprintf(outfile,"data:") ;
	for( loop = 0 ; loop < size ; loop++, data++ )
	{
		fprintf(outfile," %02x",data[0]) ;
	}
	fprintf(outfile,"\n") ;
}

#define MAX_ARGC 256
#define MAX_LINE MAX_ARGC*3+1

int make_node_hasi(struct node *list, struct config *config, FILE *fd)		// creates a new node for hasi block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_hasi ;
	// read a line of unformatted data
	char line[MAX_LINE] ;
	if( fgets(line,MAX_LINE-1,fd) == 0 )
	{
		if( Debug ) { fprintf(stderr,"debug: make_node_hasi: fgets returned 0\n") ; }
		return -1 ;
	}
	line[MAX_LINE-1] = '\0' ;
	chomp(line,MAX_LINE) ;
	// split the line into space-separated args, storing the argv
	char sep[] = " " ;
	char *argv[MAX_ARGC] ;
	int argc = 0 ;
	argv[0] = strtok(line,sep) ;	// the first arg is "data:"
	argv[0] = strtok(NULL,sep) ;	// get the first data byte
	while( argv[argc] && argc < MAX_ARGC )
	{
		argv[++argc] = strtok(NULL,sep) ;
	}
	if( argc <= 0 )
	{
		if( Debug ) { fprintf(stderr,"debug: make_node_hasi: argc <= 0\n") ; }
		return 1 ;
	}
	if( argc == MAX_ARGC ) return 1 ;
	if( Debug ) { fprintf(stderr,"debug: make_node_hasi: counted %d bytes of data\n",argc) ; }
	// use the arg count to malloc space for the data
	unsigned char *hasi = malloc(argc) ;
	if( hasi == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(hasi,0,argc) ;
	newnode->data = hasi ;
	newnode->size = argc ;
	// TODO: read a line of ascii hex
	for( int i = 0 ; i < argc ; i++ )
	{
		unsigned char c = (unsigned char )strtoul(argv[i],NULL,16) ;
		hasi[i] = c ;
		if( Debug ) { fprintf(stderr,"debug: make_node_hasi: data[%d]=%02x\n",i,c) ; }
	}
	return 0 ;
}

int gen_block_hasi(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	unsigned int actualsize = node->size ;
	if( Debug ) { fprintf(stderr,"debug: gen_block_hasi: actualsize=%u\n",actualsize) ; }
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	if( fwrite(node->data,1,actualsize,outfile) != actualsize )	// no data structure, just write size bytes
	{
		if( Debug ) { fprintf(stderr,"debug: gen_block_hasi: error on fwrite for node->data\n") ; }
		return 1 ;
	}
	return 0 ;
}


struct block_swep
{
	int32_t samplespersweep ;	// number of samples per sweep/channels (normally 2048)
	double sweepstart ;		// sweep start frequency in Hertz
	double sweepbandwidth ;		// sweep bandwidth in Hertz
	double sweeprate ;		// sweep rate in Hertz
	int32_t rangeoffset ;		// rangeoffset (not used)
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_swep(struct node *node)
{
	if( node->size < sizeof(struct block_swep) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_swep)) ;
		return 1 ;
	}
	struct block_swep *swep = (struct block_swep *)(node->data) ;
	endian_fixup(&(swep->samplespersweep),sizeof(swep->samplespersweep)) ;
	endian_fixup(&(swep->sweepstart),sizeof(swep->sweepstart)) ;
	endian_fixup(&(swep->sweepbandwidth),sizeof(swep->sweepbandwidth)) ;
	endian_fixup(&(swep->sweeprate),sizeof(swep->sweeprate)) ;
	endian_fixup(&(swep->rangeoffset),sizeof(swep->rangeoffset)) ;
	return 0 ;
}

int dump_block_swep(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_swep) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_swep)) ;
		return 1 ;
	}
	struct block_swep *swep = (struct block_swep *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_swep)) ;
	fprintf(outfile,"samplespersweep:%d\n",swep->samplespersweep) ;
	fprintf(outfile,"sweepstart:%.20lf\n",swep->sweepstart) ;
	fprintf(outfile,"sweepbandwidth:%.20lf\n",swep->sweepbandwidth) ;
	fprintf(outfile,"sweeprate:%.20lf\n",swep->sweeprate) ;
	fprintf(outfile,"rangeoffset:%d\n",swep->rangeoffset) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_swep(struct node *list, struct config *config, FILE *fd)		// creates a new node for swep block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_swep ;
	struct block_swep *swep = malloc(sizeof(struct block_swep)) ;
	if( swep == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(swep,0,sizeof(struct block_swep)) ;
	newnode->data = (unsigned char *)swep ;
	newnode->size = sizeof(struct block_swep) ;
	if( read_parameter(fd,"samplespersweep:%d",(void *)&(swep->samplespersweep)) ) return 1 ;
	if( read_parameter(fd,"sweepstart:%lf",(void *)&(swep->sweepstart)) ) return 1 ;
	if( read_parameter(fd,"sweepbandwidth:%lf",(void *)&(swep->sweepbandwidth)) ) return 1 ;
	if( read_parameter(fd,"sweeprate:%lf",(void *)&(swep->sweeprate)) ) return 1 ;
	if( read_parameter(fd,"rangeoffset:%d",(void *)&(swep->rangeoffset)) ) return 1 ;
	return 0 ;
}

int gen_block_swep(struct node *node, FILE *outfile)
{
	struct block_swep *swep = (struct block_swep *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->samplespersweep),sizeof(swep->samplespersweep)) ;
	if( fwrite(&(swep->samplespersweep),sizeof(swep->samplespersweep),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->sweepstart),sizeof(swep->sweepstart)) ;
	if( fwrite(&(swep->sweepstart),sizeof(swep->sweepstart),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->sweepbandwidth),sizeof(swep->sweepbandwidth)) ;
	if( fwrite(&(swep->sweepbandwidth),sizeof(swep->sweepbandwidth),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->sweeprate),sizeof(swep->sweeprate)) ;
	if( fwrite(&(swep->sweeprate),sizeof(swep->sweeprate),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(swep->rangeoffset),sizeof(swep->rangeoffset)) ;
	if( fwrite(&(swep->rangeoffset),sizeof(swep->rangeoffset),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_fbin
{
	fourcc bin_format ;	// format (normally 'cviq')
	fourcc bin_type ;	// type of ALVL data ('flt8','flt4','fix4','fix3','fix2')
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_fbin(struct node *node)
{
	if( node->size < sizeof(struct block_fbin) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_fbin)) ;
		return 1 ;
	}
	struct block_fbin *fbin = (struct block_fbin *)node->data ;
	endian_fixup(&(fbin->bin_format),sizeof(fbin->bin_format)) ;
	endian_fixup(&(fbin->bin_type),sizeof(fbin->bin_type)) ;
	return 0 ;
}

int dump_block_fbin(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_fbin) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_fbin)) ;
		return 1 ;
	}
	struct block_fbin *fbin = (struct block_fbin *)(node->data) ;
	config->bin_format = fbin->bin_format ;	// remember this for body data blocks
	config->bin_type = fbin->bin_type ;	// remember this for body data blocks
	fprintf(outfile,"%s\n",strkey(KEY_fbin)) ;
	fprintf(outfile,"format:%s\n",strkey(fbin->bin_format)) ;
	fprintf(outfile,"type:%s\n",strkey(fbin->bin_type)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_fbin(struct node *list, struct config *config, FILE *fd)		// creates a new node for fbin block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_fbin ;
	struct block_fbin *fbin = malloc(sizeof(struct block_fbin)) ;
	if( fbin == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(fbin,0,sizeof(struct block_fbin)) ;
	newnode->data = (unsigned char *)fbin ;
	newnode->size = sizeof(struct block_fbin) ;
	char format[16] ;
	sprintf(format,"format:%%%lus",sizeof(fbin->bin_format)) ;
	if( read_parameter(fd,format,(void *)&(fbin->bin_format)) ) return 1 ;	// read as a 4 byte string
	endian_fixup(&(fbin->bin_format),sizeof(fbin->bin_format)) ;		// then endian correct to 4 bytes int
	sprintf(format,"type:%%%lus",sizeof(fbin->bin_type)) ;
	if( read_parameter(fd,format,(void *)&(fbin->bin_type)) ) return 1 ;	// read as a 4 byte string
	endian_fixup(&(fbin->bin_type),sizeof(fbin->bin_type)) ;		// then endian correct to 4 bytes int
	config->bin_format = fbin->bin_format ;	// remember this for afft blocks
	config->bin_type = fbin->bin_type ;	// remember this for afft blocks
	if( Debug ) { fprintf(stderr,"debug: make_node_fbin: fbin->bin_format=%s fbin->bin_type=%s\n",strkey(fbin->bin_format),strkey(fbin->bin_type)) ; }
	return 0 ;
}

int gen_block_fbin(struct node *node, FILE *outfile)
{
	struct block_fbin *fbin = (struct block_fbin *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(fbin->bin_format),sizeof(fbin->bin_format)) ;
	if( fwrite(&(fbin->bin_format),sizeof(fbin->bin_format),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(fbin->bin_type),sizeof(fbin->bin_type)) ;
	if( fwrite(&(fbin->bin_type),sizeof(fbin->bin_type),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


int fixup_data_body(struct node *node)
{
	return 0 ;
}

int dump_block_body(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_BODY)) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_body(struct node *list, struct config *config, FILE *fd)		// creates a new node for body block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_BODY ;
	// fixup size at the very end
	// no explicit data block, it's composed of sub blocks
	return 0 ;
}

int gen_block_body(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_rtag
{
	uint32_t rtag ;		// receiver position (?)
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_rtag(struct node *node)
{
	if( node->size < sizeof(struct block_rtag) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_rtag)) ;
		return 1 ;
	}
	struct block_rtag *rtag = (struct block_rtag *)node->data ;
	endian_fixup(&(rtag->rtag),sizeof(rtag->rtag)) ;			// fixup the endian order
	return 0 ;
}

int dump_block_rtag(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_rtag) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_rtag)) ;
		return 1 ;
	}
	struct block_rtag *rtag = (struct block_rtag *)(node->data) ;
	fprintf(outfile,"%s\n",strkey(KEY_rtag)) ;
	fprintf(outfile,"rtag:%u\n",rtag->rtag) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_rtag(struct node *list, struct config *config, FILE *fd)		// creates a new node for rtag block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_rtag ;
	struct block_rtag *rtag = malloc(sizeof(struct block_rtag)) ;
	if( rtag == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(rtag,0,sizeof(struct block_rtag)) ;
	newnode->data = (unsigned char *)rtag ;
	newnode->size = sizeof(struct block_rtag) ;
	if( read_parameter(fd,"rtag:%u",(void *)&(rtag->rtag)) ) return 1 ;
	return 0 ;
}

int gen_block_rtag(struct node *node, FILE *outfile)
{
	struct block_rtag *rtag = (struct block_rtag *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(rtag->rtag),sizeof(rtag->rtag)) ;
	if( fwrite(&(rtag->rtag),sizeof(rtag->rtag),1,outfile) != 1 ) return 1 ;
	return 0 ;
}

struct block_gps1
{
	double lat ;				// GPS latitude in radians
	double lon ;				// GPS longitude in radians
	double alt ;				// GPS altitude in meters
	int32_t gpstimestamp ;			// GPS timestamp in appletime
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_gps1(struct node *node)
{
	if( node->size < sizeof(struct block_gps1) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_gps1)) ;
		return 1 ;
	}
	struct block_gps1 *gps1 = (struct block_gps1 *)node->data ;
	endian_fixup(&(gps1->lat),sizeof(gps1->lat)) ;				// fixup the endian order
	endian_fixup(&(gps1->lon),sizeof(gps1->lon)) ;				// fixup the endian order
	endian_fixup(&(gps1->alt),sizeof(gps1->alt)) ;				// fixup the endian order
	endian_fixup(&(gps1->gpstimestamp),sizeof(gps1->gpstimestamp)) ;	// fixup the endian order
	return 0 ;
}

int dump_block_gps1(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_gps1) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_gps1)) ;
		return 1 ;
	}
	struct block_gps1 *gps1 = (struct block_gps1 *)(node->data) ;
	uint32_t gpstimestamp = gps1->gpstimestamp ;
	time_t t = gpstimestamp ;
	fprintf(outfile,"%s\n",strkey(KEY_gps1)) ;
	fprintf(outfile,"lat:%.6lf\n",gps1->lat) ;
	fprintf(outfile,"lon:%.6lf\n",gps1->lon) ;
	fprintf(outfile,"alt:%.6lf\n",gps1->alt) ;
	if( t != 0 )
	{
		t -= 2082844800 ;	// move epoc from 1904-01-01 00:00:00 to 1970-01-01 00:00:00 
		fprintf(outfile,"gpstimestamp:%lu (NB: seconds since 1970) (%.24s)\n",t,ctime(&t)) ;
	}
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_gps1(struct node *list, struct config *config, FILE *fd)		// creates a new node for gps1 block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_gps1 ;
	struct block_gps1 *gps1 = malloc(sizeof(struct block_gps1)) ;
	if( gps1 == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(gps1,0,sizeof(struct block_gps1)) ;
	newnode->data = (unsigned char *)gps1 ;
	newnode->size = sizeof(struct block_gps1) ;
	if( read_parameter(fd,"lat:%lf",(void *)&(gps1->lat)) ) return 1 ;
	if( read_parameter(fd,"lon:%lf",(void *)&(gps1->lon)) ) return 1 ;
	if( read_parameter(fd,"alt:%lf",(void *)&(gps1->alt)) ) return 1 ;
	if( read_parameter(fd,"gpstimestamp:%u",(void *)&(gps1->gpstimestamp)) ) return 1 ;
	gps1->gpstimestamp += 2082844800 ;	// move epoc from 1970-01-01 00:00:00 to 1904-01-01 00:00:00 
	return 0 ;
}

int gen_block_gps1(struct node *node, FILE *outfile)
{
	struct block_gps1 *gps1 = (struct block_gps1 *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(gps1->lat),sizeof(gps1->lat)) ;
	if( fwrite(&(gps1->lat),sizeof(gps1->lat),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(gps1->lon),sizeof(gps1->lon)) ;
	if( fwrite(&(gps1->lon),sizeof(gps1->lon),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(gps1->alt),sizeof(gps1->alt)) ;
	if( fwrite(&(gps1->alt),sizeof(gps1->alt),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(gps1->gpstimestamp),sizeof(gps1->gpstimestamp)) ;
	if( fwrite(&(gps1->gpstimestamp),sizeof(gps1->gpstimestamp),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_indx
{
	uint32_t index ;		// index
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_indx(struct node *node)
{
	if( node->size < sizeof(struct block_indx) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_indx)) ;
		return 1 ;
	}
	struct block_indx *indx = (struct block_indx *)node->data ;
	endian_fixup(&(indx->index),sizeof(indx->index)) ;			// fixup the endian order
	return 0 ;
}

int dump_block_indx(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_indx) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_indx)) ;
		return 1 ;
	}
	struct block_indx *indx = (struct block_indx *)(node->data) ;
	config->index = indx->index ;
	fprintf(outfile,"%s\n",strkey(KEY_indx)) ;
	fprintf(outfile,"index:%u\n",indx->index) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_indx(struct node *list, struct config *config, FILE *fd)		// creates a new node for indx block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_indx ;
	struct block_indx *indx = malloc(sizeof(struct block_indx)) ;
	if( indx == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(indx,0,sizeof(struct block_indx)) ;
	newnode->data = (unsigned char *)indx ;
	newnode->size = sizeof(struct block_indx) ;
	if( read_parameter(fd,"index:%u",(void *)&(indx->index)) ) return 1 ;
	return 0 ;
}

int gen_block_indx(struct node *node, FILE *outfile)
{
	struct block_indx *indx = (struct block_indx *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(indx->index),sizeof(indx->index)) ;
	if( fwrite(&(indx->index),sizeof(indx->index),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_scal
{
	double scalar_one ;		// scaling value for I samples
	double scalar_two ;		// scaling value for Q samples
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_scal(struct node *node)
{
	if( node->size < sizeof(struct block_scal) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_scal)) ;
		return 1 ;
	}
	struct block_scal *scal = (struct block_scal *)node->data ;
	endian_fixup(&(scal->scalar_one),sizeof(scal->scalar_one)) ;
	endian_fixup(&(scal->scalar_two),sizeof(scal->scalar_two)) ;
	return 0 ;
}

int dump_block_scal(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_scal) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_scal)) ;
		return 1 ;
	}
	struct block_scal *scal = (struct block_scal *)(node->data) ;
	config->scalar_one = scal->scalar_one ;	// remember this for iqdata blocks
	config->scalar_two = scal->scalar_two ;	// remember this for iqdata blocks
	fprintf(outfile,"%s\n",strkey(KEY_scal)) ;
	fprintf(outfile,"scalar_one:%.20lf\n",scal->scalar_one) ;
	fprintf(outfile,"scalar_two:%.20lf\n",scal->scalar_two) ;
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_scal(struct node *list, struct config *config, FILE *fd)		// creates a new node for scal block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_scal ;
	struct block_scal *scal = malloc(sizeof(struct block_scal)) ;
	if( scal == NULL )
	{
		fprintf(stderr,"Malloc error on data block\n") ;
		return 1 ;
	}
	memset(scal,0,sizeof(struct block_scal)) ;
	newnode->data = (unsigned char *)scal ;
	newnode->size = sizeof(struct block_scal) ;
	if( read_parameter(fd,"scalar_one:%lf",(void *)&(scal->scalar_one)) ) return 1 ;
	if( read_parameter(fd,"scalar_two:%lf",(void *)&(scal->scalar_two)) ) return 1 ;
	config->scalar_one = scal->scalar_one ;	// remember this for iqdata blocks
	config->scalar_two = scal->scalar_two ;	// remember this for iqdata blocks
	return 0 ;
}

int gen_block_scal(struct node *node, FILE *outfile)
{
	struct block_scal *scal = (struct block_scal *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(scal->scalar_one),sizeof(scal->scalar_one)) ;
	if( fwrite(&(scal->scalar_one),sizeof(scal->scalar_one),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(scal->scalar_two),sizeof(scal->scalar_two)) ;
	if( fwrite(&(scal->scalar_two),sizeof(scal->scalar_two),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


struct block_iqdata_float		// hardcoded type
{
	float isample ;		// I sample		// hardcoded type
	float qsample ;		// Q sample		// hardcoded type
} __attribute__((packed)) ;	// make sure there's no padding

int fixup_data_afft(struct node *node)
{
	if( node->size < sizeof(struct block_iqdata_float) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_afft)) ;
		return 1 ;
	}
	struct block_iqdata_float *afft = (struct block_iqdata_float *)node->data ;		// hardcoded type
	int nsamples = (node->size)/sizeof(struct block_iqdata_float) ;		// hardcoded type
	for( int loop = 0 ; loop < nsamples ; loop++, afft++ )
	{
		endian_fixup(&(afft->isample),sizeof(afft->isample)) ;
		endian_fixup(&(afft->qsample),sizeof(afft->qsample)) ;
	}
	return 0 ;
}

int dump_block_afft(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_iqdata_float) )		// hardcoded type
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_afft)) ;
		return 1 ;
	}
	if( config->bin_format != BINFORMAT_CVIQ )	// TODO check if the typecasting is needed
	{
		fprintf(stderr,"Cannot handle BINFORMAT %u\n",(uint32_t )config->bin_format) ;		// TODO check if the typecasting is needed
		return 1 ;
	}
	if( (uint32_t )config->bin_type != BINTYPE_FLT4 )	// TODO check typecast
	{
		fprintf(stderr,"Cannot handle BINTYPE %u\n",(uint32_t )config->bin_type) ;
		return 1 ;
	}
	fprintf(outfile,"%s\n",strkey(KEY_afft)) ;
	struct block_iqdata_float *afft = (struct block_iqdata_float *)(node->data) ;		// hardcoded type
	int nsamples = (node->size)/sizeof(struct block_iqdata_float) ;		// hardcoded type
	for( int loop = 0 ; loop < nsamples ; loop++, afft++ )
	{
		double isample = afft->isample ;		// hardcoded type
		double qsample = afft->qsample ;		// hardcoded type
		fprintf(outfile,"%3d % .20lf % .20lf\n",loop,isample,qsample) ;
	}
	fprintf(outfile,"\n") ;
	return 0 ;
}

int read_iqdata_samples(struct block_iqdata_float *, int, struct config *, FILE *) ;	// declared here because it depends on the struct definition		// hardcoded type

int make_node_afft(struct node *list, struct config *config, FILE *fd)		// creates a new node for afft block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_afft ;
	int afft_lines = count_iqdata_lines(fd) ;		// count lines, 1 line per sample (i and q), use this to malloc space for the entire block
	if( afft_lines <= 0 )
	{
		fprintf(stderr,"Error counting lines in '%s' block\n",strkey(KEY_afft)) ;
		return 1 ;
	}
	if( afft_lines % 3 != 0 )
	{
		fprintf(stderr,"Bad number of lines: %d, reading '%s' block. Lines must be a multiple of 3\n",afft_lines,strkey(KEY_afft)) ;
		return 1 ;
	}
	if( Debug ) { fprintf(stderr,"debug: make_node_afft: afft_lines=%u\n",afft_lines) ; }
	int afft_samples = afft_lines ;	// a sample is a line of I,Q values
	size_t afft_size = afft_samples * sizeof(struct block_iqdata_float) ;		// hardcoded type
	if( Debug ) { fprintf(stderr,"debug: make_node_afft: afft_samples=%d afft_size=%zu\n",afft_samples,afft_size) ; }
	struct block_iqdata_float *afft_data = malloc(afft_size) ;		// hardcoded type
	if( afft_data == NULL )
	{
		fprintf(stderr,"Malloc error on '%s' data block\n",strkey(KEY_afft)) ;
		return 1 ;
	}
	memset(afft_data,0,afft_size) ;
	newnode->data = (unsigned char *)afft_data ;
	newnode->size = afft_size ;
	if( read_iqdata_samples(afft_data,afft_samples,config,fd) )	// read lines of i, q values as float and store them in afft_data
	{
		fprintf(stderr,"Error reading '%s' block\n",strkey(KEY_afft)) ;
		return 1 ;
	}
	return 0 ;
}

int count_iqdata_lines(FILE *fd)		// count how many lines of i,q data there are before we hit a blank line (end of block)
{
	char line[SIZE_LINE] ;
	int count = 0 ;
	unsigned long afft_start = ftell(fd) ;
	while( fgets(line,SIZE_LINE,fd) != NULL )
	{
		chomp(line,SIZE_LINE) ;				// remove newline
		if( strlen(line) < 1 )
			break ;
		count++ ;
	}
	fseek(fd,afft_start,SEEK_SET) ;
	return count ;
}

int read_iqdata_samples(struct block_iqdata_float *iqdata, int iqsamples, struct config *config, FILE *fd)		// hardcoded type
{
	char line[SIZE_LINE] ;
	if( (uint32_t )config->bin_format != BINFORMAT_CVIQ )
	{
		fprintf(stderr,"Cannot handle BINFORMAT %u\n",(uint32_t )config->bin_format) ;
		return 1 ;
	}
	if( (uint32_t )config->bin_type != BINTYPE_FLT4 )
	{
		fprintf(stderr,"Cannot handle BINTYPE %u\n",(uint32_t )config->bin_type) ;
		return 1 ;
	}
	unsigned long iqstart = ftell(fd) ;
	for( int sample_count = 0 ; sample_count < iqsamples ; sample_count++, iqdata++ )
	{
		if( fgets(line,SIZE_LINE,fd) == NULL ) return 1 ;
		chomp(line,SIZE_LINE) ;
		if( strlen(line) == 0 ) return 1 ;
		int count ;
		double i ;
		double q ;
		int convert_count = sscanf(line,"%d %lf %lf",&count,&i,&q) ;
		if( convert_count != 3 )
		{
			fprintf(stderr,"Failed to read iqdata %d from line %s\n",sample_count,line) ;
			return 1 ;
		}
		iqdata->isample = (float )i ;			// hardcoded type
		iqdata->qsample = (float )q ;			// hardcoded type
		//if( Debug && (sample_count == 0) ) { fprintf(stderr,"debug: read_iqdata_samples: double i=%lf q=%lf, scalar_one=%lf scalar_two=%lf, factor=%lf int i=%d q=%d\n",i,q,config->scalar_one,config->scalar_two,factor,iqdata->isample,iqdata->qsample) ; }
	}
	//fseek(fd,iqstart,SEEK_SET) ;	// not sure why the file descriptor needs to be rewound, so commented out for now
	return 0 ;
}

int gen_block_afft(struct node *node, FILE *outfile)
{
	struct block_iqdata_float *afft = (struct block_iqdata_float *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	int actual_size = node->size ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	int sample_count = actual_size/sizeof(struct block_iqdata_float) ;
	if( Debug ) { fprintf(stderr,"debug: gen_block_afft: actual size %d, sample_count %d\n",actual_size,sample_count) ; }
	for( int count = 0 ; count < sample_count ; count++, afft++ )
	{
		//if( Debug && (count == 0) ) { fprintf(stderr,"debug: gen_block_afft: sample 0 i=%d q=%d\n",afft->isample,afft->qsample) ; }
		endian_fixup(&(afft->isample),sizeof(afft->isample)) ;
		if( fwrite(&(afft->isample),sizeof(afft->isample),1,outfile) != 1 ) return 1 ;
		endian_fixup(&(afft->qsample),sizeof(afft->qsample)) ;
		if( fwrite(&(afft->qsample),sizeof(afft->qsample),1,outfile) != 1 ) return 1 ;
	}
	return 0 ;
}


int fixup_data_ifft(struct node *node)
{
	if( node->size < sizeof(struct block_iqdata_float) )
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_ifft)) ;
		return 1 ;
	}
	struct block_iqdata_float *ifft = (struct block_iqdata_float *)node->data ;
	int nsamples = (node->size)/sizeof(struct block_iqdata_float) ;
	for( int loop = 0 ; loop < nsamples ; loop++, ifft++ )
	{
		endian_fixup(&(ifft->isample),sizeof(ifft->isample)) ;
		endian_fixup(&(ifft->qsample),sizeof(ifft->qsample)) ;
	}
	return 0 ;
}

int dump_block_ifft(struct node *node, struct config *config, FILE *outfile)
{
	if( node->size < sizeof(struct block_iqdata_float) )		// hardcoded type
	{
		fprintf(stderr,"Block '%s' is truncated\n",strkey(KEY_ifft)) ;
		return 1 ;
	}
	if( (uint32_t )config->bin_format != BINFORMAT_CVIQ )
	{
		fprintf(stderr,"Cannot handle BINFORMAT %u\n",(uint32_t )config->bin_format) ;
		return 1 ;
	}
	if( (uint32_t )config->bin_type != BINTYPE_FLT4 )
	{
		fprintf(stderr,"Cannot handle BINTYPE %u\n",(uint32_t )config->bin_type) ;
		return 1 ;
	}
	fprintf(outfile,"%s\n",strkey(KEY_ifft)) ;
	struct block_iqdata_float *ifft = (struct block_iqdata_float *)(node->data) ;		// hardcoded type
	int nsamples = (node->size)/sizeof(struct block_iqdata_float) ;		// hardcoded type
	for( int loop = 0 ; loop < nsamples ; loop++, ifft++ )
	{
		double isample = ifft->isample ;		// hardcoded type
		double qsample = ifft->qsample ;		// hardcoded type
		fprintf(outfile,"%3d % .16lf % .16lf\n",loop,isample,qsample) ;
	}
	fprintf(outfile,"\n") ;
	return 0 ;
}

int make_node_ifft(struct node *list, struct config *config, FILE *fd)		// creates a new node for ifft block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error on list node\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_ifft ;
	int ifft_lines = count_iqdata_lines(fd) ;		// count lines, 1 line per sample (i and q), use this to malloc space for the entire block
	if( ifft_lines <= 0 )
	{
		fprintf(stderr,"Error counting lines in '%s' block\n",strkey(KEY_ifft)) ;
		return 1 ;
	}
	if( ifft_lines % 3 != 0 )
	{
		fprintf(stderr,"Bad number of lines: %d, reading '%s' block. Lines must be a multiple of 3\n",ifft_lines,strkey(KEY_ifft)) ;
		return 1 ;
	}
	if( Debug ) { fprintf(stderr,"debug: make_node_ifft: ifft_lines=%u\n",ifft_lines) ; }
	int ifft_samples = ifft_lines ;	// a sample is a line of I,Q values
	size_t ifft_size = ifft_samples * sizeof(struct block_iqdata_float) ;		// hardcoded type
	if( Debug ) { fprintf(stderr,"debug: make_node_ifft: ifft_samples=%d ifft_size=%zu\n",ifft_samples,ifft_size) ; }
	struct block_iqdata_float *ifft_data = malloc(ifft_size) ;		// hardcoded type
	if( ifft_data == NULL )
	{
		fprintf(stderr,"Malloc error on '%s' data block\n",strkey(KEY_ifft)) ;
		return 1 ;
	}
	memset(ifft_data,0,ifft_size) ;
	newnode->data = (unsigned char *)ifft_data ;
	newnode->size = ifft_size ;
	if( read_iqdata_samples(ifft_data,ifft_samples,config,fd) )	// read lines of i, q values as float and store them in ifft_data
	{
		fprintf(stderr,"Error reading '%s' block\n",strkey(KEY_ifft)) ;
		return 1 ;
	}
	return 0 ;
}

int gen_block_ifft(struct node *node, FILE *outfile)
{
	struct block_iqdata_float *ifft = (struct block_iqdata_float *)(node->data) ;
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	int actual_size = node->size ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	int sample_count = actual_size/sizeof(struct block_iqdata_float) ;
	if( Debug ) { fprintf(stderr,"debug: gen_block_ifft: actual size %d, sample_count %d\n",actual_size,sample_count) ; }
	for( int count = 0 ; count < sample_count ; count++, ifft++ )
	{
		//if( Debug && (count == 0) ) { fprintf(stderr,"debug: gen_block_ifft: sample 0 i=%d q=%d\n",ifft->isample,ifft->qsample) ; }
		endian_fixup(&(ifft->isample),sizeof(ifft->isample)) ;
		if( fwrite(&(ifft->isample),sizeof(ifft->isample),1,outfile) != 1 ) return 1 ;
		endian_fixup(&(ifft->qsample),sizeof(ifft->qsample)) ;
		if( fwrite(&(ifft->qsample),sizeof(ifft->qsample),1,outfile) != 1 ) return 1 ;
	}
	return 0 ;
}


int fixup_data_end(struct node *node)
{
	return 0 ;
}

int dump_block_end(struct node *node, struct config *config, FILE *outfile)
{
	fprintf(outfile,"%s\n",strkey(KEY_END)) ;
	return 0 ;
}

int make_node_end(struct node *list, struct config *config, FILE *fd)		// creates a new node for end block
{
	struct node *newnode = malloc(sizeof(struct node)) ;
	if( newnode == NULL )
	{
		fprintf(stderr,"Malloc error\n") ;
		return 1 ;
	}
	memset(newnode,0,sizeof(struct node)) ;
	list->next = newnode ;
	newnode->key = KEY_END ;
	newnode->size = 0 ;
	newnode->data = NULL ;
	return 0 ;
}

int gen_block_end(struct node *node, FILE *outfile)
{
	endian_fixup(&(node->key),sizeof(node->key)) ;
	if( fwrite(&(node->key),sizeof(node->key),1,outfile) != 1 ) return 1 ;
	endian_fixup(&(node->size),sizeof(node->size)) ;
	if( fwrite(&(node->size),sizeof(node->size),1,outfile) != 1 ) return 1 ;
	return 0 ;
}


// end of block-specific functions


void free_all_nodes(struct node *list)
{
	while( list != NULL )
	{
		struct node *next = list->next ;
		free(list) ;
		list = next ;
	}
}

void free_all_nodes_and_data(struct node *list)
{
	while( list != NULL )
	{
		if( list->data != NULL )
			free(list->data) ;
		struct node *next = list->next ;
		free(list) ;
		list = next ;
	}
}

//END
