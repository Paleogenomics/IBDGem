#ifndef FILEIO_H
#define FILEIO_H
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>

#define MAX_LINE_LEN 12288
#define MAX_FN_LEN 1024

/* File_Src is a structure that holds information about
   a file. The file can be gz compressed or not. Each
   File_Src knows its type using the is_gz boolean.
*/
typedef struct file_src {
  char fn[MAX_FN_LEN+1];
  int is_gz;
  gzFile fgz;
  FILE* f;
} File_Src;
  
/* is_gz
   Takes filename as argument
   Returns true IFF the argument filename ends in .gz
   Returns: 1 => name ends in .gz
            0 => name does not end in .gz
*/
int is_gz( const char* fn );

/* get_line_FS
   Get the next line from the argument File_Src*
   and put it in the char* line.
   Returns the return value given by the appropriate call
   to fgets or gzgets */
char* get_line_FS( File_Src* fs, char* line );

File_Src* init_FS( const char* fn );

File_Src* reset_FS( File_Src* fs );

int destroy_FS( File_Src* fs );

FILE* fileOpen( const char* name, char access_mode[] );
#endif
