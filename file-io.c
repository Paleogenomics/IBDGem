#include "file-io.h"

/* is_gz
   Takes filename as argument
   Returns true IFF the argument filename ends in .gz
   Returns: 1 => name ends in .gz
            0 => name does not end in .gz
*/
int is_gz( const char* fn ) {
  size_t fn_len;
  fn_len = strlen( fn );
  if ( (fn[fn_len-3] == '.') &&
       (fn[fn_len-2] == 'g') &&
       (fn[fn_len-1] == 'z') ) {
    return 1;
  }
  return 0;
}

char* get_line_FS( File_Src* fs, char* line ) {
  if (fs->is_gz) {
    return gzgets( fs->fgz, line, MAX_LINE_LEN );
  }
  else {
    return fgets( line, MAX_LINE_LEN, fs->f );
  }
}

File_Src* init_FS( const char* fn ) {
  File_Src* fs;
  if ( fn == NULL ) {
    return NULL;
  }

  fs = (File_Src*)malloc(sizeof( File_Src ));
  strcpy( fs->fn, fn );

  if ( is_gz( fn ) ) {
    fs->is_gz = 1;
    fs->fgz = gzopen( fs->fn, "r" );
    if ( fs->fgz == NULL ) {
      free( fs );
      return NULL;
    }
  }
  else {
    fs->is_gz = 0;
    fs->f = fileOpen( fs->fn, "r" );
    if ( fs->f == NULL ) {
      free( fs );
      return NULL;
    }
  }
  return fs;
}

File_Src* reset_FS( File_Src* fs ) {
  char fn[MAX_FN_LEN + 1];
  
  if ( fs == NULL ) {
    return NULL;
  }

  strcpy( fn, fs->fn );
  destroy_FS( fs );
  return init_FS( fn );
}

int destroy_FS( File_Src* fs ) {
  if ( fs == NULL ) {
    return 0;
  }
  if ( fs->is_gz ) {
    gzclose( fs->fgz );
  }
  else {
    fclose( fs->f );
  }
  free( fs );
  return 0;
}

/* fileOpen */
FILE* fileOpen( const char* name, char access_mode[] ) {
  FILE* f;
  f = fopen( name, access_mode );
  if ( f == NULL ) {
    fprintf( stderr, "%s\n", name );
    perror( "Cannot open file" );
    return NULL;
  }
  return f;
}
