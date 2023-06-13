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

int rewind_FS( File_Src* fs ) {
  int res;
  if ( fs == NULL ) {
    res = 0;
  }
  if ( fs->is_gz ) {
    res = gzrewind( fs->fgz );
  }
  else {
    res = fseek( fs->f , 0, SEEK_SET );
  }
  return res;
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
    fprintf( stderr, "Failed to open %s.\n", name );
    perror("Error");
    return NULL;
  }
  return f;
}



/* UC Santa Cruz (UCSC) Noncommercial License

ACCEPTANCE
In order to get any license under these terms, you must agree to them as both strict obligations and conditions to all your licenses.

COPYRIGHT LICENSE
The licensor grants you a copyright license for the software to do everything you might do with the software that would otherwise infringe the licensor's copyright in it for any permitted purpose.
However, you may only distribute the software according to Distribution License and make changes or new works based on the software according to Changes and New Works License.

DISTRIBUTION LICENSE
The licensor grants you an additional copyright license to distribute copies of the software. Your license to distribute covers distributing the software with changes and new works permitted by Changes and New Works License.

NOTICES
You must ensure that anyone who gets a copy of any part of the software from you also gets a copy of these terms, as well as the following copyright notice:
This software is Copyright ©2020-2022. The Regents of the University of California (“Regents”). All Rights Reserved.

CHANGES AND NEW WORKS LICENSE
The licensor grants you an additional copyright license to make changes and new works based on the software for any permitted purpose.

PATENT LICENSE
The licensor grants you the right to use the software as described in any patent applications or issued patents resulting from UCSC Case Number 2022-808.

NONCOMMERCIAL PURPOSES
Any noncommercial purpose is a permitted purpose.

COMMERCIAL PURPOSES
Contact Innovation Transfer, UC Santa Cruz, innovation@ucsc.edu , https://officeofresearch.ucsc.edu/iatc/ , for any commercial purpose.

PERSONAL USES
Personal use for research, experiment, and testing for the benefit of public knowledge, personal study, private entertainment, hobby projects, amateur pursuits, or religious observance, without any anticipated commercial application, is use for a permitted purpose.

NONCOMMERCIAL ORGANIZATIONS
Use by any charitable organization, educational institution, public research organization, public safety or health organization, environmental protection organization, or government institution is use for a permitted purpose regardless of the source of funding or obligations resulting from the funding.

FAIR USE
You may have "fair use" rights for the software under the law. These terms do not limit them.

NO OTHER RIGHTS
These terms do not allow you to sublicense or transfer any of your licenses to anyone else, or prevent the licensor from granting licenses to anyone else.  These terms do not imply any other licenses.

PATENT DEFENSE
If you make any written claim that the software infringes or contributes to infringement of any patent, all your licenses for the software granted under these terms end immediately. If your company makes such a claim, all your licenses end immediately for work on behalf of your company.

VIOLATIONS
The first time you are notified in writing that you have violated any of these terms, or done anything with the software not covered by your licenses, your licenses can nonetheless continue if you come into full compliance with these terms, and take practical steps to correct past violations, 
within 32 days of receiving notice.  Otherwise, all your licenses end immediately.

NO LIABILITY
As far as the law allows, the software comes as is, without any warranty or condition, and the licensor will not be liable to you for any damages arising out of these terms or the use or nature of the software, under any kind of legal claim.

DEFINITIONS
The "licensor" is Regents, and the "software" is the software the licensor makes available under these terms.
"You" refers to the individual or entity agreeing to these terms.
"Your company" is any legal entity, sole proprietorship, or other kind of organization that you work for, plus all organizations that have control over, are under the control of, or are under common control with that organization.  
"Control" means ownership of substantially all the assets of an entity, or the power to direct its management and policies by vote, contract, or otherwise.  Control can be direct or indirect.
"Your licenses" are all the licenses granted to you for the software under these terms.
"Use" means anything you do with the software requiring one of your licenses. */
