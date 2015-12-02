#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* ------------------------------------ */
/* ********** BEGIN OF HEADER STATMENTS */
/* ------------------------------------ */

/* Structures for multi-sequence table */

struct _tmb_ {
  char *name; /* pointer to the name of the sequence (aftere '>') */
  char *seq;  /* pointer to the sequence content */
  size_t size; /* size of the sequence */
};

typedef struct {
  char *data; /* sequence file data */
  size_t len; /* size of the sequence file data */
  struct _tmb_ *tab; /* table of name/sequence/size */
  size_t nb; /* number of sequence */

} Tmultiseq;

/* Functions declarations */
void *chk_calloc(size_t nmemb, size_t size);
size_t get_filesize (const char *filename);
int read_multiseq(char *filename, Tmultiseq *mseqdata);
size_t fline_count(char *filename, size_t *m);
static int cmpname(const void *t1, const void *t2);

/* ---------------------------------- */
/* END OF HEADER STATMENTS ********** */
/* ---------------------------------- */

/**
   Safe function to allocate memory
*/
void *chk_calloc(size_t nmemb, size_t size)
{
  void *mem = calloc(nmemb, size);

  if (mem == NULL)
  {
    fprintf(stderr, "Fatal error: failed to allocate %zu bytes of memory",
             nmemb * size);
    exit(EXIT_FAILURE);
  }

  /* Fill it with zeros, libc malloc bug workaround */
  return memset(mem, 0, nmemb * size);
}

/**
   Function to obtain the size of the file
*/
size_t get_filesize (const char *filename)
{
  int rf;
  struct stat finfo;

  rf = stat(filename, &finfo);

  if (rf != 0)
  {
    fprintf(stderr, "Error: could not obtain the size of %s", filename);
    exit(EXIT_FAILURE);
  }

  /* Check that file is not a directory */
  if (S_ISREG(finfo.st_mode))
  {
    /* Check for the zero-size file */
    if (finfo.st_size == 0)
    {
      fprintf(stderr, "Error: the file %s is zero-size", filename);
      exit(EXIT_FAILURE);
    }
    return finfo.st_size;
  }
  else
  {
    fprintf(stderr, "Error: %s is not a regular file", filename);
    exit(EXIT_FAILURE);
  }
}

/**
   Function to read a multi-sequence file
*/
int read_multiseq(char *filename, Tmultiseq *mseqdata)
{
  FILE *fp = NULL;
  int r;
  char *ptr, cdeb, c;
  size_t file_sz, i, k, m;

  /* Obtain the size of the file*/
  file_sz = get_filesize(filename);

  /* Read open the file */
  if ((fp = fopen(filename ,"r")) == NULL)
  {
    fprintf(stderr, "Error: could not open FASTA file %s", filename);
    exit(EXIT_FAILURE);
  }

  /* Allocate memory */
  mseqdata->data = (char *)chk_calloc(file_sz + 1, sizeof(char));

  /* Read the whole file */
  ptr = mseqdata->data;
  mseqdata->nb = 0;
  mseqdata->len = 0;
  cdeb = c = '\0';
  while (file_sz > 0)
  {
    r = fgetc(fp);
    /* stop reading if end-of-file */
    if (r == EOF)
    {
      *ptr = '\0';
      break;
    }

    if (c == '\n' || c == '\0')
    {
      cdeb = (char)r;
    }

    c = (char)r;
    /* raise char-flag when start-of-sequence-name and
       increment the sequence number counter */
    if (c == '>')
    {
      ++ mseqdata->nb;
      if (ptr != mseqdata->data)
      {
        *ptr = '\0';
        ++ mseqdata->len;
        ++ ptr;
        -- file_sz;
      }
    }

    /* replace 'end of line' by 'end of string' for the sequence name */
    if (c == '\n' && cdeb == '>')
    {
      c = '\0';
    }

    /* store the char if not 'end of line' */
    if (c != '\n')
    {
      *ptr = c;
      ++ mseqdata->len;
      ++ ptr;
      -- file_sz;
    }
  }

  /* close the file */
  fclose(fp);

  /* allocate memory for the multiseq table */
  mseqdata->tab = (struct _tmb_ *)chk_calloc(mseqdata->nb,
                                             sizeof(struct _tmb_));

  /* assign to the table element pointer their values */
  k = i = 0;
  while (i < mseqdata->len)
  {
    switch (mseqdata->data[i]) {
      case '>':
        mseqdata->tab[k].name = &(mseqdata->data[i + 1]);
        break;
      case '\0':
        mseqdata->tab[k].seq = &(mseqdata->data[i + 1]);
        m = i + 1;
        while (mseqdata->data[m] != '\0') { ++ m;};
        mseqdata->tab[k].size = m - (i + 1);
        i = m;
        ++ k;
        break;
      default:
        break;
    }
    ++ i;
  }

  return EXIT_SUCCESS;
}

/**
   Free the multi-sequence file structure
*/
void free_multiseq(Tmultiseq *mseqdata)
{
  free(mseqdata->data);
  free(mseqdata->tab);
}

/**
   Function to count the number of line in a file
   (same as the UNIX command line wc -l)
   Give also the size of the longest line in m
*/
size_t fline_count(char *filename, size_t *m)
{
  FILE *fp = NULL;
  int r;
  size_t n, k;
 
  if ((fp = fopen(filename ,"r")) == NULL)
  {
    fprintf(stderr, "Error: could not open FASTA file %s", filename);
    exit(EXIT_FAILURE);
  }
  
  *m = n = k = 0;
  while ((r = fgetc(fp)) != EOF)
  {
    ++ k;
    if ((char)r == '\n')
    {
      ++ n;
      if (k > *m)
      {
        *m = k;
      }
      k = 0;
    }
  }
  fclose(fp);

  return n;
}

/**
   Function to compare two sequence names (used in qsort and bsearch)
*/

static int cmpname(const void *t1, const void *t2)
{
  struct _tmb_ *tt1 = (struct _tmb_ *)t1;
  struct _tmb_ *tt2 = (struct _tmb_ *)t2;

  return strcmp(tt1->name, tt2->name);
}


/**
    -----------------------------------------
    -------------- MAIN PROGRAM -------------
    -----------------------------------------

    argv[1] : name of the sequence file
    argv[2] : name of the selected range file
              having the following line format
              <seq_name> <start_idx> <stop_idx>
*/

int main (int argc, char *argv[])
{
  FILE *fp = NULL;
  char *lline = NULL, *lname = NULL;
  Tmultiseq msdat;
  size_t nl, nsel, maxl, lc;
  size_t *tab_start = NULL, *tab_stop = NULL;
  struct _tmb_ memb, **tab_nss = NULL, *ppos = NULL;
  size_t i, j;

  if (argc != 3)
  {
     fputs("Usage:\n", stderr);
     fprintf(stderr,
        "  %s <multi-sequence file name> <ranges selection file name>\n",
        argv[0]);
     exit(EXIT_FAILURE);
  }

  /* Read the sequences file */
  read_multiseq(argv[1], &msdat);

  /* Sort the table of name/seq/size along the sequence name */
  qsort(msdat.tab, msdat.nb, sizeof(struct _tmb_), cmpname);

  /* Get the number of line in the selection file */
  nl = fline_count(argv[2], &maxl);

  /* Allocate memory for the selection table */
  tab_start = (size_t *)chk_calloc(nl, sizeof(size_t));  
  tab_stop = (size_t *)chk_calloc(nl, sizeof(size_t));  
  tab_nss = (struct _tmb_ **)chk_calloc(nl, sizeof(struct _tmb_ *));  
  /* Allocate temporary memory */
  maxl += 3; /* 3 -> '\0' & '\n' */
  lname = (char *)chk_calloc(maxl, sizeof(char));
  lline = (char *)chk_calloc(maxl, sizeof(char));

  /* Read each line of the selection file */
  nsel = lc = 0;
  fp = fopen(argv[2], "r");
  while (fgets(lline, maxl, fp) != NULL)
  {
    ++ lc;
    sscanf(lline,"%s %zd %zd\n", lname, &i, &j);
    if (i <= j)
    {
      /* binary search of the sequence in the table from its name */
      memb.name = lname;
      ppos = bsearch(&memb, msdat.tab, msdat.nb,
                     sizeof(struct _tmb_), cmpname);
      if (ppos == NULL)
      {
        fprintf(stderr, "Warning: Unknown sequence name [line %zd]: %s\n",
                lc, lname);
      }
      else
      {
        tab_start[nsel] = i;
        tab_stop[nsel] = j;
        tab_nss[nsel] = ppos;
        ++nsel;
      }
    }
    else
    {
      fprintf(stderr,
                "Warning: inconsistent start=%zd > stop=%zd [line %zd]\n",
                 i, j, lc);
    }
  }

  /* free temporary memory ... */
  free(lname);
  free(lline);
  /* ... and close file */
  fclose(fp);

  /* Print the selected sequences (in selection order) */
  {
    size_t sssz; /* size of the sub-sequence */
    char *ssp; /* pointer to the beginning of the sub-sequence */

    for (i=0; i<nsel; ++i)
    {
      /* print the sequence name */
      printf(">%s\n", tab_nss[i]->name);

      /* print the sequence from start to stop */
      sssz = tab_stop[i] - tab_start[i] + 1;
         /* IMPORTANT: if in the selection file, the sequence
             index is 1-2-3...-n  and not 0-1-2...-n
             therefore replace the following line by
             ssp = tab_nss[i].seq + tab_start[i] - 1; */ 
      ssp = tab_nss[i]->seq + tab_start[i];
      /* For safety test out-of-range indices, IMPORTANT holds here also */
      if (tab_start[i] > tab_nss[i]->size || tab_stop[i] > tab_nss[i]->size)
      {
        fprintf(stderr, "Warning: start=%zd or stop=%zd is out of range\n",
                     tab_start[i], tab_stop[i]);
      }
      else
      {
        for (j=0; j<sssz; ++j)
        {
          putchar(ssp[j]);
        }
        putchar('\n');
      }
    }
  }

  /* free all allocated memory */
  free_multiseq(&msdat);
  free(tab_start);
  free(tab_stop);
  free(tab_nss);

  return EXIT_SUCCESS;
}

