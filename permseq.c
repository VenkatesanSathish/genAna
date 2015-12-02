#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mtwist.h"
#include "randistrs.h"
#include "fact_algo.h"
#include "permalign.h"

int main(int argc, char *argv[])
{
  Taainfo *aa_info_tab = NULL;
  Tlocus *loci_tab = NULL;
  Tlocus *loci_tab_bck = NULL;
  char **seqs = NULL;
  int nb_seq, i, j, iof;
  int *grp_tab = NULL, nb_grp;
  int nb_perm;
  size_t *sssz = NULL, tt, nn, x;
  FILE *ofp1 = NULL, *ofp2 = NULL;
  char ofname1[1024] = { '\0' };
  char ofname2[1024] = { '\0' };

  if (argc < 3)
  {
    fprintf(stderr, "Usage: \n");
    fprintf(stderr, " %s K N file1 file2 ... fileN [with N >= 2]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  aa_info_tab = (Taainfo *)chk_calloc(256, sizeof(Taainfo));
  nb_perm = atoi(argv[1]);
  nb_seq = atoi(argv[2]);
  seqs = (char **)chk_calloc(nb_seq, sizeof(char*));
  sssz = (size_t *)chk_calloc(nb_seq, sizeof(size_t));

  if (nb_seq < 2)
  {
    fprintf(stderr, "Error: Need at least two files\n");
    exit(EXIT_FAILURE);
  }

  iof = 3;
  /* Read each sequence */
  for (i=iof; i<nb_seq+iof; ++i)
  {
    seqs[i - iof] = read_nucs_fasta(argv[i], sssz + i - iof);
  }

  if (argc > iof + nb_seq)
  {
    sprintf(ofname1, "%s_s1", argv[iof + nb_seq]);
    sprintf(ofname2, "%s_s2", argv[iof + nb_seq]);
    ofp1 = fopen(ofname1, "a");
    if (ofp1 == NULL)
    {
      fprintf(stderr, "Error: Could not open file %s\n", ofname1);
      exit(EXIT_FAILURE);
    }
    ofp2 = fopen(ofname2, "a");
    if (ofp2 == NULL)
    {
      fprintf(stderr, "Error: Could not open file %s\n", ofname2);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    ofp1 = stdout;
    ofp2 = stdout;
  }

  /* Generate a seed for the random number generator */
//  mt_goodseed(); /* better but slow */
  mt_seed(); /* faster but not so good */

  /* Choose the smallest size (for safety) */
  tt = sssz[0];
  for (i=0; i<nb_seq; ++i)
  {
    if ((sssz[i] % 3) != 0)
    {
      fprintf(stderr,
       "Warning: length of %s is not multiple of three !\n",
       argv[i+2]);
    }
    if (sssz[i] < tt)
    {
      tt = sssz[i];
    }
  }

  nn = tt / 3;

  loci_tab = (Tlocus *)chk_calloc(2 * nn, sizeof(Tlocus));
  loci_tab_bck = loci_tab + nn;
  nb_grp = read_seqs(seqs, loci_tab, aa_info_tab, nb_seq, nn);
  /* make a copy of the original sequence */
  memcpy(loci_tab_bck, loci_tab, nn * sizeof(Tlocus));

  /* do not need the sequences any more */
  for (i=0; i<nb_seq; ++i)
  {
    free(seqs[i]);
  }
  free(seqs);

  /* sort the sequence along its match order */  
  qsort(loci_tab_bck, nn, sizeof(Tlocus), cmplocus_mkey);
  count_match(loci_tab_bck, aa_info_tab, nn);
  cal_iperm(aa_info_tab);

  /* form the list of groups of tuples to permute 
     but only keeping those that can permute at least twice */
  grp_tab = (int *)chk_calloc(nb_grp, sizeof(int));

  j = 0;
  for (i=0; i<256; ++i)
  {
    if (aa_info_tab[i].nb_mt > 1 && aa_info_tab[i].pn > 1)
    {
      grp_tab[j] = i;
      ++ j;
    }
  }
  nb_grp = j;

  /* get the total number of permutation */
  x = 1;
  for(i=0; i<256; ++i)
  {
    if (aa_info_tab[i].pn > 1)
    {
      x *= aa_info_tab[i].pn;
    }
  }

  /* get back the sequence in original order */
  qsort(loci_tab_bck, nn, sizeof(Tlocus), cmplocus_pos);

  for (i=0; i<nb_perm; ++i)
  {
    gen_rand_permutuple(grp_tab, nb_grp, loci_tab, aa_info_tab, nn);
    output_seqs_ver3(ofp1, loci_tab, 0, nn, i);
    output_seqs_ver3(ofp2, loci_tab, 1, nn, i);
  }
  
  if (ofp1 != stdout && ofp2 != stdout)
  {
    fclose(ofp1);
    fclose(ofp2);
  }

  /* freeing memory */
  free_lt(loci_tab, nn);
  free(sssz);
  for (i=0; i<256; ++i)
  {
    if (aa_info_tab[i].nb_mt > 1)
    {
      free(aa_info_tab[i].mttab);
      free(aa_info_tab[i].top);
    }
  }

  free(aa_info_tab);
  free(grp_tab);

  return EXIT_SUCCESS;
}

