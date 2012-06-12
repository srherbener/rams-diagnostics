#include <stdio.h>
#include <malloc/malloc.h>
#include <time.h>
#include <sys/time.h>

#define BIG_NX 3000
//#define BIG_NX 6
#define BIG_NY 200
#define BIG_NZ 65
#define BIG_NITEMS 10
//#define BIG_NITEMS 5000

double my_time();

main(int argc, char *argv[], char **envp)
  {
  float *data;
  FILE *fid;
  int i;
  int buf_size;
  int total_elems;
  double tstart;
  double ttot;

  // fill up the buffer
  total_elems = BIG_NX * BIG_NY * BIG_NZ;
  buf_size = total_elems * sizeof(float);
  data = (float *) malloc(buf_size);
  printf("Total elements = %d, size = %d\n", total_elems, buf_size);
  for (i=0; i<total_elems; i++)
    {
    data[i] = i;
    }

  fid = fopen("test_big_file", "w");

  
  printf("Writing file:\n");
  tstart = my_time();
  for (i=0; i<BIG_NITEMS; i++)
    {
    fwrite(data, total_elems, sizeof(float), fid);
    }
  ttot = my_time() - tstart;

  fclose(fid);

  printf("\n");
  printf("  Total elapsed time (s): %f, Average time per write (s): %f\n", ttot, (ttot/(double)BIG_NITEMS));
  }

double my_time()
  {
  struct timeval tstruct; // two elements: tv_sec, tv_usec
  double dbl_time;
 
  gettimeofday(&tstruct,NULL);
  dbl_time = (double)tstruct.tv_sec + ((double)tstruct.tv_usec * 1.0e-6);

  return dbl_time;
  }
