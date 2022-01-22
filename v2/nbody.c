//
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//
typedef float              f32;
typedef double             f64;
typedef unsigned long long u64;

//
typedef struct particle_s {

  f32 *x, *y, *z;
  f32 *vx, *vy, *vz;
  f32 *fx, *fy, *fz;
  u64 n;
  
} particle_t;

#ifdef VERBOSE
//
void printFinal(particle_t p) {
    FILE *fp;
    fp  =fopen("Final_pos2.txt", "w");

    for (u64 i = 0;i < p.n;i++) {
        fprintf(fp,"%lld:{%f %f %f}\n",i,p.x[i],p.y[i],p.z[i]);
    }

    fclose(fp);
}
#endif

//
particle_t alloc(particle_t p) { 
  p.x = malloc(sizeof(f32) * p.n);
  p.y = malloc(sizeof(f32) * p.n);
  p.z = malloc(sizeof(f32) * p.n);

  p.vx = malloc(sizeof(f32) * p.n);
  p.vy = malloc(sizeof(f32) * p.n);
  p.vz = malloc(sizeof(f32) * p.n);

  p.fy = malloc(sizeof(f32) * p.n);
  p.fz = malloc(sizeof(f32) * p.n);
  p.fx = malloc(sizeof(f32) * p.n);

  return p;
}

//
void freeParticle(particle_t p) {
  free(p.x);
  free(p.y);
  free(p.z);

  free(p.vx);
  free(p.vy);
  free(p.vz);

  free(p.fx);
  free(p.fy);
  free(p.fz);
}

//
void init(particle_t p)
{
  for (u64 i = 0; i < p.n; i++)
    {
      //
      u64 r1 = (u64)rand();
      u64 r2 = (u64)rand();
      f32 sign = (r1 > r2) ? 1 : -1;
      
      //
      p.x[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p.y[i] = (f32)rand() / (f32)RAND_MAX;
      p.z[i] = sign * (f32)rand() / (f32)RAND_MAX;

      //
      p.vx[i] = (f32)rand() / (f32)RAND_MAX;
      p.vy[i] = sign * (f32)rand() / (f32)RAND_MAX;
      p.vz[i] = (f32)rand() / (f32)RAND_MAX;

    }
}

//
void move_particles(particle_t p, const f32 dt)
{
  //
  const f32 softening = 1e-20;

  //
  for (u64 i = 0; i < p.n; i++)
    {
      clock_t st = clock();

      p.fx[i] = 0.0;
      p.fy[i] = 0.0;
      p.fz[i] = 0.0;

      f32 fx = 0, fy = 0, fz = 0;


       // 25 floating-point operations
      for (u64 j = i; j < p.n; j++)
      {
          //Newton's law
          const f32 dx = p.x[j] - p.x[i]; //1
          const f32 dy = p.y[j] - p.y[i]; //2
          const f32 dz = p.z[j] - p.z[i]; //3
          const f32 d_2 = (dx * dx) + (dy * dy) + (dz * dz) + softening; //9
          const f32 d_3_over_2 = pow(d_2, 1.5); //10

          //Net force
          fx = dx / d_3_over_2; //11
          fy = dy / d_3_over_2; //12
          fz = dz / d_3_over_2; //13

          p.fx[i] += fx;    //14
          p.fy[i] += fy;    //15
          p.fz[i] += fz;    //16

          p.fx[j] += fx;    //17
          p.fy[j] += fy;    //18
          p.fz[j] += fz;    //19
      }

      //
      p.vx[i] += dt * p.fx[i]; //21
      p.vy[i] += dt * p.fy[i]; //23
      p.vz[i] += dt * p.fz[i]; //25

      clock_t en = clock();

      f32 time_passed = (float)(en - st) / CLOCKS_PER_SEC;

      //if (!(i%100)) printf("Time ellapsed for loop %lld : %f\n", i, time_passed);

    }

  //3 floating-point operations
  for (u64 i = 0; i < p.n; i++)
    {
      p.x[i] += dt * p.vx[i];
      p.y[i] += dt * p.vy[i];
      p.z[i] += dt * p.vz[i];
    }
}

//
int main(int argc, char **argv)
{
  //
  const u64 n = (argc > 1) ? atoll(argv[1]) : 16384;
  const u64 steps= 10;
  const f32 dt = 0.01;

  //
  f64 rate = 0.0, drate = 0.0;

  //Steps to skip for warm up
  const u64 warmup = 3;
  
  //
  particle_t p;
  p.n = n;
  p = alloc(p);

#ifdef VERBOSE
  srand(69);
#endif

  //
  init(p);

  const u64 s = 9 * sizeof(f32) * n + sizeof(u64);
  
  printf("\n\033[1mTotal memory size:\033[0m %llu B, %llu KiB, %llu MiB\n\n", s, s >> 10, s >> 20);
  
  //
  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);
  
  //
  f32 avrg = 0;

  //
  for (u64 i = 0; i < steps; i++)
    {
      //Measure
      const f64 start = omp_get_wtime();

      move_particles(p, dt);

      const f64 end = omp_get_wtime();

      avrg += end-start;

      //Number of interactions/iterations
      const f32 h1 = ( (f32)(n) * (f32)(n - 1) ) / (f32)(2.0);

      //GFLOPS
      const f32 h2 = (25.0 * h1 + 3.0 * (f32)n) * 1e-9;
      
      if (i >= warmup)
	{
	  rate += h2 / (end - start);
	  drate += (h2 * h2) / ((end - start) * (end - start));
	}

      //
      printf("%5llu %10.3e %10.3e %8.1f %s\n",
	     i,
	     (end - start),
	     h1 / (end - start),
	     h2 / (end - start),
	     (i < warmup) ? "*" : "");
      
      fflush(stdout);
    }

  //
  rate /= (f64)(steps - warmup);
  drate = sqrt(drate / (f64)(steps - warmup) - (rate * rate));

  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.1lf +- %.1lf GFLOP/s\033[0m\n",
	 "Average performance:", "", rate, drate);
  printf("-----------------------------------------------------\n");
  printf("Temps moyen d'éxécutions d'une étapes: %f s\n",avrg/(f32)steps);
  printf("-----------------------------------------------------\n");
  
#ifdef VERBOSE
  //
  printFinal(p);
#endif

  //
  freeParticle(p);

  //
  return 0;
}
