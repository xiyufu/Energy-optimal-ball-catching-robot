/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) solver_ ## ID
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef casadi_real
#define casadi_real double
#endif

#define to_double(x) (double) x
#define to_int(x) (int) x
#define CASADI_CAST(x,y) (x) y

/* Pre-c99 compatibility */
#if __STDC_VERSION__ < 199901L
  #define fmin CASADI_PREFIX(fmin)
  casadi_real fmin(casadi_real x, casadi_real y) { return x<y ? x : y;}
  #define fmax CASADI_PREFIX(fmax)
  casadi_real fmax(casadi_real x, casadi_real y) { return x>y ? x : y;}
#endif

/* CasADi extensions */
#define sq CASADI_PREFIX(sq)
casadi_real sq(casadi_real x) { return x*x;}
#define sign CASADI_PREFIX(sign)
casadi_real CASADI_PREFIX(sign)(casadi_real x) { return x<0 ? -1 : x>0 ? 1 : x;}
#define twice CASADI_PREFIX(twice)
casadi_real twice(casadi_real x) { return x+x;}
#define if_else CASADI_PREFIX(if_else)
casadi_real if_else(casadi_real c, casadi_real x, casadi_real y) { return c!=0 ? x : y;}

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)

/* Printing routine */
#define PRINTF printf

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

static const int casadi_s0[22] = {18, 1, 0, 18, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
static const int casadi_s1[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};
static const int casadi_s2[4] = {0, 1, 0, 0};

/* integrator:(x0[18],p[18],z0[6],rx0[0],rp[0],rz0[0])->(xf[18],qf[0],zf[6],rxf[0],rqf[0],rzf[0]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem) {
  #error Code generation not supported for IdasInterface
  return 0;
}

CASADI_SYMBOL_EXPORT int integrator(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void integrator_incref(void) {
}

CASADI_SYMBOL_EXPORT void integrator_decref(void) {
}

CASADI_SYMBOL_EXPORT int integrator_n_in(void) { return 6;}

CASADI_SYMBOL_EXPORT int integrator_n_out(void) { return 6;}

CASADI_SYMBOL_EXPORT const char* integrator_name_in(int i){
  switch (i) {
    case 0: return "x0";
    case 1: return "p";
    case 2: return "z0";
    case 3: return "rx0";
    case 4: return "rp";
    case 5: return "rz0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* integrator_name_out(int i){
  switch (i) {
    case 0: return "xf";
    case 1: return "qf";
    case 2: return "zf";
    case 3: return "rxf";
    case 4: return "rqf";
    case 5: return "rzf";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* integrator_sparsity_in(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    case 4: return casadi_s2;
    case 5: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* integrator_sparsity_out(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s2;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    case 4: return casadi_s2;
    case 5: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int integrator_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 13;
  if (sz_res) *sz_res = 12;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 1101;
  return 0;
}

int main_integrator(int argc, char* argv[]) {
  int *iw = 0;
  casadi_real w[1167];
  const casadi_real* arg[13] = {w+0, w+18, w+36, w+42, w+42, w+42};
  casadi_real* res[12] = {w+42, w+60, w+60, w+66, w+66, w+66};
  int j;
  casadi_real* a = w;
  for (j=0; j<42; ++j) scanf("%lf", a++);
  int flag = integrator(arg, res, iw, w+66, 0);
  if (flag) return flag;
  const casadi_real* r = w+42;
  for (j=0; j<24; ++j) PRINTF("%g ", *r++);
  PRINTF("\n");
  return 0;
}


int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "integrator")==0) {
    return main_integrator(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'integrator'\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif
