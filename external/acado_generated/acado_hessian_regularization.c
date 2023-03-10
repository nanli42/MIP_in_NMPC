/*
 *    This file was auto-generated using the ACADO Toolkit.
 *    
 *    While ACADO Toolkit is free software released under the terms of
 *    the GNU Lesser General Public License (LGPL), the generated code
 *    as such remains the property of the user who used ACADO Toolkit
 *    to generate this code. In particular, user dependent data of the code
 *    do not inherit the GNU LGPL license. On the other hand, parts of the
 *    generated code that are a direct copy of source code from the
 *    ACADO Toolkit or the software tools it is based on, remain, as derived
 *    work, automatically covered by the LGPL license.
 *    
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *    
 */


#include "acado_common.h"

#define ACADO_EPS 1e-12
#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#include <math.h>

#define DIM 11

static real_t hypot2(real_t x, real_t y) {
  return sqrt(x*x+y*y);
}

/* Symmetric Householder reduction to tridiagonal form. */

static void acado_tred2(real_t *V, real_t *d, real_t *e) {

/* This is derived from the Algol procedures tred2 by
   Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   Fortran subroutine in EISPACK. */

  int i,j,k;
  real_t f,g,h,hh;
  for (j = 0; j < DIM; j++) {
    d[j] = V[(DIM-1)*DIM+j];
  }

/* Householder reduction to tridiagonal form. */

  for (i = DIM-1; i > 0; i--) {

    /* Scale to avoid under/overflow. */

    real_t scale = 0.0;
    real_t h = 0.0;
    for (k = 0; k < i; k++) {
      scale = scale + fabs(d[k]);
    }
    if (scale == 0.0) {
      e[i] = d[i-1];
      for (j = 0; j < i; j++) {
        d[j] = V[(i-1)*DIM+j];
        V[i*DIM+j] = 0.0;
        V[j*DIM+i] = 0.0;
      }
    } else {

      /* Generate Householder vector. */

      for (k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      f = d[i-1];
      g = sqrt(h);
      if (f > 0) {
        g = -g;
      }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (j = 0; j < i; j++) {
        e[j] = 0.0;
      }

      /* Apply similarity transformation to remaining columns. */

      for (j = 0; j < i; j++) {
        f = d[j];
        V[j*DIM+i] = f;
        g = e[j] + V[j*DIM+j] * f;
        for (k = j+1; k <= i-1; k++) {
          g += V[k*DIM+j] * d[k];
          e[k] += V[k*DIM+j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      hh = f / (h + h);
      for (j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (k = j; k <= i-1; k++) {
          V[k*DIM+j] -= (f * e[k] + g * d[k]);
        }
        d[j] = V[(i-1)*DIM+j];
        V[i*DIM+j] = 0.0;
      }
    }
    d[i] = h;
  }

  /* Accumulate transformations. */

  for (i = 0; i < DIM-1; i++) {
    V[(DIM-1)*DIM+i] = V[i*DIM+i];
    V[i*DIM+i] = 1.0;
    h = d[i+1];
    if (h != 0.0) {
      for (k = 0; k <= i; k++) {
        d[k] = V[k*DIM+i+1] / h;
      }
      for (j = 0; j <= i; j++) {
        g = 0.0;
        for (k = 0; k <= i; k++) {
          g += V[k*DIM+i+1] * V[k*DIM+j];
        }
        for (k = 0; k <= i; k++) {
          V[k*DIM+j] -= g * d[k];
        }
      }
    }
    for (k = 0; k <= i; k++) {
      V[k*DIM+i+1] = 0.0;
    }
  }
  for (j = 0; j < DIM; j++) {
    d[j] = V[(DIM-1)*DIM+j];
    V[(DIM-1)*DIM+j] = 0.0;
  }
  V[(DIM-1)*DIM+DIM-1] = 1.0;
  e[0] = 0.0;
} 

/* Symmetric tridiagonal QL algorithm. */

static void acado_tql2(real_t *V, real_t *d, real_t *e) {

/*  This is derived from the Algol procedures tql2, by
   Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
   Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
   Fortran subroutine in EISPACK. */

  int i,j,m,l,k;
  real_t g,p,r,dl1,h,f,tst1,eps;
  real_t c,c2,c3,el1,s,s2;

  for (i = 1; i < DIM; i++) {
    e[i-1] = e[i];
  }
  e[DIM-1] = 0.0;

  f = 0.0;
  tst1 = 0.0;
  eps = pow(2.0,-52.0);
  for (l = 0; l < DIM; l++) {

    /* Find small subdiagonal element */

    tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
    m = l;
    while (m < DIM) {
      if (fabs(e[m]) <= eps*tst1) {
        break;
      }
      m++;
    }

    /* If m == l, d[l] is an eigenvalue,
       otherwise, iterate. */

    if (m > l) {
      int iter = 0;
      do {
        iter = iter + 1;
        /* Compute implicit shift */

        g = d[l];
        p = (d[l+1] - g) / (2.0 * e[l]);
        r = hypot2(p,1.0);
        if (p < 0) {
          r = -r;
        }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        dl1 = d[l+1];
        h = g - d[l];
        for (i = l+2; i < DIM; i++) {
          d[i] -= h;
        }
        f = f + h;

        /* Implicit QL transformation. */

        p = d[m];
        c = 1.0;
        c2 = c;
        c3 = c;
        el1 = e[l+1];
        s = 0.0;
        s2 = 0.0;
        for (i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = hypot2(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          /* Accumulate transformation. */

          for (k = 0; k < DIM; k++) {
            h = V[k*DIM+i+1];
            V[k*DIM+i+1] = s * V[k*DIM+i] + c * h;
            V[k*DIM+i] = c * V[k*DIM+i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        /* Check for convergence. */

      } while (fabs(e[l]) > eps*tst1 && iter < 20);  /* (Check iteration count here.) */
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
}

void acado_eigen_decomposition(real_t *A, real_t *V, real_t *d) {
  int i,j;
  real_t e[DIM];
  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      V[i*DIM+j] = A[i*DIM+j];
    }
  }
  acado_tred2(V, d, e);
  acado_tql2(V, d, e);
}

void acado_reconstruct_A(real_t *A, real_t *V, real_t *d) {
  int i, j, k;
  for( i = 0; i < DIM; i++ ) {
	for( j = 0; j <= i; j++ ) {
	  A[i*DIM+j] = 0.0;
	  for( k = 0; k < DIM; k++ ) {
	    A[i*DIM+j] += V[i*DIM+k]*d[k]*V[j*DIM+k];
	  }
	  A[j*DIM+i] = A[i*DIM+j];
	}
  }
}

/* cutting regularization */
/*void acado_regularize(real_t *A) {
  int i;
  real_t V[DIM*DIM];
  real_t d[DIM];
  
  acado_eigen_decomposition(A, V, d);
  
  for (i = 0; i < DIM; i++) {
    if( d[i] <= ACADO_EPS ) d[i] = ACADO_EPS;
  }
  
  acado_reconstruct_A(A, V, d);
}*/

/* mirroring regularization */
void acado_regularize(real_t *A) {
  int i;
  real_t V[DIM*DIM];
  real_t d[DIM];
  
  acado_eigen_decomposition(A, V, d);
  
  for (i = 0; i < DIM; i++) {
    if( d[i] >= -ACADO_EPS && d[i] <= ACADO_EPS ) d[i] = ACADO_EPS;
    else if( d[i] < 0 ) d[i] = -d[i];
  }
  
  acado_reconstruct_A(A, V, d);
}
