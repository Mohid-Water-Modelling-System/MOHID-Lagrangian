/*
 Copyright 2004, Magnus Hagdorn
 
 This file is part of proj4.
 
 proj4 is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 proj4 is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with proj4; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <projects.h>
#include "cfortran.h"

/*
 * error string
 */
FCALLSCFUN1(STRING,pj_strerrno,PRJF_STRERRNO,prjf_strerrno,INT);

/*
 * initialise projection structure
 */
#define prjf_init_STRV_A4 NUM_ELEM_ARG(2)
int cfort_pj_init(int *prj, int nargs, int largs, char *args)
{
  int i;
  char **cargs;

  /* copying arguments */
  cargs = (char **) malloc(sizeof(char *)*nargs);
  for (i=0;i<nargs;i++)
    cargs[i] = (args+i*(largs+1));

  *prj = (int) pj_init(nargs,cargs);
  free(cargs);
  if (!*prj)
    return pj_errno;
  else
    return 0;
}
FCALLSCFUN4(INT,cfort_pj_init,PRJF_INIT,prjf_init, PINT, INT, INT, STRINGV);

/*
 * free projection structure
 */
int cfort_pj_free(int prj)
{
  pj_free((PJ *) prj);
  return 0;
}
FCALLSCFUN1(INT,cfort_pj_free, PRJF_FREE, prjf_free, INT);

/*
 * forward transform
 */
int cfort_pj_fwd(int prj, double lam, double phi, double *x, double *y)
{
  projUV data;

  data.u = lam* DEG_TO_RAD;
  data.v = phi* DEG_TO_RAD;
  
  data = pj_fwd(data, (PJ *) prj);

  *x = data.u;
  *y = data.v;

  if (data.u==HUGE_VAL && data.v==HUGE_VAL)
    return  pj_errno;
  else
    return 0;
}
FCALLSCFUN5(INT,cfort_pj_fwd,PRJF_FWD,prjf_fwd,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE);

/*
 * inverse transform
 */
int cfort_pj_inv(int prj, double x, double y, double *lam, double *phi)
{
  projUV data;
  
  data.u = x;
  data.v = y;
  
  data = pj_inv(data, (PJ *) prj);

  *lam = data.u * RAD_TO_DEG;
  *phi = data.v * RAD_TO_DEG;

  if (data.u==HUGE_VAL && data.v==HUGE_VAL)
    return  pj_errno;
  else
    return 0;
}
FCALLSCFUN5(INT,cfort_pj_inv,PRJF_INV,prjf_inv,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE);
