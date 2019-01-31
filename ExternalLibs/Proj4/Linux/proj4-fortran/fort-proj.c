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


/*
 * error string
 */
int cfort_pj_strerrno(int err, char *err_msg)
{
    char* msg = pj_strerrno(err);
    if (msg != NULL)
        strcpy(err_msg, msg);
    return pj_ctx_get_errno(pj_get_default_ctx());
}

/*
 * initialise projection structure
 */
int cfort_pj_init_plus(projPJ *prjdefn, char *args)
{
    *prjdefn = pj_init_plus(args);
	return pj_ctx_get_errno(pj_get_default_ctx());
}

/*
 * free projection structure
 */
int cfort_pj_free(projPJ *prjdefn)
{
    pj_free(*prjdefn);
	return pj_ctx_get_errno(pj_get_default_ctx());
}

/*
 * transform projection
 */
int cfort_pj_transform(projPJ *srcdefn, projPJ *dstdefn,
                         long point_count, int point_offset,
                         double *x, double *y, double *z)
{
    return  pj_transform(*srcdefn, *dstdefn,
                         point_count, point_offset, x, y, z);
}
