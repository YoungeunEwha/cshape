#include "crust.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
//headers from qhull
#include "qhull.h"
#include "poly.h"
#include "qset.h"
#include "geom.h"

/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

/*-------------------------------------------------------------------*/
//
// checking each vertex is on the convexhull and computing their normal 

void compute_normals(bool * is_on_convexhull, coordT * normals)
{
	tVertex  ptr_v;
	tVertex * all_v = NULL;
	int vsize = 0;
	int id = 0;

	//global varibles for qhull
	static char * options = (char *)"qhull QJ Pp";
	int curlong, totlong;
	coordT * pt = NULL;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	int vid = 0, verid[3];
	tsFace sface;
	double length;

	//count number of points
	ptr_v = vertices;
	do {
		vsize++;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//allocate memory
	pt = (coordT*)calloc(vsize * 3, sizeof(coordT)); //each point will have three coord
	all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
	assert(pt && all_v);

	vid = 0;
	//copy points
	ptr_v = vertices;
	do {
		pt[id++] = ptr_v->v[0];
		pt[id++] = ptr_v->v[1];
		pt[id++] = ptr_v->v[2];
		all_v[ptr_v->vnum] = ptr_v;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	// initialzation
	for (id = 0; id < vsize; id++)
	{
		is_on_convexhull[id] = false;
		normals[3*id] = 0.0;
		normals[3 * id + 1] = 0.0;
		normals[3 * id+2] = 0.0;
	}

	//using qhull

	qh_init_A(stdin, stdout, stderr, 0, NULL);
	qh_initflags(options);
	qh_init_B(pt, vsize, 3, false);
	qh_qhull();
	qh_check_output();
	qh_produce_output();

	//loop through all faces
	FORALLfacets
	{
		//get vertices of facet
		//loop through each vertex
		vid = 0;
		FOREACHvertex_(facet->vertices)
		{
			is_on_convexhull[qh_pointid(vertex->point)] = true;
			verid[vid++] = qh_pointid(vertex->point);
		}
		for (id = 0; id < 3; id++)
		{
			normals[verid[0] * 3 + id] += facet->normal[id];
			normals[verid[1] * 3 + id] += facet->normal[id];
			normals[verid[2] * 3 + id] += facet->normal[id];
		}
	}

		//normalize
	for (id = 0; id < vsize; id++)
	{
		if (is_on_convexhull[id])
		{
			length = sqrt(normals[id * 3] * normals[id * 3] + normals[id * 3 + 1] * normals[id * 3 + 1] + normals[id * 3 + 2] * normals[id * 3 + 2]);
			normals[id * 3] /= length;
			normals[id * 3 + 1] /= length;
			normals[id * 3 + 2] /= length;
		}
	}

	//not used
	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;

	//free mem
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
}

void remove_nonsamples(int psize, coordT * all_poles) //sample size
{
	tVertex  ptr_v;
	tVertex * all_v = NULL;
	int vsize = 0, fsize = 0, esize = 0, tsize = 0;
	int id = 0;

	//global varibles for qhull
	static char * options = (char *)"delaunay QJ Pp";
	int curlong, totlong;
	coordT * pt = NULL;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	//facetT *neighbor, **neighborp;
	int vid = 0;
	tsFace sface;
	tsVertex svertex1, svertex2, svertex3, svertex4;
	double volume;
	int vids[4];


	//count number of points
	ptr_v = vertices;
	do {
		vsize++;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//allocate memory
	pt = (coordT*)calloc((vsize + psize) * 4, sizeof(coordT)); //each point will have three coord
	all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
	assert(pt && all_v);

	//
	// create points in 4D (x,y,z,x^2+y^2+z^2)
	//

	ptr_v = vertices;
	do {
		pt[id++] = ptr_v->v[0];
		pt[id++] = ptr_v->v[1];
		pt[id++] = ptr_v->v[2];
		pt[id++] = ptr_v->v[0] * ptr_v->v[0] + ptr_v->v[1] * ptr_v->v[1] + ptr_v->v[2] * ptr_v->v[2];
		all_v[ptr_v->vnum] = ptr_v;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	for (id = 0; id < psize; id++)
	{
		pt[4 * (vsize + id)] = all_poles[3 * id];
		pt[4 * (vsize + id) + 1] = all_poles[3 * id + 1];
		pt[4 * (vsize + id) + 2] = all_poles[3 * id + 2];
		pt[4 * (vsize + id) + 3] = all_poles[3 * id] * all_poles[3 * id] + all_poles[3 * id + 1] * all_poles[3 * id + 1] + all_poles[3 * id + 2] * all_poles[3 * id + 2];
	}

	//
	// compute convex hull in 4D by calling qhull
	// use flags: static char * options=(char *)"delaunay QJ Pp";
	//

	qh_init_A(stdin, stdout, stderr, 0, NULL);

	qh DELAUNAY = True;     /* 'd'   */
	//qh SCALElast = True;    /* 'Qbb' */
	//qh KEEPcoplanar = True; /* 'Qc', to keep coplanars in 'p' */

	qh_initflags(options);
	qh_init_B(pt, vsize + psize, 4, false);
	qh_qhull();
	qh_check_output();


	//loop through all faces
	FORALLfacets
	{

		//get vertices of tetra
		//loop through each vertex
		vid = 0;
		FOREACHvertex_(facet->vertices)
		{
			//get the id of the vertex
			vids[vid++] = qh_pointid(vertex->point);
		}

		// compute the volume of the tetrahedron

		for (id = 0; id<3; id++)
		{
			svertex1.v[id] = pt[4 * vids[0] + id];
			svertex2.v[id] = pt[4 * vids[1] + id];
			svertex3.v[id] = pt[4 * vids[2] + id];
			svertex4.v[id] = pt[4 * vids[3] + id];
		}
		sface.vertex[0] = &svertex1;
		sface.vertex[1] = &svertex2;
		sface.vertex[2] = &svertex3;
		volume = Volumei(&sface, &svertex4);
		if (!facet->upperdelaunay && abs(volume) > 1.0e-5) //remove volume zero tetras
		{

			if (vids[0] < vsize && vids[1] < vsize && vids[2] < vsize)
			{

				MakeFace(all_v[vids[0]], all_v[vids[1]], all_v[vids[2]], NULL);
			}
			if (vids[0] < vsize && vids[1] < vsize && vids[3] < vsize)
			{
				MakeFace(all_v[vids[0]], all_v[vids[1]], all_v[vids[3]], NULL);
			}
			if (vids[1] < vsize && vids[2] < vsize && vids[3] < vsize)
			{
				MakeFace(all_v[vids[1]], all_v[vids[2]], all_v[vids[3]], NULL);
			}
			if (vids[2] < vsize && vids[0] < vsize && vids[3] < vsize)
			{
				MakeFace(all_v[vids[2]], all_v[vids[0]], all_v[vids[3]], NULL);
			}
		}
	}

		//not used
	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;

	//free mem
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
}


void	Crust( void )
{
	tVertex  ptr_v;
	tVertex * all_v = NULL;
	int vsize = 0, fsize = 0, esize = 0, tsize = 0;
	int id = 0;

	//global varibles for qhull
	static char * options = (char *)"qhull v Qz Qbb QJ Pp";
	int curlong, totlong;
	coordT * pt = NULL;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	facetT *neighbor, **neighborp;
	int vid = 0, i;
	tsFace sface;
	double volume;
	tTetra tetra;
	pointT s_radius;

	//variables for poles
	bool *is_on_convexhull;
	coordT *voronoi_vertex, *normals, *all_poles;
	tVertex vert;
	double dist, max_dist, inner, min_inner;
	int psize;
	double threshold = 1000;
	coordT pole[3], antipole[3];


	//count number of points
	ptr_v = vertices;
	do {
		vsize++;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//allocate memory
	pt = (coordT*)calloc(vsize * 4, sizeof(coordT)); //each point will have three coord
	all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
	assert(pt && all_v);

	is_on_convexhull = (bool*)calloc(vsize, sizeof(bool));
	normals = (coordT*)calloc(vsize * 3, sizeof(coordT));
	assert(is_on_convexhull && normals);

	all_poles = (coordT*)calloc(vsize * 6, sizeof(coordT));
	assert(all_poles);

	//
	// checking vertices are on the convexhulll and computing their normal normals
	//
	compute_normals(is_on_convexhull, normals);

	//
	// compute voronoi vertics using qhull
	//
	ptr_v = vertices;
	do {
		pt[id++] = ptr_v->v[0];
		pt[id++] = ptr_v->v[1];
		pt[id++] = ptr_v->v[2];
		pt[id++] = ptr_v->v[0] * ptr_v->v[0] + ptr_v->v[1] * ptr_v->v[1] + ptr_v->v[2] * ptr_v->v[2];
		all_v[ptr_v->vnum] = ptr_v;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	qh_init_A(stdin, stdout, stderr, 0, NULL);

	qh_initflags(options);
	qh_init_B(pt, vsize, 4, false);
	qh_qhull();
	qh_setvoronoi_all();


    //
	//compute poles and antipoles
	//
	
	psize = 0;
	FORALLvertices{
		vid = qh_pointid(vertex->point);

		//compute poles
		if (!is_on_convexhull[vid])
		{
			
			max_dist = -1.0;
			FOREACHneighbor_(vertex)
			{
				if (!neighbor->upperdelaunay)
				{
					voronoi_vertex = neighbor->center;
					dist = sqrt(pow(all_v[vid]->v[0] - voronoi_vertex[0], 2) + pow(all_v[vid]->v[1] - voronoi_vertex[1], 2) + pow(all_v[vid]->v[2] - voronoi_vertex[2], 2));
					if (dist > max_dist)
					{
						max_dist = dist;
						pole[0] = voronoi_vertex[0];
						pole[1] = voronoi_vertex[1];
						pole[2] = voronoi_vertex[2];
					}
				}

			}
			normals[vid * 3] = (pole[0] - all_v[vid]->v[0]) / max_dist;
			normals[vid * 3 + 1] = (pole[1] - all_v[vid]->v[1]) / max_dist;
			normals[vid * 3 + 2] = (pole[2] - all_v[vid]->v[2]) / max_dist;
			if (max_dist < threshold) //remove too far poles
			{
				all_poles[3 * psize] = pole[0]; all_poles[3 * psize + 1] = pole[1]; all_poles[3 * psize + 2] = pole[2];
				psize++;
			}
			
		}
		
	
			//compute antipoles

			min_inner = 100.0; //compute the negative projection
			FOREACHneighbor_(vertex)
			{
				if (!neighbor->upperdelaunay)
				{
					voronoi_vertex = neighbor->center;
					inner = (voronoi_vertex[0] - all_v[vid]->v[0])*normals[vid * 3] + (voronoi_vertex[1] - all_v[vid]->v[1])*normals[vid * 3 + 1] + (voronoi_vertex[2] - all_v[vid]->v[2])*normals[vid * 3 + 2];
					if (inner < min_inner)
					{
						min_inner = inner;
						antipole[0] = voronoi_vertex[0];
						antipole[1] = voronoi_vertex[1];
						antipole[2] = voronoi_vertex[2];
					}
				}
			}

			dist = sqrt(pow(antipole[0] - all_v[vid]->v[0], 2) + pow(antipole[1] - all_v[vid]->v[1], 2) + pow(antipole[2] - all_v[vid]->v[2], 2));
			if (dist < threshold) //remove too far poles
			{
				all_poles[3 * psize] = antipole[0]; all_poles[3 * psize + 1] = antipole[1]; all_poles[3 * psize + 2] = antipole[2];
				psize++;
			}

	}
	//not used
	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;

	

	//free mem
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);

	remove_nonsamples(psize, all_poles);

	free(normals);
	free(is_on_convexhull);
	free(all_poles);
	normals = NULL;
	is_on_convexhull = NULL;
	all_poles = NULL;
	
}

