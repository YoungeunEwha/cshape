#include "alpha-shape.h"
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
void AlphaShape( unsigned int alpha )
{
	//implement 3d Alpha Shape alogorithm here
	//
	// put results (vertices, edges, faces, and tetrahedra)
	// in the global structures above
	// see main.c for detail on how qhull does it
	//

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
	int vid = 0, i;
	tTetra  tetra;
	tsFace face;
	double volume;
	pointT *center, s_radius; //center and squared radius of the tetrahedron


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

	//
	// compute convex hull in 4D by calling qhull
	// use flags: static char * options=(char *)"delaunay QJ Pp";
	//

	qh_init_A(stdin, stdout, stderr, 0, NULL);

	qh_initflags(options);
	qh_init_B(pt, vsize, 4, false);
	qh_qhull();
	qh_check_output();

	//
	//loop through all faces and call MakeNullTetra() to make a tetra 
	//remember that this is in 4D, so each face is a tetrahedron
	//
	//to fill the teta: get vertices of facet and loop through each vertex
	//use FOREACHvertex_()
	//

	//loop through all faces
	FORALLfacets
	{
		// compute the volume
		vid = 0;
		FOREACHvertex_(facet->vertices)
		{
			//get the vertex of the facet to compute volume
			face.vertex[vid++] = all_v[qh_pointid(vertex->point)];
		}

		volume = Volumei(&face, face.vertex[3]);

		//get the delaunay triangulation
		if (facet->normal[3] < -1.0e-5 && abs(volume) > 1.0e-5) //normal direction for the new axis && volume
		{
			tetra = tetra = MakeNullTetra();

			//get vertices of tetra
			//loop through each vertex
			vid = 0;
			FOREACHvertex_(facet->vertices)
			{
				//get the id of the vertex
				tetra->vertex[vid++] = all_v[qh_pointid(vertex->point)];
			}

			//get the cener of the circumsphere
			center = qh_facetcenter(facet->vertices);

			s_radius = 0.0;
			for (i = 0; i < 3; i++)
			{
				s_radius += (center[i] - tetra->vertex[0]->v[i]) * (center[i] - tetra->vertex[0]->v[i]);
			}
		
			//the circumsphere of the tetra has squared radius larger than alpha
			if (s_radius <= alpha)
			{
				MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[2], NULL);
				MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[3], NULL);
				MakeFace(tetra->vertex[1], tetra->vertex[2], tetra->vertex[3], NULL);
				MakeFace(tetra->vertex[2], tetra->vertex[0], tetra->vertex[3], NULL);
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