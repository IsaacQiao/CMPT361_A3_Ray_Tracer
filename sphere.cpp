#include "sphere.h"
#include <stdlib.h>
#include <math.h>

/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection, 
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/
float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit) {
	float A = pow(u.x , 2) + pow(u.y , 2) + pow(u.z , 2);
	float B = 2 * (u.x * (o.x - sph->center.x) + u.y * (o.y - sph->center.y) + u.z * (o.z - sph->center.z));
	float C = pow(o.x - sph->center.x , 2) + pow(o.y - sph->center.y , 2) + pow(o.z - sph->center.z , 2) - pow(sph->radius , 2);

	float sqr = pow(B , 2) - 4 * A * C;

	if (sqr < 0){
		return -1; //sqrt less than 0, no intersec
	}

	else{
		float t1 = (-B + sqrt(sqr)) / (2*A);
		float t2 = (-B - sqrt(sqr)) / (2*A);

		if (t2 < 0.001){
			return -1; // self shadow, no intersec
		}

		if (t1 >= 0.001){
			// set hit
			hit->x = o.x + t1 * u.x;
	    	hit->y = o.y + t1 * u.y;
	    	hit->z = o.z + t1 * u.z;
    	}

    	if (t2 >= 0.001){
			// set hit
			hit->x = o.x + t2 * u.x;
	    	hit->y = o.y + t2 * u.y;
	    	hit->z = o.z + t2 * u.z;
    	}

    	// calculate v and its length to be used in intersect_scene func
    	Vector v = {hit->x - o.x, hit->y - o.y, hit->z - o.z};
        
        return vec_len(v);
	}
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres *intersect_scene(Point o, Vector u, Spheres *sphs, Point *hit) {
	//return NULL;
    Spheres *closest = NULL;

    float min_dis = 10000000000000000;

    while (sphs != NULL) {// find the closest sph if there is sph intersec with 0->u
        float dis = intersect_sphere(o, u, sphs, hit);

        if ((dis != -1.0) && (dis < min_dis)) {
            min_dis = min_dis;
            closest = sphs;
        }

        sphs = sphs->next;
    }

    return closest;
}

/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}
