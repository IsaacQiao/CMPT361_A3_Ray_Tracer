#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"

//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];  

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;

// light 1 position and color
extern Point light1;
extern float light1_ambient[3];
extern float light1_diffuse[3];
extern float light1_specular[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int shadow_on;
extern int step_max;

/////////////////////////////////////////////////////////////////////

/*********************************************************************
 * Phong illumination - you need to implement this!
 *********************************************************************/
RGB_float phong(Point q, Vector v, Vector surf_norm, Spheres *sph) {


  RGB_float C = {0, 0, 0}; //initialize

  // Phong's local illumination model
  // illumination = Global_ambient + ambient + (decay * diffiuse) + (decay * specular)
  // C = (1/(a+bd+cd^2)) * (Id * Kd *(n*l)) + (1/(a+bd+cd^2)) * (Is * Ks *(r*v)^N)
  // A = Iga * Kga + Ia * Ka
  // I = C + A

  // turn q and light1 into vector l
  Vector l = get_vec(q, light1);
  normalize(&l);

  // d is the distance between the light source and the point on the object
  float d = vec_len(l);

  // 1/(a+bd+cd^2)
  float abcd = 1 / (decay_a + decay_b * d + decay_c * pow(d,2));

  // decay * Diffuse
  // (1/(a+bd+cd^2)) * (Id * Kd *(n*l))
  float nl = vec_dot(surf_norm, l); // (n*l)
  C.r += abcd * (light1_diffuse[0] * sph->mat_diffuse[0] * nl);
  C.g += abcd * (light1_diffuse[1] * sph->mat_diffuse[1] * nl);
  C.b += abcd * (light1_diffuse[2] * sph->mat_diffuse[2] * nl);

  // compute r first
  float angle = vec_dot(surf_norm, l);
  if (angle < 0) angle = 0;
  Vector scaled_surf_norm = vec_scale(surf_norm, 2*angle);
  Vector r = vec_minus(scaled_surf_norm, l);
  normalize(&r);

  // N is the shininess parameter for the object
  int N = sph->mat_shineness;

  // compute (r*v)^N
  float rv = vec_dot(r, v);
  float rvn = pow(rv, N);

  // decay*Specular
  // (1/(a+bd+cd^2)) * (Is * Ks *(r*v)^N)
  C.r += abcd * (light1_specular[0] * sph->mat_specular[0] * rvn);
  C.g += abcd * (light1_specular[1] * sph->mat_specular[1] * rvn);
  C.b += abcd * (light1_specular[2] * sph->mat_specular[2] * rvn);


  /*RGB_float ambient = {0, 0, 0};

  // Global ambient: Iga * Kga
  // Apply thevalues contained in the global ambient array to the sphere
  ambient.r += global_ambient[0] * sph->reflectance;
  ambient.g += global_ambient[1] * sph->reflectance;
  ambient.b += global_ambient[2] * sph->reflectance;


  // Ambient: Ia * Ka
  ambient.r += light1_ambient[0] * sph->mat_ambient[0];
  ambient.g += light1_ambient[1] * sph->mat_ambient[1];
  ambient.b += light1_ambient[2] * sph->mat_ambient[2];

  // Check if shadows are enabled
  bool intersectShadow = intersect_sphere_shadow(q, l, scene);
  if (shadow_on && intersectShadow){
    color = ambient;
  }*/


  return C;
}

/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace() {
//
// do your thing here
//
	RGB_float color;
	return color;
}

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color;
  Point cur_pixel_pos;
  Vector ray;

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;

  for (i=0; i<win_height; i++) {
    for (j=0; j<win_width; j++) {
      ray = get_vec(eye_pos, cur_pixel_pos);

      //
      // You need to change this!!!
      //
      // ret_color = recursive_ray_trace();
      // ret_color = background_clr; // just background for now
      Point *p = new Point;
      Spheres *sph = intersect_scene(eye_pos, ray, scene, p);

      if (sph != NULL){ // sph is the closest sphere intersec
        Vector v = get_vec(*p, eye_pos);
        normalize(&v);
        Vector surf_norm = sphere_normal(*p, sph);
        ret_color = phong(*p, v, surf_norm, sph);
      }

      else{//no intersection
        ret_color = background_clr;
      }
      
      // Parallel rays can be cast instead using below
      //
      // ray.x = ray.y = 0;
      // ray.z = -1.0;
      // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

// Checkboard for testing
//RGB_float clr = {float(i/32), 0, float(j/32)};
//ret_color = clr;

      frame[i][j][0] = GLfloat(ret_color.r);
      frame[i][j][1] = GLfloat(ret_color.g);
      frame[i][j][2] = GLfloat(ret_color.b);

      cur_pixel_pos.x += x_grid_size;
    }

    cur_pixel_pos.y += y_grid_size;
    cur_pixel_pos.x = x_start;
  }
}
