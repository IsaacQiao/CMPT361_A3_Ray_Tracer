#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include "global.h"
#include "sphere.h"

using namespace std;
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
extern int reflect_on;
extern int step_max;
extern int chess_on;

/////////////////////////////////////////////////////////////////////

// this function check if the light ray intersect with shadows
bool check_shadow(Point o, Vector u, Spheres *sph)
{
  normalize(&u);

  while (sph) {
  float A = pow(u.x , 2) + pow(u.y , 2) + pow(u.z , 2);
  float B = 2 * (u.x * (o.x - sph->center.x) + u.y * (o.y - sph->center.y) + u.z * (o.z - sph->center.z));
  float C = pow(o.x - sph->center.x , 2) + pow(o.y - sph->center.y , 2) + pow(o.z - sph->center.z , 2) - pow(sph->radius , 2);

  float sqr = pow(B , 2) - 4 * A * C;

  float t1 = (-B + sqrt(sqr)) / (2*A);
  float t2 = (-B - sqrt(sqr)) / (2*A);

    if (sqr > 0 && t1 > 0.01 && t2 > 0.01) {
      return true;
    }

    sph = sph->next;
  }
  return false;
}

/*********************************************************************
 * Phong illumination - you need to implement this!
 *********************************************************************/
RGB_float phong(Point Intersect_point, Vector View_v, Vector surf_norm, Spheres *sph) {

  RGB_float C = {0, 0, 0}; //initialize

  // Phong's local illumination model
  // illumination = Global_ambient + ambient + (decay * diffiuse) + (decay * specular)
  // C = (1/(a+bd+cd^2)) * (Id * Kd *(n*l)) + (1/(a+bd+cd^2)) * (Is * Ks *(r*v)^N)
  // A = Iga * Kga + Ia * Ka
  // I = C + A

  // turn p and light1 into vector l
  Vector l = get_vec(Intersect_point, light1);
  normalize(&l);

  // d is the distance between the light source and the point on the object
  float d = vec_len(l);

  // 1/(a+bd+cd^2)
  float abcd = 1 / (decay_a + decay_b * d + decay_c * pow(d,2));

  // Diffuse with attenuation
  // (1/(a+bd+cd^2)) * (Id * Kd *(n*l))
  float nl = vec_dot(surf_norm, l); // (n*l)
  C.r += abcd * (light1_diffuse[0] * sph->mat_diffuse[0] * nl);
  C.g += abcd * (light1_diffuse[1] * sph->mat_diffuse[1] * nl);
  C.b += abcd * (light1_diffuse[2] * sph->mat_diffuse[2] * nl);

  // / get refected ray r
  float angle = vec_dot(surf_norm, l);
  if (angle < 0) angle = 0;
  Vector r = vec_minus(vec_scale(surf_norm, 2*angle), l);
  normalize(&r);

  // N is the shininess parameter for the object
  int N = sph->mat_shineness;

  // compute (r*v)^N
  float rv = vec_dot(r, View_v);
  float rvn = pow(rv, N);

  // Specular with attenuation
  // (1/(a+bd+cd^2)) * (Is * Ks *(r*v)^N)
  C.r += abcd * (light1_specular[0] * sph->mat_specular[0] * rvn);
  C.g += abcd * (light1_specular[1] * sph->mat_specular[1] * rvn);
  C.b += abcd * (light1_specular[2] * sph->mat_specular[2] * rvn);

  RGB_float A = {0, 0, 0};

  // Global ambient: Iga * Kga
  // Apply thevalues contained in the global ambient array to the sphere
  A.r += global_ambient[0] * sph->reflectance;
  A.g += global_ambient[1] * sph->reflectance;
  A.b += global_ambient[2] * sph->reflectance;

  // Ambient: Ia * Ka
  A.r += light1_ambient[0] * sph->mat_ambient[0];
  A.g += light1_ambient[1] * sph->mat_ambient[1];
  A.b += light1_ambient[2] * sph->mat_ambient[2];

  C.r += A.r;
  C.g += A.g;
  C.b += A.b;

  // Check if shadows are enabled
  // if so, change I to A
  if (shadow_on && check_shadow(Intersect_point, l, scene)){
    C = A;
  }
  // otherwise I is C

  return C;
}

/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(Point eye, Vector ray, int step_now) {
//
// do your thing here
//
	RGB_float color = background_clr;

  Point *Intersect_point = new Point;
  Spheres *sph = intersect_scene(eye, ray, scene, Intersect_point);
  
  Point board_hit;
  if (chess_on == 1 && intersect_chessboard(eye, ray, &board_hit))
  {
    Vector eye_vec = get_vec(board_hit, eye_pos);
    normalize(&eye_vec);

    color = board_color(board_hit);

    Vector shadow_v = get_vec(board_hit, light1);

    if (shadow_on && in_board(board_hit) && check_shadow(board_hit, shadow_v, scene)) {
      color = clr_scale(color, 0);
    }

    if (reflect_on == 1 && in_board(board_hit) && step_now < step_max) {
      Vector board_norm = {0, -1, 0};
      normalize(&board_norm);

      normalize(&ray);
      Vector reflected_ray = vec_plus(vec_scale(board_norm, -2 * vec_dot(board_norm, ray)), ray);

      RGB_float reflected_color = recursive_ray_trace(board_hit, reflected_ray, step_now + 1);

      color = clr_add(color, clr_scale(reflected_color, -0.35));
    }
  }

  if (sph != NULL){ // sph is the closest sphere intersec
    Vector View_v = get_vec(*Intersect_point, eye);
    normalize(&View_v);

    Vector surf_norm = sphere_normal(*Intersect_point, sph);
    normalize(&surf_norm);

    color = phong(*Intersect_point, View_v, surf_norm, sph);

    if (reflect_on == 1 && step_now < step_max){
      // get refected ray r
      float angle = vec_dot(surf_norm, vec_scale(ray, -1));
      if (angle < 0) angle = 0;
      Vector r = vec_plus(vec_scale(surf_norm, 2*angle), ray);
      normalize(&r);

      RGB_float reflected_color = clr_scale(recursive_ray_trace(*Intersect_point, r, step_now + 1), sph->reflectance);

      color = clr_add(color, reflected_color);
    }
  }

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
      
      ret_color = recursive_ray_trace(eye_pos, ray, 0);

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
