#include<bits/stdc++.h>
#include <GL/glut.h>  // GLUT, include glu.h and gl.h

#include "1805020_bitmap_image.hpp"

#include <cmath>
using namespace std;
#define pi acos(-1)
#define EPS 1e-9
class point {
public: 
        double x, y, z;
        double n;
    
    point() {
        x = 0;
        y = 0;
        z = 0;
        n = 1;
    }

    point(double _x, double _y, double _z) {
        x = _x;
        y = _y;
        z = _z;
        n = 1;
    }

    point(double _x, double _y, double _z, double _n) {
        x = _x;
        y = _y;
        z = _z;
        n = _n;
    }

    point(const point& p) {
        x = p.x;
        y = p.y;
        z = p.z;
        n = p.n;
    }

    point operator+(const point& p) {
        return point(x + p.x, y + p.y, z + p.z);
    }

    point operator-(const point& p) {
        return point(x - p.x, y - p.y, z - p.z);
    }

    point operator*(const double& d) {
        return point(x * d, y * d, z * d);
    }

    point operator/(const double& d) {
        return point(x / d, y / d, z / d);
    }

    double operator*(const point& p) {
        return x * p.x + y * p.y + z * p.z;
    }

    point operator-() {
        return point(-x, -y, -z);
    }

    point operator^(const point& p) {
        return point(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
    }
    
    double find_len() {
        double len = sqrt(x * x + y * y + z * z);
        return len;
    }

    void normalize() {
        double len = sqrt(x * x + y * y + z * z);
        x /= len;
        y /= len;
        z /= len;
    }
    void clip() {
        x = min(x, 1.0);
        y = min(y, 1.0);
        z = min(z, 1.0);

        x = max(x, 0.0);
        y = max(y, 0.0);
        z = max(z, 0.0);
    }
};
/* Initialize OpenGL Graphics */
void initGL() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   // Black and opaque
    glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
}


// ----------------------------defined variables-----------------------------//
// Global variables
bitmap_image tile1("texture_b.bmp");
bitmap_image tile2("texture_w.bmp");
point camera_target(0, 0, 0);
point pos(0, 160, 60);   // position of the eye
point l;     // look/forward direction
point r;     // right direction
point u(0, 0, 1);     // up direction


double far_plane, near_plane;

double field_of_vision;
double aspect_ratio;

double level_of_recursion;
double number_of_pixels;

double cell_width_board;

double ambient_checker, diffuse_checker, reflection_checker;

int number_of_objects;

double window_width;
double window_height;

double image_width; // number_of_pixels
double image_height; // number_of_pixels

int texture_mode;

int number_of_lights;

int number_of_spotlight;
bitmap_image image;

//----------------------------defined variables-----------------------------//

/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
        glColor3f(1,0,0);   // Red
        // X axis
        glVertex3f(0,0,0);
        glVertex3f(1,0,0);

        glColor3f(0,1,0);   // Green
        // Y axis
        glVertex3f(0,0,0);
        glVertex3f(0,1,0);

        glColor3f(0,0,1);   // Blue
        // Z axis
        glVertex3f(0,0,0);
        glVertex3f(0,0,1);
    glEnd();
}

void drawSphere(double radius, int stacks, int slices, point color) {
    struct point points[stacks+1][slices+1];
    for (int j = 0; j <= stacks; j++) {
        double phi = -M_PI / 2.0 + j * M_PI / stacks;
        double r = radius * cos(phi);
        double h = radius * sin(phi);
        for (int i = 0; i < slices+1; i++) {
            double theta = i * 2.0 * M_PI / slices;
            points[j][i].x = r * cos(theta);
            points[j][i].y = r * sin(theta);
            points[j][i].z = h;
        }
    }

    glBegin(GL_QUADS);
        for (int j = 0; j < stacks; j++) {
            for (int i = 0; i < slices; i++) {
                GLfloat c = (2+cos((i+j) * 2.0 * M_PI / slices)) / 3;
                glColor3f(color.x, color.y, color.z);
                glVertex3f(points[j][i].x, points[j][i].y, points[j][i].z);
                glVertex3f(points[j][i+1].x, points[j][i+1].y, points[j][i+1].z);

                glVertex3f(points[j+1][i+1].x, points[j+1][i+1].y, points[j+1][i+1].z);
                glVertex3f(points[j+1][i].x, points[j+1][i].y, points[j+1][i].z);
            }
        }
    glEnd();
}

class Ray {
    public:
    point start, dir;
    Ray() {
        start = point(0, 0, 0);
        dir = point(0, 0, 0);
    }
    Ray(point _start, point _dir) {
        start = _start;
        _dir.normalize();
        dir = _dir;
    }

};


class Light {
    public:
        point p; 
        point color;
        double fall_off_param;
        void draw() {
            // glPointSize(5);
            // glBegin(GL_POINTS);
            // glColor3f(0.5, 0.5, 0.5);
            // glVertex3f(p.x, p.y, p.z);
            // glEnd();
            glPushMatrix();
            {
                glColor3f(0.5, 0.5, 0.5);
                glTranslatef(p.x, p.y, p.z);
                //drawSphere(1, 20, 20, color);
                glutSolidSphere(6, 20, 20);
            }
            glPopMatrix();
        }
};

class SpotLight{
public:
    point p;
    point dir;
    point looking_point;
    double fall_off_param;
    double cutoffAngle; // this is different from the spotlight

    // void draw() {
    //     point color = pointLight.color;
    //     point p = pointLight.p;

    //     glPointSize(15);
    //     glBegin(GL_POINTS);
    //     glColor3f(color.x, color.y, color.z);
    //     glVertex3f(p.x, p.y, p.z);
    //     glEnd();
    // }
    void draw() {
        // glPointSize(5);
        // glBegin(GL_POINTS);
        // glColor3f(0.5, 0.5, 0.5);
        // glVertex3f(p.x, p.y, p.z);
        // glEnd();
        glPushMatrix();
        {
            glColor3f(0.5, 0.5, 0.5);
            glTranslatef(p.x, p.y, p.z);
            //drawSphere(1, 20, 20, color);
            glutSolidSphere(6, 20, 20);
        }
        glPopMatrix();
    }

};

void set_camera() {
    point up(0, 0, 1);
    point dir = camera_target - pos;
    dir.normalize();
    l = dir;
    r = up ^ dir;
    r.normalize();
    u = dir ^ r;
    u.normalize();

}

class TriangularIntersection{
    public:
    point v1, v2, v3;
    point normal;
    TriangularIntersection() {
        v1 = point(0, 0, 0);
        v2 = point(0, 0, 0);
        v3 = point(0, 0, 0);
        normal = point(0, 0, 0);
    }
    TriangularIntersection(point _v1, point _v2, point _v3) {
        v1 = _v1;
        v2 = _v2;
        v3 = _v3;
        normal = (v2 - v1) ^ (v3 - v1);
        normal.normalize();
    }

    double getArea(point a, point b) {
        point v = a ^ b;

        return fabs(v.find_len());
    }

    bool checkInTriangle(point pt) {
        double area = getArea(v2 - v1, pt - v1) + getArea(v3 - v1, pt - v1) + getArea(v3 - v2, pt - v2);
        double area2 = getArea(v2 - v1, v3 - v1);
        if(abs(area - area2) < 1e-9) return true;
        return false;
    }
    double getIntersection(Ray ray) {
        double denom = ray.dir * normal;
        //denom = fabs(denom);
        if(abs(denom) < 1e-9) return -1;
        double t = (normal * (v1 - ray.start)) / denom;
        if(t < 0) return -1;
        //t = fabs(t);
        if(!checkInTriangle(ray.start + ray.dir * t)) return -1; 
        return t;
    }
};


class Floor {
    public: 
        point ref_floor;
        double width;
        double height;
        Floor() {
            ref_floor = point(-50, -50, 0);
            width = 1;
            height = 1;
        }

        Floor(double _width) {
            height = _width;
            width = _width;
            ref_floor = point(-(width * 100) / 2, -(width * 100) / 2, 0);
        }
        void drawFloor() {
            for(int i = 0; i < 100; i++) {
                for(int j = 0; j < 100; j++) {
                    if((i + j) % 2 == 0) {
                        glColor3f(0, 0, 0);
                    } else {
                        glColor3f(1, 1, 1);
                    }
                    glBegin(GL_QUADS);
                    {
                        glVertex3f(ref_floor.x + i * height, ref_floor.y + j * width, 0);
                        glVertex3f(ref_floor.x + (i + 1) * height, ref_floor.y + j * width, 0);
                        glVertex3f(ref_floor.x + (i + 1) * height, ref_floor.y + (j + 1) * width, 0);
                        glVertex3f(ref_floor.x + i * height, ref_floor.y + (j + 1) * width, 0);
                    
                    }
                    glEnd();
                }
            }
        }

        point getcolor(point pt) {
            int x = (pt.x - ref_floor.x) / width;
            int y = (pt.y - ref_floor.y) / height;
            if(x < 0 or x >= 100 * 100 or y < 0 or y >= 100 * 100) return point(0, 0, 0);
            
            if(!texture_mode) {
                if((x + y) % 2 == 0) {
                    return point(0, 0, 0);
                } else {
                    return point(1, 1, 1);
                }
            }
            double i = (pt.x - ref_floor.x) - x * width;
            double j = (pt.y - ref_floor.y) - y * width;
            if((x + y) % 2 == 0) {
                int px = (int)((i / width) * tile1.width()) % tile1.width();
                int py = (int)((j / height)* tile1.height()) % tile1.height();
                unsigned char r, g, b;
                tile1.get_pixel(px, py, r, g, b);
                return point(r/ 255.0, g / 255.0, b / 255.0);
            } else {
                int px = (int)((i / width) * tile2.width()) % tile2.width();
                int py = (int)((j / height)* tile2.height()) % tile2.height();
                unsigned char r, g, b;
                tile2.get_pixel(px, py, r, g, b);
                return point(r/ 255.0, g / 255.0, b / 255.0);
            }
        }
        double intersect(Ray ray, point& color_c, int level) {
            double t = -1;
            double a = ray.dir * point(0, 0, 1);
            if(abs(a) < 1e-9) return -1;
            t = (point(0, 0, 1) * (ref_floor - ray.start)) / a;

            point p = ray.start + ray.dir * t;
            if(p.x <= ref_floor.x) return -1;
            if(p.x >= abs(ref_floor.x) and p.y <= ref_floor.y and p.y >= abs(ref_floor.y)) return -1;

            if(level == 1) color_c = getcolor(p);
            return t;
        }
};



class sphere {
    public:
        point center;
        double radius;
        int slices, stacks;
        point color;
        double ambient, diffuse, specular, reflection;
        double shineness;
        sphere() {
            center = point(0, 0, 0.5);
            radius = 1;
            slices = 32;
            stacks = 32;
            color = point(1, 1, 1);
        }
        sphere(point _center, double _radius, int _slices, int _stacks) {
            center = _center;
            radius = _radius;
            slices = _slices;
            stacks = _stacks;
        }
        void drawSph() {
            glPushMatrix();
            {
                glTranslatef(center.x, center.y, center.z);
                drawSphere(radius, stacks, slices, color);
            }
            glPopMatrix();
        }

        double sphere_intersect(Ray ray, point& norm) {
            point x = ray.start - center;
            double a = ray.dir * ray.dir;
            double b = 2 * (ray.dir * x);
            double c = (x * x) - radius * radius;
            double nishchayok = b * b - 4 * a * c;
            double t = -1;
            if(nishchayok < EPS) return -1;
            else {
                if(abs(a) < EPS) {
                    t = -c / b;
                    return t;
                }
                double t1 = (-b + sqrt(nishchayok)) / (2 * a);
                double t2 = (-b - sqrt(nishchayok)) / (2 * a);
                if(t2 < t1) swap(t1, t2);
                if(t1 > EPS) t = t1;
                else if(t2 > EPS) t = t2;

                point p = ray.start + ray.dir * t;
                norm =  p - center;
                norm.normalize();
                return t;
            }
        } 
};

void drawCubeBySachin(point ref, double length, point color) {
    glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
        // Top face (y = 1.0f)
        // Define vertices in counter-clockwise (CCW) order with normal pointing out
        //glColor3f(0.0f, 1.0f, 0.0f);     // Green
        glColor3f(color.x, color.y, color.z);
        glVertex3f(ref.x, ref.y, ref.z);
        glVertex3f(ref.x + length, ref.y, ref.z);
        glVertex3f(ref.x + length, ref.y + length, ref.z);
        glVertex3f(ref.x, ref.y + length, ref.z);

        //glColor3f(1.0f, 0.5f, 0.0f);     // Orange
        glColor3f(color.x, color.y, color.z);
        glVertex3f(ref.x, ref.y, ref.z + length);
        glVertex3f(ref.x + length, ref.y, ref.z + length);
        glVertex3f(ref.x + length, ref.y + length, ref.z + length);
        glVertex3f(ref.x, ref.y + length, ref.z + length);
        
        //glColor3f(1.0f, 0.0f, 0.0f);     // Red
        
        glColor3f(color.x, color.y, color.z);
        glVertex3f(ref.x, ref.y, ref.z);
        glVertex3f(ref.x + length, ref.y, ref.z);
        glVertex3f(ref.x + length, ref.y, ref.z + length);
        glVertex3f(ref.x, ref.y, ref.z + length);

        //glColor3f(1.0f, 1.0f, 0.0f);     // Yellow
        glColor3f(color.x, color.y, color.z);
        glVertex3f(ref.x, ref.y + length, ref.z);
        glVertex3f(ref.x + length, ref.y + length, ref.z);
        glVertex3f(ref.x + length, ref.y + length, ref.z + length);
        glVertex3f(ref.x, ref.y + length, ref.z + length);

        //glColor3f(0.0f, 0.0f, 1.0f);     // Blue
        glColor3f(color.x, color.y, color.z);
        glVertex3f(ref.x, ref.y, ref.z);
        glVertex3f(ref.x, ref.y + length, ref.z);
        glVertex3f(ref.x, ref.y + length, ref.z + length);
        glVertex3f(ref.x, ref.y, ref.z + length);

        //glColor3f(1.0f, 0.0f, 1.0f);     // Magenta
        glColor3f(color.x, color.y, color.z);
        glVertex3f(ref.x + length, ref.y, ref.z);
        glVertex3f(ref.x + length, ref.y + length, ref.z);
        glVertex3f(ref.x + length, ref.y + length, ref.z + length);
        glVertex3f(ref.x + length, ref.y, ref.z + length);
        
        
    glEnd();  // End of drawing color-cube
}


class cube {
    public:
    point center;
    double side;
    point color;
    double ambient, diffuse, specular, reflection;
    double shineness;

    vector<point> points; 
    cube() {
        center = point(0, 0, 0);
        side = 1;
        color = point(1, 1, 1);
    }
    cube(point _center, double _side) {
        center = _center;
        side = _side;


    }
    void draw_cube() {
        glPushMatrix();
        {
            glTranslatef(center.x, center.y, center.z);
            drawCubeBySachin(center, side, color);
        }
        glPopMatrix();
    }

    void set_points() {
        /// first face

        points.push_back(point(center.x, center.y, center.z));
        points.push_back(point(center.x + side, center.y, center.z));
        points.push_back(point(center.x + side, center.y + side, center.z));
        points.push_back(point(center.x, center.y + side, center.z));


        /// second face

        points.push_back(point(center.x, center.y, center.z + side));
        points.push_back(point(center.x + side, center.y, center.z + side));
        points.push_back(point(center.x + side, center.y + side, center.z + side));
        points.push_back(point(center.x, center.y + side, center.z + side));

        /// third face

        points.push_back(point(center.x, center.y, center.z));
        points.push_back(point(center.x + side, center.y, center.z));
        points.push_back(point(center.x + side, center.y, center.z + side));
        points.push_back(point(center.x, center.y, center.z + side));

        /// fourth face
        // glVertex3f(ref.x, ref.y + length, ref.z);
        // glVertex3f(ref.x + length, ref.y + length, ref.z);
        // glVertex3f(ref.x + length, ref.y + length, ref.z + length);
        // glVertex3f(ref.x, ref.y + length, ref.z + length);
        points.push_back(point(center.x, center.y + side, center.z));
        points.push_back(point(center.x + side, center.y + side, center.z));
        points.push_back(point(center.x + side, center.y + side, center.z + side));
        points.push_back(point(center.x, center.y + side, center.z + side));


        /// fifth face
        // glVertex3f(ref.x, ref.y, ref.z);
        // glVertex3f(ref.x, ref.y + length, ref.z);
        // glVertex3f(ref.x, ref.y + length, ref.z + length);
        // glVertex3f(ref.x, ref.y, ref.z + length);
        points.push_back(point(center.x, center.y, center.z));
        points.push_back(point(center.x, center.y + side, center.z));
        points.push_back(point(center.x, center.y + side, center.z + side));
        points.push_back(point(center.x, center.y, center.z + side));


        /// sixth face
        // glColor3f(color.x, color.y, color.z);
        // glVertex3f(ref.x + length, ref.y, ref.z);
        // glVertex3f(ref.x + length, ref.y + length, ref.z);
        // glVertex3f(ref.x + length, ref.y + length, ref.z + length);
        // glVertex3f(ref.x + length, ref.y, ref.z + length);
        points.push_back(point(center.x + side, center.y, center.z));
        points.push_back(point(center.x + side, center.y + side, center.z));
        points.push_back(point(center.x + side, center.y + side, center.z + side));
        points.push_back(point(center.x + side, center.y, center.z + side));
    }
    double cube_intersect(Ray ray, point& norm) {
        double t_min = 1000000;
        //cout << "printinig points size " << points.size() << "\n";  
        for(int i = 0; i < points.size(); i += 4) {
            point p1 = points[i];
            point p2 = points[i + 1];
            point p3 = points[i + 2];
            point p4 = points[i + 3];

            TriangularIntersection triangle(p1, p2, p3);
            double t = triangle.getIntersection(ray);
            //cout << t << " ";
            if(t > EPS and t < t_min) {
                t_min = min(t_min, t);
                point temp = ray.start + ray.dir * t;
                if(temp * triangle.normal < 0) {
                    norm = triangle.normal;
                } else {
                    norm = triangle.normal * -1;
                }
            }

            triangle = TriangularIntersection(p1, p3, p4);
            t = triangle.getIntersection(ray);
            //cout << t << "\n";
            if(t > EPS and t < t_min) {
                t_min = min(t_min, t);
                point temp = ray.start + ray.dir * t;
                if(temp * triangle.normal < 0) {
                    norm = triangle.normal;
                } else {
                    norm = triangle.normal * -1;
                }
            }
        }
        return t_min;
    }
};

class pyramid {

public:
    
    point bottom_point;
    double width, height;
    point color;
    double ambient, diffuse, specular, reflection;
    double shineness;
    vector<point> points; 


    pyramid() {
        bottom_point = point(0, 0, 0);
        width = 1;
        height = 1;
        color = point(1, 1, 1);
    }
    
    pyramid(point _bottom_point, double _width, double _height) {
        bottom_point = _bottom_point;
        width = _width;
        height = _height;
    }
    
    void draw_pyramid() {
        point p[5];
        p[0] = bottom_point;
        p[1] = bottom_point + point(width, 0, 0);
        p[2] = bottom_point + point(0, width, 0);
        p[3] = bottom_point + point(width, width, 0);
        p[4] = bottom_point + point(width / 2, width / 2, height);

        for(int i = 0; i < 5; i++) {
            points.push_back(p[i]);
        }

        glBegin(GL_QUADS);
        {
            glColor3f(color.x, color.y, color.z);
            //glColor3f(1, 0, 0);
            glVertex3f(p[0].x, p[0].y, p[0].z);
            glVertex3f(p[1].x, p[1].y, p[1].z);
            glVertex3f(p[3].x, p[3].y, p[3].z);
            glVertex3f(p[2].x, p[2].y, p[2].z);
        }
        glEnd();

        glBegin(GL_TRIANGLES);
        {
            glColor3f(color.x, color.y, color.z);
            //glColor3f(0, 1, 0);
            glVertex3f(p[0].x, p[0].y, p[0].z);
            glVertex3f(p[1].x, p[1].y, p[1].z);
            glVertex3f(p[4].x, p[4].y, p[4].z);
        }
        glEnd();

        glBegin(GL_TRIANGLES);
        {
            glColor3f(color.x, color.y, color.z);
            //glColor3f(0, 0, 1);
            glVertex3f(p[1].x, p[1].y, p[1].z);
            glVertex3f(p[3].x, p[3].y, p[3].z);
            glVertex3f(p[4].x, p[4].y, p[4].z);
        }
        glEnd();

        glBegin(GL_TRIANGLES);
        {
            glColor3f(color.x, color.y, color.z);
            //glColor3f(1, 1, 0);
            glVertex3f(p[3].x, p[3].y, p[3].z);
            glVertex3f(p[2].x, p[2].y, p[2].z);
            glVertex3f(p[4].x, p[4].y, p[4].z);
        }
        glEnd();

        glBegin(GL_TRIANGLES);
        {
            glColor3f(color.x, color.y, color.z);
            //glColor3f(0, 1, 1);
            glVertex3f(p[2].x, p[2].y, p[2].z);
            glVertex3f(p[0].x, p[0].y, p[0].z);
            glVertex3f(p[4].x, p[4].y, p[4].z);
        }
        glEnd();
    }

    double pyramid_intersect(Ray ray, point& norm) {
        double t_min = 1000000;
        TriangularIntersection triangle(points[0], points[1], points[3]);
        double t = triangle.getIntersection(ray);
        if(t > EPS and t < t_min) {
            t_min = t;
            point temp = ray.start + ray.dir * t;
            if(temp * triangle.normal < 0) {
                norm = triangle.normal;
            } else {
                norm = triangle.normal * (-1);
            }
        }
        triangle = TriangularIntersection(points[0], points[3], points[2]);
        t = triangle.getIntersection(ray);
        if(t > EPS and t < t_min) {
            t_min = t;
            point temp = ray.start + ray.dir * t;
            if(temp * triangle.normal < 0) {
                norm = triangle.normal;
            } else {
                norm = triangle.normal * (-1);
            }
        }
        triangle = TriangularIntersection(points[0], points[1], points[4]);
        t = triangle.getIntersection(ray);
        if(t > EPS and t < t_min) {
            t_min = t;
            point temp = ray.start + ray.dir * t;
            if(temp * triangle.normal < 0) {
                norm = triangle.normal;
            } else {
                norm = triangle.normal * (-1);
            }
        }
        triangle = TriangularIntersection(points[1], points[2], points[4]);
        t = triangle.getIntersection(ray);
        if(t > EPS and t < t_min) {
            t_min = t;
            point temp = ray.start + ray.dir * t;
            if(temp * triangle.normal < 0) {
                norm = triangle.normal;
            } else {
                norm = triangle.normal * (-1);
            }
        }
        triangle = TriangularIntersection(points[2], points[3], points[4]);
        t = triangle.getIntersection(ray);
        if(t > EPS and t < t_min) {
            t_min = t;
            point temp = ray.start + ray.dir * t;
            if(temp * triangle.normal < 0) {
                norm = triangle.normal;
            } else {
                norm = triangle.normal * (-1);
            }
        }
        triangle = TriangularIntersection(points[3], points[0], points[4]);
        t = triangle.getIntersection(ray);
        if(t > EPS and t < t_min) {
            t_min = t;
            point temp = ray.start + ray.dir * t;
            if(temp * triangle.normal < 0) {
                norm = triangle.normal;
            } else {
                norm = triangle.normal * (-1);
            }
        }
        return t_min;
    }
};


vector<sphere> spheres;
vector<cube> cubes;
vector<pyramid> pyramids;
vector<Light> lights;
vector<SpotLight> spotlights;
Floor ground;



point get_color(Ray ray, int level) {
    if(level == 0) return point(0, 0, 0);

    double t_min = 1000000;
    point curr_normal;

    int nearest_object_index = -1;
    point color_c(0, 0, 0);

    double t = ground.intersect(ray, color_c, 0);
    if(t > (double)1e-9) {
        t_min = t;
        nearest_object_index = 100;
        curr_normal = point(0, 0, 1);
        point pt = ray.start + ray.dir * t;
        if(pt * curr_normal < 0) {
            curr_normal = curr_normal;
        } else {
            curr_normal = curr_normal * (-1);
        }
    }

    for(int k = 0; k < cubes.size(); k++) {
        point norm;
        t = cubes[k].cube_intersect(ray, norm);
        //if(t > (double)(1e5)) continue;
        //if (t > EPS) cout << "WHY ? " << t << endl;
        //cout << t << "\n";
        if(t > (double)EPS and t < t_min) {
            t_min = t;
            nearest_object_index = k;
            curr_normal = norm;
        }
    }

    for(int k = 0; k < pyramids.size(); k++) {
        point norm;
        t = pyramids[k].pyramid_intersect(ray, norm);
        //cout << t << "\n";
        if(t > (double)EPS and t < t_min) {
            t_min = t;
            nearest_object_index = k + 1000;
            curr_normal = norm;
        }
    }

    for(int k = 0; k < spheres.size(); k++) {
        point norm;
        t = spheres[k].sphere_intersect(ray, norm);
        if(t > (double)EPS and t < t_min) {
            t_min = t;
            nearest_object_index = k + 10000;
            curr_normal = norm;
        }
    }
    
    if(nearest_object_index == -1) return point(0, 0, 0);

    point fina = ray.start + ray.dir * t_min;

    point ini = ray.start;

    point ray_v = fina - ini;

    curr_normal = curr_normal * ((ini - fina) * curr_normal);
    
    point temp = ray_v + curr_normal * 2;
    //temp = temp * 2;
    point next_dir = temp;
    curr_normal.normalize();
    next_dir.normalize();

    Ray new_ray;
    new_ray.start = fina;
    new_ray.dir = next_dir;
    
    long double lambert = 0;
    long double phong = 0;

    for(int i = 0; i < lights.size(); i++) {
        point ray_vec = lights[i].p - fina;
        double len = ray_vec.find_len();
        double distance = len;
        ray_vec.normalize();
        int flag = 0;
        for(int j = 0; j < cubes.size() and flag == 0; j++) {
            point norm;
            Ray ray2(fina, ray_vec);
            double t = cubes[j].cube_intersect(ray2, norm);
            if(t > (double)EPS and t < len) {
                flag = 1;
            }
        }
        for(int j = 0; j < pyramids.size() and flag == 0; j++) {
            point norm;
            Ray ray2(fina, ray_vec);
            double t = pyramids[j].pyramid_intersect(ray2, norm);
            if(t > (double)EPS and t < len) {
                flag = 1;
            }
        }
        for(int j = 0; j < spheres.size() and flag == 0; j++) {
            point norm;
            Ray ray2(fina, ray_vec);
            double t = spheres[j].sphere_intersect(ray2, norm);
            if(t > (double)EPS and t < len) {
                flag = 1;
            }
        }
        if(flag == 1) continue;
        long double shine = 0;
        if(nearest_object_index == 100) {
            shine = 0;
        }
        else if(nearest_object_index >= 0 and nearest_object_index < 100) {
            shine = cubes[nearest_object_index].shineness;
        }
        else if(nearest_object_index >= 1000 and nearest_object_index < 10000) {
            shine = pyramids[nearest_object_index - 1000].shineness;
        }
        else if(nearest_object_index >= 10000) {
            shine = spheres[nearest_object_index - 10000].shineness;
        }
        long double scaling_factor = exp(-1 * distance * distance * lights[i].fall_off_param);
        lambert += (ray_vec * curr_normal) * scaling_factor;
        phong += pow(next_dir * ray_vec, shine) * scaling_factor;

    }


    for(int i = 0; i < spotlights.size(); i++) {
        point ray_vec = spotlights[i].p - fina;
        double len = ray_vec.find_len();
        double distance = len;
        ray_vec.normalize();
        int flag = 0;
        point direction_of_source = spotlights[i].looking_point - spotlights[i].p;
        direction_of_source.normalize();
        double angle = acos(direction_of_source * ray_vec);
        angle *= (180 / pi);
        if(angle > spotlights[i].cutoffAngle) flag = 1;
        for(int j = 0; j < cubes.size() and flag == 0; j++) {
            point norm;
            Ray ray2(fina, ray_vec);
            double t = cubes[j].cube_intersect(ray2, norm);
            if(t > (double)EPS and t < len) {
                flag = 1;
            }
        }
        for(int j = 0; j < pyramids.size() and flag == 0; j++) {
            point norm;
            Ray ray2(fina, ray_vec);
            double t = pyramids[j].pyramid_intersect(ray2, norm);
            if(t > (double)EPS and t < len) {
                flag = 1;
            }
        }
        for(int j = 0; j < spheres.size() and flag == 0; j++) {
            point norm;
            Ray ray2(fina, ray_vec);
            double t = spheres[j].sphere_intersect(ray2, norm);
            if(t > (double)EPS and t < len) {
                flag = 1;
            }
        }
        if(flag == 1) continue;
        long double shine = 0;
        if(nearest_object_index == 100) {
            shine = 0;
        }
        else if(nearest_object_index >= 0 and nearest_object_index < 100) {
            shine = cubes[nearest_object_index].shineness;
        }
        else if(nearest_object_index >= 1000 and nearest_object_index < 10000) {
            shine = pyramids[nearest_object_index - 1000].shineness;
        }
        else if(nearest_object_index >= 10000) {
            shine = spheres[nearest_object_index - 10000].shineness;
        }
        long double scaling_factor = exp(-1 * distance * distance * lights[i].fall_off_param);
        lambert += (ray_vec * curr_normal) * scaling_factor;
        phong += pow(next_dir * ray_vec, shine) * scaling_factor;

    }
    

    if(nearest_object_index == 100) {
        color_c = {0, 0, 0};
        ground.intersect(ray, color_c, 1);
        //color_c.clip();
        color_c = color_c * ambient_checker + get_color(new_ray, level - 1) * reflection_checker + color_c * diffuse_checker * lambert;
        return color_c;
    }
    else if(nearest_object_index >= 0 and nearest_object_index < 100) {
        color_c = cubes[nearest_object_index].color;
        //color_c.clip();
        color_c = color_c * cubes[nearest_object_index].ambient + get_color(new_ray, level - 1) * cubes[nearest_object_index].reflection + color_c * cubes[nearest_object_index].diffuse * lambert + color_c * cubes[nearest_object_index].specular * phong;
        //color_c = color_c * lambert * cubes[nearest_object_index].diffuse + color_c * phong * cubes[nearest_object_index].specular;
        return color_c;
    }
    else if(nearest_object_index >= 1000 and nearest_object_index < 10000) {
        color_c = pyramids[nearest_object_index - 1000].color;
        //color_c.clip();
        color_c = color_c * pyramids[nearest_object_index - 1000].ambient + get_color(new_ray, level - 1) * pyramids[nearest_object_index - 1000].reflection + color_c * pyramids[nearest_object_index - 1000].diffuse * lambert + color_c * pyramids[nearest_object_index - 1000].specular * phong;
        return color_c;
    }
    else if(nearest_object_index >= 10000) {
        color_c = spheres[nearest_object_index - 10000].color;
        //color_c.clip();
        color_c = color_c * spheres[nearest_object_index - 10000].ambient + get_color(new_ray, level - 1) * spheres[nearest_object_index - 10000].reflection + color_c * spheres[nearest_object_index - 10000].diffuse * lambert + color_c * spheres[nearest_object_index - 10000].specular * phong;
        return color_c;
    }
    return point(0, 0, 0);
}

void capture() {
    cout << "capturing the image\n";

    for(int i = 0; i < image_width; i++) {
        for(int j = 0; j < image_height; j++) {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    double screen_height = 2 * near_plane * tan((field_of_vision * pi) / 360);
    double screen_width = screen_height / aspect_ratio;

    point screen_mid = pos + l * near_plane;
    double du = screen_width / (image_width);
    double dv = screen_height / (image_height);

    // point top_left = mid + (r * (screen_width / 2)) - (u * (screen_height / 2));
    point top_left = screen_mid + (u * (screen_height / 2)) - (r * (screen_width / 2));
    // cout << screen_mid.x << " " << screen_mid.y << " " << screen_mid.z << "\n";
    // cout << top_left.x << " " << top_left.y << " " << top_left.z << "\n";
    // cout << du << " " << dv << "\n";

    int nearest_object_index = -1;

    //point point_buffer[(int)image_width][(int)image_height];

    for(int i = 0; i < image_height; i++) {
        for(int j = 0; j < image_width; j++) {
            //point_buffer[i][j] = top_left + (r * (i * du)) - (u * (j * dv));
            point pixel = top_left + (r * (j * du)) - (u * (i * dv));
            Ray ray(pixel, pixel - pos);
            point color = get_color(ray, level_of_recursion);
            color.clip();
            image.set_pixel(j, i, color.x * 255, color.y * 255, color.z * 255);
        }
    }
    cout << "Saving the image\n";
    image.save_image("out.bmp");
    cout << "Reached here\n";
}


void take_input() {
    ifstream input_file("input.txt");
    if(input_file.is_open()) {
        input_file >> near_plane >> far_plane;
        input_file >> field_of_vision;
        input_file >> aspect_ratio;
        input_file >> level_of_recursion;
        input_file >> number_of_pixels;
        input_file >> cell_width_board;
        input_file >> ambient_checker >> diffuse_checker >> reflection_checker;
        input_file >> number_of_objects;

        cout << number_of_objects << "\n";
        for(int i = 0; i < number_of_objects; i++) {
            string type; input_file >> type;
            if(type == "cube") {
                cube c;
                input_file >> c.center.x >> c.center.y >> c.center.z;
                input_file >> c.side;
                input_file >> c.color.x >> c.color.y >> c.color.z;
                input_file >> c.ambient >> c.diffuse >> c.specular >> c.reflection;
                input_file >> c.shineness;

                cout << c.center.x << " " << c.center.y << " " << c.center.z << "\n";
                cout << c.side << "\n";
                cout << c.color.x << " " << c.color.y << " " << c.color.z << "\n";
                cout << c.ambient << " " << c.diffuse << " " << c.reflection << "\n";
                cout << c.shineness << "\n";

                c.set_points();
                cubes.push_back(c);
            }
            else if(type == "sphere") {
                sphere s;
                input_file >> s.center.x >> s.center.y >> s.center.z;
                input_file >> s.radius;
                input_file >> s.color.x >> s.color.y >> s.color.z;
                input_file >> s.ambient >> s.diffuse >> s.specular >> s.reflection;
                input_file >> s.shineness;

                cout << s.center.x << " " << s.center.y << " " << s.center.z << "\n";
                spheres.push_back(s);
            }
            else {
                pyramid p;
                input_file >> p.bottom_point.x >> p.bottom_point.y >> p.bottom_point.z;
                input_file >> p.width >> p.height;
                input_file >> p.color.x >> p.color.y >> p.color.z;
                input_file >> p.ambient >> p.diffuse >> p.specular >> p.reflection;
                input_file >> p.shineness;

                cout << p.bottom_point.x << " " << p.bottom_point.y << " " << p.bottom_point.z << "\n";
                pyramids.push_back(p);
            }
        }


        input_file >> number_of_lights;
        for(int i = 0; i < number_of_lights; i++) {
            Light l;
            input_file >> l.p.x >> l.p.y >> l.p.z;
            input_file >> l.fall_off_param;
            lights.push_back(l);
        }

        input_file >> number_of_spotlight;
        for(int i = 0; i < number_of_spotlight; i++) {
            SpotLight sp; 
            input_file >> sp.p.x >> sp.p.y >> sp.p.z;
            input_file >> sp.fall_off_param;
            input_file >> sp.looking_point.x >> sp.looking_point.y >> sp.looking_point.z;
            input_file >> sp.cutoffAngle;
            spotlights.push_back(sp);
        }
        cout << "Successfully taken all the inputs\n";
        for(int i = 0; i < cubes.size(); i++) {
            cout << "Cube " << i << ": ";
            cout << cubes[i].center.x << " " << cubes[i].center.y << " " << cubes[i].center.z << " ";
            cout << "\n";
        }
        for(int i = 0; i < spheres.size(); i++) {
            cout << "Sphere " << i << ": ";
            cout << spheres[i].center.x << " " << spheres[i].center.y << " " << spheres[i].center.z << " ";
            cout << "\n";
        }
        for(int i = 0; i < pyramids.size(); i++) {
            cout << "Pyramid " << i << ": ";
            cout << pyramids[i].bottom_point.x << " " << pyramids[i].bottom_point.y << " " << pyramids[i].bottom_point.z << " ";
            cout << "\n";
        }
        image_height = image_width = number_of_pixels;

    } else {
        cout << "Unable to open file" << endl;
    }
    //input_file.open("input.txt");
}


/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display() {
    
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);             // To operate on Model-View matrix
    glLoadIdentity();                       // Reset the model-view matrix

    // default arguments of gluLookAt
    // gluLookAt(0,0,0, 0,0,-100, 0,1,0);

    // control viewing (or camera)
    gluLookAt(pos.x,pos.y,pos.z,
              pos.x+l.x,pos.y+l.y,pos.z+l.z,
              u.x,u.y,u.z);
    
    drawAxes();
    // draw
    // drawSphere(5,100,100);
    // Floor floor;
    // floor.drawFloor();
    
    ground = Floor(cell_width_board);
    ground.drawFloor();

    // sphere sph(point(0, 0, 0.5), 0.5, 32, 32);
    // sph.drawSph();
    
    // cube c(point(3, 3, 0), 1);
    // c.draw_cube();

    // pyramid p(point(3, 0, 0), 1, 1);
    // p.draw_pyramid();
    
    for(int i = 0; i < spheres.size(); i++) {
        spheres[i].drawSph();
    }
    for(int i = 0; i < cubes. size(); i++) {
        cubes[i].draw_cube();
    }
    for(int i = 0; i < pyramids.size(); i++) {
        pyramids[i].draw_pyramid();
    }
    for(int i = 0; i < lights.size(); i++) {
        //cout << "light - i: " << i << " " << lights[i].p.x << " " << lights[i].p.y << " " << lights[i].p.z << "\n";
        lights[i].draw();
    }

    for(int i = 0; i < spotlights.size(); i++) {
        spotlights[i].draw();
    }
    glutSwapBuffers();  // Render now
}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    /*if (width >= height) {
        // aspect >= 1, set the height from -1 to 1, with larger width
        gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
    } else {
        // aspect < 1, set the width to -1 to 1, with larger height
        gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
    }*/
    // Enable perspective projection with fovy, aspect, zNear and zFar
    //gluPerspective(45.0f, aspect, 0.1f, 100.0f);
    gluPerspective(field_of_vision, aspect_ratio, near_plane, far_plane);
}




void keyboardListener(unsigned char key, int xx,int yy){
    double rate = 0.01;
	switch(key){

        case '0':
            capture();
            break;

		case '1':
			r.x = r.x*cos(rate)+l.x*sin(rate);
			r.y = r.y*cos(rate)+l.y*sin(rate);
			r.z = r.z*cos(rate)+l.z*sin(rate);

			l.x = l.x*cos(rate)-r.x*sin(rate);
			l.y = l.y*cos(rate)-r.y*sin(rate);
			l.z = l.z*cos(rate)-r.z*sin(rate);
            
			break;

        case '2':
			r.x = r.x*cos(-rate)+l.x*sin(-rate);
			r.y = r.y*cos(-rate)+l.y*sin(-rate);
			r.z = r.z*cos(-rate)+l.z*sin(-rate);

			l.x = l.x*cos(-rate)-r.x*sin(-rate);
			l.y = l.y*cos(-rate)-r.y*sin(-rate);
			l.z = l.z*cos(-rate)-r.z*sin(-rate);
            
			break;

        case '3':
			l.x = l.x*cos(rate)+u.x*sin(rate);
			l.y = l.y*cos(rate)+u.y*sin(rate);
			l.z = l.z*cos(rate)+u.z*sin(rate);

			u.x = u.x*cos(rate)-l.x*sin(rate);
			u.y = u.y*cos(rate)-l.y*sin(rate);
			u.z = u.z*cos(rate)-l.z*sin(rate);
            
			break;

        case '4':
			l.x = l.x*cos(-rate)+u.x*sin(-rate);
			l.y = l.y*cos(-rate)+u.y*sin(-rate);
			l.z = l.z*cos(-rate)+u.z*sin(-rate);

			u.x = u.x*cos(-rate)-l.x*sin(-rate);
			u.y = u.y*cos(-rate)-l.y*sin(-rate);
			u.z = u.z*cos(-rate)-l.z*sin(-rate);
            
			break;

        case '5':
			u.x = u.x*cos(rate)+r.x*sin(rate);
			u.y = u.y*cos(rate)+r.y*sin(rate);
			u.z = u.z*cos(rate)+r.z*sin(rate);

			r.x = r.x*cos(rate)-u.x*sin(rate);
			r.y = r.y*cos(rate)-u.y*sin(rate);
			r.z = r.z*cos(rate)-u.z*sin(rate);
            
			break;

        case '6':
			u.x = u.x*cos(-rate)+r.x*sin(-rate);
			u.y = u.y*cos(-rate)+r.y*sin(-rate);
			u.z = u.z*cos(-rate)+r.z*sin(-rate);

			r.x = r.x*cos(-rate)-u.x*sin(-rate);
			r.y = r.y*cos(-rate)-u.y*sin(-rate);
			r.z = r.z*cos(-rate)-u.z*sin(-rate);
			
            break;
        
        case ' ':
            texture_mode = 1 - texture_mode;

		default:
			break;
	}
	glutPostRedisplay();
}


void specialKeyListener(int key, int x,int y)
{
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			pos.x -= l.x;
			pos.y -= l.y;
            camera_target.x -= l.x;
            camera_target.y -= l.y;
			//pos.z += l.z;
			break;
		case GLUT_KEY_UP:		// up arrow key
			pos.x += l.x;
			pos.y += l.y;
			//pos.z-=l.z;
            camera_target.x += l.x;
            camera_target.y += l.y;
            break;

		case GLUT_KEY_RIGHT:
			// pos.x+=r.x;
			// pos.y+=r.y;
			// pos.z+=r.z;
			pos = pos + r;
            camera_target = camera_target + r;
            break;
		case GLUT_KEY_LEFT :
			// pos.x-=r.x;
			// pos.y-=r.y;
			// pos.z-=r.z;
			pos = pos - r;
            camera_target = camera_target - r;
            break;

		case GLUT_KEY_PAGE_UP:
		    pos.x+=u.x;
			pos.y+=u.y;
			pos.z+=u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
            pos.x-=u.x;
			pos.y-=u.y;
			pos.z-=u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
	glutPostRedisplay();
}

/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
    // pos.x=0;pos.y=0;pos.z=20;

    // l.x=0;l.y=0;l.z=-1;
    // u.x=0;u.y=1;u.z=0;
    // r.x=1;r.y=0;r.z=0;
    take_input();
    set_camera();
    l.x = 0, l.y = 1, l.z = 0;
    r.x = 1, r.y = 0, r.z = 0;
    u.x = 0, u.y = 0, u.z = 1;
    image = bitmap_image(image_width, image_height);
    glutInit(&argc, argv);                  // Initialize GLUT
    glutInitWindowSize(number_of_pixels, number_of_pixels);           // Set the window's initial width & height
    glutInitWindowPosition(50, 50);         // Position the window's initial top-left corner
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
    glutCreateWindow("OpenGL 3D Drawing 2");          // Create a window with the given title
    glutDisplayFunc(display);               // Register display callback handler for window re-paint
    glutReshapeFunc(reshape);               // Register callback handler for window re-shape

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);

    initGL();                               // Our own OpenGL initialization
    glutMainLoop();                         // Enter the event-processing loop
    return 0;
}