#include<bits/stdc++.h>
using namespace std;
#define PI acos(-1)
class point {
    public: 
        double x, y, z; 
        double n; 
        point() {
            x = y = z = 0;
            n = 1.0;
        }
        point(double _x, double _y, double _z) {
            x = _x; y = _y; z = _z;
            n = 1.0;
        }
        point(double _x, double _y, double _z, double _n) {
            x = _x; y = _y; z = _z; n = _n;
        }
        point operator + (point p) {
            return point(x + p.x, y + p.y, z + p.z);
        }
        point operator - (point p) {
            return point(x - p.x, y - p.y, z - p.z);
        }
        point operator * (double k) {
            return point(x * k, y * k, z * k);
        }
        point operator / (double k) {
            return point(x / k, y / k, z / k);
        }
        double operator *(point b)  {
            return x * b.x + y* b.y + z * b.z;
        }
        point operator ^(point p) {
            return point(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
        }
        void normalize() {
            double r = sqrtl(x * x + y * y + z * z);
            x /= r; y /= r; z /= r;
        }
        double dot(point p) {
            return x * p.x + y * p.y + z * p.z;
        }
        point cross(point p) {
            return point(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
        }
        void scale(double n, int type) {
            if(type == 1) {
                x *= n; y *= n; z *= n;
            }
            else {
                x /= n; y /= n; z /= n;
            }
        }
        void scaleDown() {
            x /= n; y /= n; z /= n;
        }
};
class Triangle {
    public:
        point pt1, pt2, pt3;
        Triangle(point _pt1, point _pt2, point _pt3) {
            pt1 = _pt1; pt2 = _pt2; pt3 = _pt3;
        }
};

class Matrix {
    public: 
        vector<vector<double>> mat;
        int dimension_cnt;
        Matrix() {
            dimension_cnt = 4;
            mat.resize(dimension_cnt);
            for(int i = 0; i < dimension_cnt; i++) {
                mat[i].resize(dimension_cnt);
            }
        }
        Matrix(int _dimension_cnt) {
            dimension_cnt = _dimension_cnt;
            mat.resize(dimension_cnt);
            for(int i = 0; i < dimension_cnt; i++) {
                mat[i].resize(dimension_cnt);
            }
        }
        Matrix(int _dimension_cnt, bool is_identitiy) {
            dimension_cnt = _dimension_cnt;
            mat.resize(dimension_cnt);
            for(int i = 0; i < dimension_cnt; i++) {
                mat[i].resize(dimension_cnt);
                for(int j = 0; j < dimension_cnt; j++) {
                    mat[i][j] = (i == j);
                }
            }
        }
        Matrix operator * (Matrix m) {
            Matrix res(dimension_cnt);
            for(int i = 0; i < dimension_cnt; i++) {
                for(int j = 0; j < dimension_cnt; j++) {
                    for(int k = 0; k < dimension_cnt; k++) {
                        res.mat[i][j] += mat[i][k] * m.mat[k][j];
                    }
                }
            }
            return res;
        }
        point operator * (point b) {
            double res[4];
            for(int i = 0; i < 4; i++) res[i] = 0;

            double vec[4] = {b.x, b.y, b.z, b.n};

            for(int i = 0; i < 4; i++) {
                res[i] = 0;
                for(int j = 0; j < 4; j++) {
                    res[i] += mat[i][j]* vec[j];
                }
            }
            point ret(res[0], res[1], res[2], res[3]);
            ret.scaleDown();
            return ret;
        }

        void print() {
            for(int i = 0; i < dimension_cnt; i++) {
                for(int j = 0; j < dimension_cnt; j++) {
                    cout << mat[i][j] << " ";
                }
                cout << endl;
            }
        }

};
Matrix createTranslation(int dim, point pt) {
        Matrix res(dim, true);
        res.mat[0][3] = pt.x;
        res.mat[1][3] = pt.y;
        res.mat[2][3] = pt.z;
        return res;
    }

Matrix createScale(int dim, point pt) {
    Matrix res(dim, true);
    res.mat[0][0] = pt.x;
    res.mat[1][1] = pt.y;
    res.mat[2][2] = pt.z;
    return res;
}
point applyRodriguesFormula(point x, point a, double theta) {   
    double rad = theta * PI / 180;
    return x * cos(rad) + a * (a * x) * (1 - cos(rad)) + (a ^ x) * sin(rad);
}
Matrix createRotation(point p, double theta) {
    Matrix res(4, true);
    point a(p.x,p.y,p.z);
    a.normalize();
    point x(1,0,0); // unit vector along x axis
    point y(0,1,0); // unit vector along y axis
    point z(0,0,1); // unit vector along z axis

    point c1 = applyRodriguesFormula(x, a, theta);
    point c2 = applyRodriguesFormula(y, a, theta);
    point c3 = applyRodriguesFormula(z, a, theta);

    res.mat[0][0] = c1.x;
    res.mat[1][0] = c1.y;
    res.mat[2][0] = c1.z;

    res.mat[0][1] = c2.x;
    res.mat[1][1] = c2.y;
    res.mat[2][1] = c2.z;

    res.mat[0][2] = c3.x;
    res.mat[1][2] = c3.y;
    res.mat[2][2] = c3.z;
    return res;
}

Matrix createView(point eye, point look, point up) {
    point z = (look-eye);
    z.normalize();
    point x = (z^up);
    x.normalize();
    point y = (x^z);
    y.normalize();

    Matrix res = Matrix(4,true);
    res.mat[0][0] = x.x;
    res.mat[0][1] = x.y;
    res.mat[0][2] = x.z;

    res.mat[1][0] = y.x;
    res.mat[1][1] = y.y;
    res.mat[1][2] = y.z;

    res.mat[2][0] = -z.x;
    res.mat[2][1] = -z.y;
    res.mat[2][2] = -z.z;

    Matrix t = createTranslation(4, point(-eye.x,-eye.y,-eye.z));

    res = t * res;
    return res;

}
Matrix createProjection(double fovY, double aspect, double near, double far) {
    Matrix res = Matrix(4, true);
    double fovX = fovY * aspect;
    double t = near * tan(((fovY/2) * PI) / 180.0);
    double r = near * tan(((fovX/2) * PI) / 180.0);
    
    res.mat[0][0] = near/r;
    res.mat[1][1] = near/t;
    res.mat[2][2] = -(far + near)/(far - near);
    res.mat[3][2] = -1;
    res.mat[2][3] = -(2.0 * far * near)/(far - near);
    res.mat[3][3] = 0;

    return res;

}

int main() {
    
    ifstream fin("scene.txt");
    ofstream fout("stage1.txt");
    point eye, look, up;
    double fovY, aspect, near, far;
    int total_triangle = 0;
    fin >> eye.x >> eye.y >> eye.z;
    fin >> look.x >> look.y >> look.z;
    fin >> up.x >> up.y >> up.z;
    fin >> fovY >> aspect >> near >> far;
    stack<Matrix> st;
    Matrix curr(4, true);
    st.push(curr);
    int flag = 0;
    while(1) {
        string txt; fin >> txt;
        if(txt == "triangle") {
            point p1, p2, p3; 
            fin >> p1.x >> p1.y >> p1.z;
            fin >> p2.x >> p2.y >> p2.z;
            fin >> p3.x >> p3.y >> p3.z;
            point p1t = st.top() * p1;
            point p2t = st.top() * p2;
            point p3t = st.top() * p3;
            fout << setprecision(7) << fixed;
            fout << p1t.x << " " << p1t.y << " " << p1t.z << endl;
            fout << p2t.x << " " << p2t.y << " " << p2t.z << endl;
            fout << p3t.x << " " << p3t.y << " " << p3t.z << endl;
            fout << endl;
            total_triangle++;
        }
        else if(txt == "translate") {
            double tx, ty, tz;
            fin >> tx >> ty >> tz;
            Matrix t = createTranslation(4, point(tx,ty,tz));
            if(flag == 0) {
                for(int i = 0; i < 4; i++) {
                    for(int j = 0; j < 4; j++) {
                        cout << t.mat[i][j] << " ";
                    }
                    cout << endl;
                }
            }
            Matrix res = st.top() * t;
            st.pop();
            st.push(res);
            if(flag == 0) {
                for(int i = 0; i < 4; i++) {
                    for(int j = 0; j < 4; j++) {
                        cout << res.mat[i][j] << " ";
                    }
                    cout << endl;
                }
            
            }
            flag = 1;
        }
        else if(txt == "scale") {
            double sx, sy, sz;
            fin >> sx >> sy >> sz;
            Matrix t = createScale(4, point(sx,sy,sz));
            Matrix res = st.top() * t;
            st.pop();
            st.push(res);
        }
        else if(txt == "rotate") {
            double angle, rx, ry, rz;
            fin >> angle >> rx >> ry >> rz;
            Matrix t = createRotation(point(rx,ry,rz), angle);
            Matrix res = st.top() * t;
            st.pop();
            st.push(res);
        }
        else if(txt == "push") {
            st.push(st.top());
        }
        else if(txt == "pop") {
            st.pop();
        }
        else if(txt == "end") {
            break;
        }
        else {
            cout << "Invalid input" << endl;
            break;
        }
    }
    fin.close();
    fout.close();

    fin.open("stage1.txt");
    fout.open("stage2.txt");

    Matrix view = createView(eye, look, up);
    for(int i = 0; i < total_triangle; i++) {
        point p1, p2, p3;
        fin >> p1.x >> p1.y >> p1.z;
        fin >> p2.x >> p2.y >> p2.z;
        fin >> p3.x >> p3.y >> p3.z;
        point p1t = view * p1;
        point p2t = view * p2;
        point p3t = view * p3;
        fout << setprecision(7) << fixed;
        fout << p1t.x << " " << p1t.y << " " << p1t.z << endl;
        fout << p2t.x << " " << p2t.y << " " << p2t.z << endl;
        fout << p3t.x << " " << p3t.y << " " << p3t.z << endl;
        fout << endl;
    }

    fin.close();
    fout.close();

    fin.open("stage2.txt");
    fout.open("stage3.txt");

    Matrix projection = createProjection(fovY, aspect, near, far);
    for(int i = 0; i < total_triangle; i++) {
        point p1, p2, p3;
        fin >> p1.x >> p1.y >> p1.z;
        fin >> p2.x >> p2.y >> p2.z;
        fin >> p3.x >> p3.y >> p3.z;
        point p1t = projection * p1;
        point p2t = projection * p2;
        point p3t = projection * p3;

        fout << setprecision(7) << fixed;
        fout << p1t.x << " " << p1t.y << " " << p1t.z << endl;
        fout << p2t.x << " " << p2t.y << " " << p2t.z << endl;
        fout << p3t.x << " " << p3t.y << " " << p3t.z << endl;
        fout << endl;
    }

    fin.close();
    fout.close();
}