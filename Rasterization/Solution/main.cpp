#include<bits/stdc++.h>
#include"bitmap_image.hpp"
using namespace std;

#define inf 1e8
#define eps 1e-8

static unsigned long int g_seed = 1;
inline int Random()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

struct Matrix{
    vector< vector<double> > mat;
    int row, col;
    Matrix(int row, int col){
        this->row = row;
        this->col = col;
        mat.resize(row);
        for(int i=0;i<row;i++){
            mat[i].resize(col,0);
        }
    }

    Matrix operator*(Matrix const& obj)
    {
        Matrix res(row, obj.col);
        for(int i=0;i<row;i++){
            for(int j=0;j<obj.col;j++){
                for(int k=0;k<col;k++){
                    res.mat[i][j] += (mat[i][k]*obj.mat[k][j]);
                }
            }
        }
        return res;
    }

};

struct Vector{
    double x, y, z;
    Vector(double x=0, double y=0, double z=0){
        this->x =x;
        this->y =y;
        this->z =z;
    }

    Vector operator+(Vector const& obj)
    {
        Vector res;
        res.x = x+obj.x;
        res.y = y+obj.y;
        res.z = z+obj.z;
        return res;
    }

    Vector operator-(Vector const& obj)
    {
        Vector res;
        res.x = x-obj.x;
        res.y = y-obj.y;
        res.z = z-obj.z;
        return res;
    }

    Vector operator*(double d)
    {
        Vector res(x*d, y*d, z*d);
        return res;
    }

    double dot(Vector const& obj)
    {
        return x*obj.x + y*obj.y + z*obj.z;
    }

    Vector cross(Vector const& obj){
        Vector res;
        res.x = y*obj.z - z*obj.y;
        res.y = z*obj.x - x*obj.z;
        res.z = x*obj.y - y*obj.x;
        return res;
    }

    void normalize(){
        double l = sqrt(x*x+y*y+z*z);
        x = x/l;
        y = y/l;
        z = z/l;
    }

};

Vector R(Vector x, Vector a, double theta){
    Vector res =  x*cos(theta)+a*(a.dot(x))*(1-cos(theta))+a.cross(x)*sin(theta);
    return res;
}

pair<Vector, Vector> findIntersections(double y, vector<Vector> p){
    Vector lf(inf,y,inf), rt(-inf,y,inf);

    for(int i=0;i<p.size();i++){
        int j=(i+1)%p.size();
        if(abs(p[i].y-p[j].y)<eps && abs(p[i].y-y)<eps){
            if(p[i].x < lf.x) lf = Vector(p[i].x,y,p[i].z);
            if(p[j].x < lf.x) lf = Vector(p[j].x,y,p[j].z);

            if(p[i].x > rt.x) rt = Vector(p[i].x,y,p[i].z);
            if(p[j].x > rt.x) rt = Vector(p[j].x,y,p[j].z);
        }else{
            double x = (p[i].x-p[j].x)*(y - p[i].y)/(p[i].y-p[j].y) + p[i].x; 
            double z = (p[i].z-p[j].z)*(y - p[i].y)/(p[i].y-p[j].y) + p[i].z; 
            if(x>=min(p[i].x,p[j].x) && x<=max(p[i].x,p[j].x)){
                if(x < lf.x) lf = Vector(x,y,z);
                if(x > rt.x) rt = Vector(x,y,z);
            }
        }
    }

    return make_pair(lf, rt);
}

Vector eye, look, up;
double fovY, aspectRatio, near, far;
ifstream config("config.txt"),in1("scene.txt"), in2("stage1.txt"), in3("stage2.txt"),in4("stage3.txt");
ofstream out1("stage1.txt"), out2("stage2.txt"), out3("stage3.txt"),out4("z_buffer.txt");

void stage1(){
    out1<<fixed<<setprecision(7);
    in1>>eye.x>>eye.y>>eye.z>>look.x>>look.y>>look.z>>up.x>>up.y>>up.z>>fovY>>aspectRatio>>near>>far;

    stack<Matrix> st;
    Matrix M(4,4);
    M.mat[0][0] = M.mat[1][1] = M.mat[2][2] = M.mat[3][3] = 1;
    while(true){        string command;
        in1>>command;
        if(command == "triangle"){
            for(int i=0;i<3;i++){
                Matrix P(4,1);
                in1>>P.mat[0][0]>>P.mat[1][0]>>P.mat[2][0];
                P.mat[3][0]=1;
                Matrix Pp = M*P;
                Pp.mat[0][0]/=Pp.mat[3][0];
                Pp.mat[1][0]/=Pp.mat[3][0];
                Pp.mat[2][0]/=Pp.mat[3][0];
                Pp.mat[3][0]/=Pp.mat[3][0];
                out1<<Pp.mat[0][0]<<' '<<Pp.mat[1][0]<<' '<<Pp.mat[2][0]<<endl;
            }
            out1<<endl;
        }
        else if(command == "translate"){
            Matrix T(4,4);
            in1>>T.mat[0][3]>>T.mat[1][3]>>T.mat[2][3];
            T.mat[0][0] = T.mat[1][1] = T.mat[2][2] = T.mat[3][3] = 1;
            M = M*T;
        }
        else if(command == "scale"){
            Matrix T(4,4);
            in1>>T.mat[0][0]>>T.mat[1][1]>>T.mat[2][2];
            T.mat[3][3] = 1;
            M = M*T;
        }
        else if(command == "rotate"){
            double angle;
            Vector a;
            in1>>angle>>a.x>>a.y>>a.z;
            angle = (angle/180.0)*M_PI;
            a.normalize();
            Vector i(1,0,0), j(0,1,0), k(0,0,1);
            Vector c1 = R(i, a, angle);
            Vector c2 = R(j, a, angle);
            Vector c3 = R(k, a, angle);

            Matrix T(4,4);
            T.mat[0][0] = c1.x, T.mat[1][0] = c1.y, T.mat[2][0]=c1.z;
            T.mat[0][1] = c2.x, T.mat[1][1] = c2.y, T.mat[2][1]=c2.z;
            T.mat[0][2] = c3.x, T.mat[1][2] = c3.y, T.mat[2][2]=c3.z;
            T.mat[3][3] = 1;
            M = M*T;
        }
        else if(command == "push"){
            st.push(M);
        }
        else if(command == "pop"){
            M = st.top();
            st.pop();
        }
        else if(command == "end"){
            break;
        }
    }

}

void stage2(){
    out2<<fixed<<setprecision(7);
    Vector l = look - eye;
    l.normalize();
    Vector r = l.cross(up);
    r.normalize();
    Vector u = r.cross(l);

    Matrix T(4,4),R(4,4);
    T.mat[0][0]=T.mat[1][1]=T.mat[2][2]=T.mat[3][3]=1;
    T.mat[0][3]=-eye.x, T.mat[1][3]=-eye.y, T.mat[2][3]=-eye.z;
    
    R.mat[0][0]=r.x, R.mat[0][1]=r.y, R.mat[0][2]=r.z,
    R.mat[1][0]=u.x, R.mat[1][1]=u.y, R.mat[1][2]=u.z,
    R.mat[2][0]=-l.x, R.mat[2][1]=-l.y, R.mat[2][2]=-l.z, R.mat[3][3]=1;
    Matrix V = R*T;
    Matrix P(4,1);
    int cnt = 0;
    while(in2>>P.mat[0][0]>>P.mat[1][0]>>P.mat[2][0]){
        P.mat[3][0]=1;
        Matrix Pp = V*P;
        Pp.mat[0][0]/=Pp.mat[3][0];
        Pp.mat[1][0]/=Pp.mat[3][0];
        Pp.mat[2][0]/=Pp.mat[3][0];
        Pp.mat[3][0]/=Pp.mat[3][0];
        out2<<Pp.mat[0][0]<<' '<<Pp.mat[1][0]<<' '<<Pp.mat[2][0]<<endl;
        cnt++;
        if(cnt%3==0) out2<<endl;
    }

}

void stage3(){
    out3<<fixed<<setprecision(7);
    double fovX = fovY * aspectRatio;
    double t = near * tan((fovY/180.0)*M_PI/2);
    double r = near * tan((fovX/180.0)*M_PI/2);
    Matrix Prj(4,4);
    Prj.mat[0][0]=near/r,Prj.mat[1][1]=near/t,Prj.mat[2][2]=-(far+near)/(far-near), 
    Prj.mat[2][3]=-(2*far*near)/(far-near), Prj.mat[3][2]=-1;
    
    Matrix P(4,1);
    int cnt = 0;
    while(in3>>P.mat[0][0]>>P.mat[1][0]>>P.mat[2][0]){
        P.mat[3][0]=1;
        Matrix Pp = Prj*P;
        Pp.mat[0][0]/=Pp.mat[3][0];
        Pp.mat[1][0]/=Pp.mat[3][0];
        Pp.mat[2][0]/=Pp.mat[3][0];
        Pp.mat[3][0]/=Pp.mat[3][0];
        out3<<Pp.mat[0][0]<<' '<<Pp.mat[1][0]<<' '<<Pp.mat[2][0]<<endl;
        cnt++;
        if(cnt%3==0) out3<<endl;
    }
}

void stage4(){
    out4<<fixed<<setprecision(6);
    double max_x = 1.00, min_x = -1.00, max_y = 1.00, min_y = -1.00, max_z = 1.00, min_z = -1.00;

    int screen_width, screen_height;
    config>>screen_width>>screen_height;

    double **z_buffer = new double*[screen_height];
    for(int i=0;i<screen_height;i++) z_buffer[i] = new double[screen_width];
    bitmap_image *image;
    image = new bitmap_image(screen_width, screen_height);

    for(int i=0;i<screen_height;i++){
        for(int j=0;j<screen_width;j++){
            z_buffer[i][j] = max_z;
            image->set_pixel(j,i,0, 0, 0);
        }
    }

    double dx = (max_x - min_x)/screen_width;
    double dy = (max_y - min_y)/screen_height;
    double top_y = max_y - dy/2.00;
    double left_x = min_x + dx/2.00;

    vector<Vector> p(3);
    while(in4>>p[0].x>>p[0].y>>p[0].z>>p[1].x>>p[1].y>>p[1].z>>p[2].x>>p[2].y>>p[2].z){
        int color1 = Random(), color2 = Random(), color3 = Random();
        double ymax = max(p[0].y, max(p[1].y, p[2].y));
        double ymin = min(p[0].y, min(p[1].y, p[2].y));

        int r1 = max(0, (int)ceil((top_y-ymax)/dy));
        int r2 = min(499, (int)floor((top_y-ymin)/dy));
        //scan through y-axis
        for(int i=r1;i<=r2;i++){
            double ys = top_y - i * dy;
            pair<Vector, Vector> intersections = findIntersections(ys, p);
            double xa = intersections.first.x , xb = intersections.second.x, za = intersections.first.z, zb = intersections.second.z; 
            int c1 = max(0,(int)round((xa - left_x)/dx)); 
            int c2 = min(499,(int)round((xb - left_x)/dx)); 
            //scan through x-axis
            for(int j=c1;j<=c2;j++){
                double x = left_x + j* dx;
                double z;
                if(abs(xa-xb)<eps) z = za;
                else{
                    z = (za-zb)*(x-xa)/(xa-xb) + za;
                }

                if(z > min_z && z < z_buffer[i][j]){
                    z_buffer[i][j] = z;
                    image->set_pixel(j,i,color1, color2, color3);
                }
            }           
        }

    }

    for(int i=0;i<screen_height;i++){
        for(int j=0;j<screen_width;j++){
            if(z_buffer[i][j] < max_z) out4<<z_buffer[i][j]<<"\t";
        }
        out4<<endl;
    }
    image->save_image("out.bmp");

    delete image;
    for(int i=0;i<screen_height;i++) delete z_buffer[i];
    delete z_buffer;

}

int main(){
    stage1();
    stage2();
    stage3();
    stage4();
}