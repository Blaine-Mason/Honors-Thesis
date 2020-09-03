#include <math.h>
#include <stdlib.h>

#include <vector>
#include <iostream>

using namespace std;

float max_angle = 30.0;
float max_length = 40.0;
float min_length = 0.8;

double uniform(double start, double end) {
  return start + (rand()/(double)(RAND_MAX))*(end - start);
}

void randomVector(float *v) {
    //si parte da 1, 1, 1 si somma x e y con lunghezza minore di tan(max_angle)
    float m = 0.5*tan(max_angle*(float)M_PI/180.0f);
    v[0] = 1.0f + uniform(-m, m);
    v[1] = 1.0f + uniform(-m, m);
    v[2] = 1.0f + uniform(-m, m);

    float length = uniform(min_length, max_length);
    //rinormalizziamo a length
    float r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])*length;
    v[0] /= r;
    v[1] /= r;
    v[2] /= r;
}

int main(int argc, char *argv[]) {


    if(argc != 4) {
        cout << "Usage: " << argv[0] << "number_of_generators max_angle max_length\n\n";
        return 0;
    }

    vector<float> points;

    int n = atoi(argv[1]);
    max_angle = atoi(argv[2]);
    max_length = atoi(argv[3]);


    points.resize(n*3);

    points.resize(n*3);
    for(int i = 0; i < n; i++) {
        float *v = &points[i*3];
        randomVector(v);
        cout << v[0] << " " << v[1] << " " << v[2] << endl;
    }
    return 1;
}
