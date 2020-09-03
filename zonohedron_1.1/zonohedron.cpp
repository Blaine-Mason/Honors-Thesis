#include <math.h>
#include <assert.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
using namespace std;

//returns -1, 0, 1 depending on sign of val
template <typename T> int sign(T val) {
	return (T(0) < val) - (val < T(0));
}

//simple 3d vector class
struct Vec3 {
	double v[3];

	Vec3() {}
	Vec3(double x, double y, double z) { v[0] = x; v[1] = y; v[2] = z; }
	double &operator[](int i) { return v[i]; }
	double operator[](int i) const { return v[i]; }

	Vec3 operator+(const Vec3 &b) const { return Vec3(v[0]+b[0], v[1] + b[1], v[2] + b[2]); }
	Vec3 operator-(const Vec3 &b) const { return Vec3(v[0]-b[0], v[1] - b[1], v[2] - b[2]); }

	Vec3 operator*(double s) const { return Vec3(v[0]*s, v[1]*s, v[2]*s); }
	double operator*(const Vec3 &a) const { return v[0]*a[0] + v[1]*a[1] + v[2]*a[2]; }

	Vec3 &operator+=(const Vec3 &a) { v[0] += a[0]; v[1] += a[1]; v[2] += a[2]; return *this; }
	Vec3 &operator-=(const Vec3 &a) { v[0] -= a[0]; v[1] -= a[1]; v[2] -= a[2]; return *this; }
	Vec3 operator-() { return Vec3(-v[0], -v[1], -v[2]); }

	Vec3 &operator*=(double s) { v[0] *= s; v[1] *= s; v[2] *= s; return *this; }

	bool operator!=(const Vec3 &a) const { return v[0] != a[0] || v[1] != a[1] || v[2] != a[2]; }

	Vec3 operator^(const Vec3 &b) const {
		Vec3 c;
		c[0] = v[1]*b[2] - v[2]*b[1];
		c[1] = v[2]*b[0] - v[0]*b[2];
		c[2] = v[0]*b[1] - v[1]*b[0];
		return c;
	}
	double squaredNorm() { return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; }
	double norm() { return sqrt(squaredNorm()); }
	void normalize() { (*this) *= 1/norm(); }

	//project on the plane 1, 1, 1 passing by n, then project on xy plane.
/*	bool operator<(const Vec3 &b) const {
		Vec3 n(1, 1, 1);
		n *= 1/n.norm();
		Vec3 A = (*this) * (1/((*this)*n)) - n;
		Vec3 B = b * (1/(b*n)) - n;
		return atan2(A[1], A[0]) < atan2(B[1], B[0]);
	}*/
	bool operator<(const Vec3 &b) const {
		if(v[0] == b[0]) {
			if(v[1] == b[1]) {
				return v[2] < b[2];
			}
			return v[1] < b[1];
		}
		return v[0] < b[0];
	}
	bool operator==(const Vec3 &b) const {
		return v[0] == b[0] && v[1] == b[1] && v[2] == b[2];
	}
};

//compute solid angle of the 3 vectors
double solidAngle(Vec3 &a, Vec3 &b, Vec3 &c) {
	double A = a.norm();
	double B = b.norm();
	double C = c.norm();
	return 2*atan( (a^b)*c / (A*B*C + (a*b)*C + (a*c)*B + (b*c)*A));
}

//compute solid angle for a fan of vertices (all must be on the boundary and sorted).
double solidAngle(vector<Vec3> &star, vector<int> &edges) {
	double solid_angle = 0;

	Vec3 n(1.0, 1.0, 1.0);
	for(int i = 0; i < edges.size(); i += 2) {
		double s = solidAngle(n, star[edges[i]], star[edges[i+1]]);
		solid_angle += s;
	}
	return  solid_angle;
}

//unique sign for degenerate vertex (careful, will fail with more than 2 identical vertices)
int triSign(int i, int j, int k,  vector<Vec3> &star) {
	Vec3 a = star[i] + Vec3(1, i, i*i)*1e05;
	Vec3 b = star[j] + Vec3(1, j, j*j)*1e05;
	Vec3 c = star[k] + Vec3(1, k, k*k)*1e05;
	return sign((a^b)*c);
}

Vec3 project(vector<int> &point, vector<Vec3> &star) {
	Vec3 v(0, 0, 0);
	int n = star.size();
	for(int i = 0; i < n; i++) {
		v += star[i]*point[i];
	}
	return v;
}

void readTxt(const char *filename, vector<Vec3> &star) {
	std::ifstream input(filename);
	Vec3 v;
	while(input >> v[0]) {
		input >> v[1];
		input >> v[2];
		star.push_back(v);
	}

}

double zonohedron(vector<Vec3> &star, vector<Vec3> &vertices, vector<int> &edges) {

	int n = star.size();

	vector<int> face(n);
	//vectors lying on the solid angle at the origin of the zonohedron
	edges.clear();

	double volume = 0;

	for(int i = 0; i < n; i++) {
		for(int j = i+1; j < n; j++) {
			Vec3 normal = star[i]^star[j];
			int count = 0; //count number of generator on une side of this plane.
			for(int k = 0; k < n; k++) {
				if(k == i || k == j)
					face[k] = 0;
				else {
					face[k] = sign(normal * star[k]);
					//in case of degenerate faces, this will decompose it correctly. (why? I just guessed)
					if(face[k] == 0) {
						cout << "Warning degenerate face. " << endl;
						face[k] = triSign(i, j, k, star); //sign(k - i)* sign(k - j);
					}
				}
				if(face[k] > 0) count++;
			}
			if(count == 0) {
				edges.push_back(j);
				edges.push_back(i);
			}
			if(count == n-2) {
				edges.push_back(i);
				edges.push_back(j);
			}

			//find 4 vertices
			face[i] = -1;
			face[j] = -1;  //-1, -1
			vertices.push_back(project(face, star));

			face[i] = 1;   //1, -1
			vertices.push_back(project(face, star));

			face[j] = 1;   //1, 1
			vertices.push_back(project(face, star));

			face[i] = -1;   //-1 1
			vertices.push_back(project(face, star));

			volume += vertices.back()*normal;

		}
	}

	//add other half of vertices
	int nvert = vertices.size();
	for(int i = nvert-1; i >= 0; i--)  //otherwise vertex order in faces get inverted
		vertices.push_back(vertices[i]*-1);

	volume /= 3; //they are all cones.
	return volume;
}


/*

#set dgrid3d
#
set pm3d

set style line 1 lt 4 lw .5
set pm3d at s hidden3d 1
unset surf
set zrange[0: 1]
set palette defined ( 0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red", 5 "red")
splot 'toytest.txt' notitle

#select view direction using mouse

set output toygini.pdf
set term pdfcairo
splot 'toytest.txt' notitle
unset output
#unset flushes buffer (stupid gnuplot)
#(or svg)

*/

void study(vector<Vec3> &star, Vec3 range, double xstep, double ystep) {

	double z = range[2];

	star.push_back(Vec3(0, 0, 0));
	for(float x = 0; x < range[0]; x += 0.5) {
		for(float y = 0; y < 10; y += 0.5) {

			star.back() = Vec3(x, y, z);

			vector<Vec3> vertices;
			vector<int> edges;

			double volume = zonohedron(star, vertices, edges);

			Vec3 diagonal(0, 0, 0);
			double total_length_squared = 0;
			for(int i = 0; i < star.size(); i++) {
				diagonal += star[i];
				total_length_squared += star[i].squaredNorm();
			}
			double a = diagonal[2] / diagonal.norm();
			double gini = volume / (diagonal[0] * diagonal[1] * diagonal[2]);
			cout << x << ", " << y << ", " << gini << endl;
		}
		cout << endl;
	}
}

//intersects edge v0-v1 at z provided by middle (where the result is stored.
bool intersects(const Vec3 &v0, const Vec3 &v1, Vec3 &middle) {
	double z = middle[2];
	if(v0[2] > z && v1[2] > z) return false;
	if(v0[2] < z && v1[2] < z) return false;
	double len = v0[2] - v1[2]; //might be negative or zero.
	if(fabs(len) < 0.00001) {
		middle = (v0 + v1)*0.5;
		middle[2] = z;
		return true;
	}
	double u = v0[2] - z;
	double w = z - v1[2];
	middle[0] = (u * v1[0] + w*v0[0])*(1/len);
	middle[1] = (u * v1[1] + w*v0[1])*(1/len);
	return true;
}


double B0(double t) {
	t = 1 - t;
	return t*t*t/6.0;
}

double B1(double t) {
	return (3.0*t*t*t - 6.0*t*t + 4.0)/6.0;
}

double B2(double t) {
	return (-3.0*t*t*t + 3.0*t*t + 3.0*t + 1.0)/6.0;
}

double B3(double t) {
	return t*t*t/6.0;
}


double binomial(int n, int k) {
	double result = 1;
	for(int i = 0; i < k; i++)
		result *= ((double)n - i) / ((double)i + 1);
	return result;
}

//t goes from = control[0] to 1 = controls[0] after a full circle.;
Vec3 bspline3(double x, std::vector<Vec3> &controls) {
	int n = controls.size();
	//find integer part of t and fractional

	int p = floor(x*n);
	double t = x*n - p;

	Vec3 P0 = controls[(p + n-1)%n];
	Vec3 P1 = controls[p];
	Vec3 P2 = controls[(p+1)%n];
	Vec3 P3 = controls[(p+2)%n];

	Vec3 pos = P0*B0(t) + P1*B1(t) + P2*B2(t) + P3*B3(t);

	return pos;
}

Vec3 movingAverage(double x, std::vector<Vec3> &controls) {
	int n = controls.size();
	//find integer part of t and fractional

	//max distance:
	int p = floor(x*n);
	double t = x*n - p;
	Vec3 target = controls[p]*(1-t) + controls[(p + 1)%n]*t;

	double maxd = 0.0;
	std::vector<double> distances;
	for(int i = 0; i < n; i++) {
		double d = (controls[i] - target).norm();
		if(d > maxd) maxd = d;
		distances.push_back(d);
	}
	double sigma = maxd/10;
	double totw = 0;
	Vec3 pos(0, 0, 0);
	for(int i = 0; i < n; i++) {
		double d = distances[i];
		double w = exp(-d*d/(2*sigma*sigma));
		pos += controls[i] * w;
		totw += w;
	}
	pos *= 1/totw;
	return pos;
}


//level is the fractional height of he zonohedron
void isoquants(char *file, double level, vector<Vec3> &vertices, vector<int> &edges) {

	//find max (min is 0)/
	double maxz = 0.0;
	double maxx = 0.0;
	double maxy = 0.0;
	for(int i = 0; i < vertices.size(); i++) {
		if(vertices[i][0] > maxx)
			maxx = vertices[i][0];
		if(vertices[i][1] > maxy)
			maxy = vertices[i][1];
		if(vertices[i][2] > maxz)
			maxz = vertices[i][2];
	}

	double z = -maxz + 2*maxz*level;
	vector<Vec3> line;
	Vec3 middle(0, 0, z);
	for(int i = 0; i < vertices.size(); i+= 4) {
		Vec3 &v0 = vertices[i+0];
		Vec3 &v1 = vertices[i+1];
		Vec3 &v2 = vertices[i+2];
		Vec3 &v3 = vertices[i+3];
		if(intersects(v0, v1, middle))
			line.push_back(middle);
		if(intersects(v1, v2, middle))
			line.push_back(middle);
		if(intersects(v2, v3, middle))
			line.push_back(middle);
		if(intersects(v3, v0, middle))
			line.push_back(middle);
	}
	for(int i = 0; i < line.size(); i++) {
		line[i][0] += maxx;
		line[i][1] += maxy;
	}

	//sort points in line using centroid and angle
	Vec3 c(0, 0, 0);
	for(int i = 0; i < line.size(); i++)
		c += line[i];
	c *= 1.0/line.size();

	sort(line.begin(), line.end(), [c](const Vec3 &a, const Vec3 &b) { return atan2(-a[1] + c[1], -a[0] + c[0]) < atan2(-b[1] + c[1], -b[0] + c[0]); });

	//remove duplicated points
	int count = 0;
	for(int i = 1; i < line.size(); i++) {
		if(line[count] != line[i]) {
			count++;
			line[count] = line[i];
		}
	}
	line.resize(count+1);

	vector<Vec3> optimal;
	int start = false;
	for(int i = 1; i < line.size(); i++) {
		if(line[i][1] <= line[i-1][1] && line[i][0] >= line[i-1][0] ) {
			start = true;
			optimal.push_back(line[i]);
		} else {
			if(start == true) break;
		}
	}

	if(optimal.size() < 5)
		return;

	vector<double> sms(optimal.size()); //saggio marginale sostituzione
	for(int i = 1; i < optimal.size()-1; i++) {
		double dx = optimal[i+1][0] - optimal[i-1][0];
		double dy = optimal[i+1][1] - optimal[i-1][1];
		assert(dx >= 0);
		assert(dy <= 0);
		sms[i] += -dy / dx;
	}

	vector<double> es(optimal.size()); //elasticita' di sostituzione
	for(int i = 2; i < optimal.size()-2; i++) {
		double ds = sms[i+1] - sms[i-1];
		double dds = ds / sms[i];
		double dx = (optimal[i+1][1]/optimal[i+1][0]) - (optimal[i-1][1]/optimal[i-1][0]);
		double ddx = dx / (optimal[i][1]/optimal[i][0]);
		es[i] = ddx/dds;
	}

	//smooth es:
	vector<double> smoothed(es.size());
	int side = 20;
	for(int i = 2; i < es.size()-2; i++) {
		double tot = 0.0;
		double x = 0.0;
		double y = 0.0;
		for(int k = i-side; k <= i+side; k++) {
			if(k < 2 || k > es.size()-2) continue;
			double w = fabs(side+1 - k)/(side+1);
			tot += w;
			x += w * es[k];
		}
		smoothed[i] = x/tot;
	}

	//save points as csv
	ofstream fileb;
	fileb.open("isoquant_" + std::to_string(z) + ".csv");
	for(int i = 2; i < optimal.size()-2; i++) {
		fileb << optimal[i][0] << "," << optimal[i][1] << "," << sms[i] << "," << smoothed[i] << "\n";
	}

	/*	std::vector<Vec3> interpolated;

	for(int i = 0; i < 1000; i++) {
		interpolated.push_back(movingAverage(i/1000.0, optimal));
	}

	ofstream file1;
	file1.open("isoint_" + std::to_string(z) + ".csv");
	for(int i = 0; i < interpolated.size(); i++) {
		file1 << interpolated[i][0] << "," << interpolated[i][1] << "\n";
		all << interpolated[i][0] << "," << interpolated[i][1] << "\n";
	} */
}

int main(int argc, char *argv[]) {

	if(argc != 2 && argc != 3 && argc != 8) {
		cout << "Usage: " << argv[0] << " <file> [output.obj]\n\n";
		cout << " The format of the file is 3 numbers (%f %f %f) for every line.\n";
		cout << " If you want to evaluate gini index when adding a new vector use:\n";
		cout << argv[0] << " <file> gini <x range> <y range> <z value> <x step> <y step>\n";
		cout << " An array of x,y,gini(x,y) will be written to terminal.\n";
		return -1;
	}

	vector<Vec3> star;

	readTxt(argv[1], star);
	if(star.size() < 3) {
		cout << "Need at least 3 vectors, found " << star.size() << endl;
		return 0;
	}

	if(argc > 2 && string(argv[2]) == "gini") {
		if(argc != 8) {
			cerr << "Expecting parameters as: <file> gini <x range> <y range> <z value> <x step> <y step>\n";
			return -1;
		}
		double x_range= atof(argv[3]);
		double y_range= atof(argv[4]);
		double z_value= atof(argv[5]);
		double x_step= atof(argv[6]);
		double y_step= atof(argv[7]);
		study(star, Vec3(x_range, y_range, z_value), x_step, y_step);
		return 0;
	}
	vector<Vec3> vertices;
	vector<int> edges;

	double volume = zonohedron(star, vertices, edges);

	isoquants("test", 0.5, vertices, edges);

	cout << "Total volume: " << volume << endl;

	Vec3 diagonal(0, 0, 0);
	double total_length_squared = 0;
	for(int i = 0; i < star.size(); i++) {
		diagonal += star[i];
		total_length_squared += star[i].squaredNorm();
	}

	cout << "Diagonal: " << diagonal[0] << " " << diagonal[1] << " " << diagonal[2] << endl;
	cout << "Diagonal norm: " << diagonal.norm() << endl;
	cout << "Cosine against output: " << diagonal[2] / diagonal.norm()  << endl;
	cout << "Tangent of angle diagonal - plane x-y: "
		 << diagonal[2]/sqrt(pow(diagonal[0], 2) + pow(diagonal[1], 2)) << endl;
	cout << "Tangent against input axes: " << diagonal[1]/diagonal[0] << endl;
	//cout << "Cosine of projection of diagonal on plane x-y with x axis: "
	//     << diagonal[0]/sqrt(diagonal[0]*diagonal[0] + diagonal[1]*diagonal[1]) << endl;
	cout << "Volume against diagonal cubed: " << volume/ pow(diagonal.norm(), 3) << endl;
	cout << "Sum of squared norms: " << total_length_squared << endl;
	//	cout << "Mystery number: " << binomial(star.size(), 3) * pow(total_length_squared/3.0, 1.5) << endl;

	double gini = volume / (diagonal[0] * diagonal[1] * diagonal[2]);
	cout << "Gini 3D index: " << gini << endl;

	cout << "Solid angle: " << solidAngle(star, edges) << endl;

	vector<Vec3> nstar = star;
	vector<Vec3> nvertices;
	vector<int> nedges;

	for(int i = 0; i < nstar.size(); i++)
		nstar[i].normalize();
	double normalized_volume = zonohedron(nstar, nvertices, nedges);
	cout << "Normalized vectors volume: " << normalized_volume << endl;

	vector<Vec3> span; //vectors on the boundary of the origin solid cone.
	for(int i = 0; i < edges.size(); i += 2)
		span.push_back(star[edges[i]]);

	vector<Vec3> e_vertices;
	float e_volume = zonohedron(span, e_vertices, edges);
	Vec3 e_diagonal(0, 0, 0);
	for(int i = 0; i < span.size(); i++)
		e_diagonal += span[i];

	cout << "Volume against diagonal cubed of boundary vectors: " << e_volume / pow(e_diagonal.norm(), 3) << endl;

	ofstream boundary;
	boundary.open ("boundary.txt");
	for(int i = 0; i < span.size(); i++) {
		Vec3 &v = span[i];
		boundary << v[0] << " " << v[1] << " " << v[2] << endl;
	}
	boundary.close();

	if(argc == 2) return 0; //no output

	ofstream file;
	file.open(argv[2]);



	int n = star.size();
	vector<int> faces(n*(n-1)*4); //4 indexes per face
	for(unsigned int i = 0; i < faces.size(); i++)
		faces[i] = i;

	//optional unification of vertices. This MIGHT workz if based on floating point computations...

	for(unsigned int i = 0; i < vertices.size(); i++) {
		Vec3 v = vertices[i];
		v += diagonal;
		v *= 0.5;
		file << "v " << v[0] << " " << v[1] << " " << v[2] << endl;
	}

	bool quads = true;
	if(quads) {
		//write down obj
		for(int i = 0; i < n*(n-1); i++)
			file  << "f " << i*4+1 << " " << i*4 + 2 << " " << i*4 + 3 << " " << i*4 + 4 << endl;
	} else {
		for(int i = 0; i < n*(n-1); i++) {
			file  << "f " << i*4+1 << " " << i*4 + 2 << " " << i*4 + 3 << endl;
			file  << "f " << i*4+3 << " " << i*4 + 4 << " " << i*4 + 1 << endl;
		}
	}

	file.close();
	return 0;
}

