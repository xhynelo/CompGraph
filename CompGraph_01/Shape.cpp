#include "stdafx.h"
#include "Shape.h"

#include <cmath>
#include <algorithm>
#include <map>

Shape::Shape()
{
	v = 0;
	e = 0;
	f = 0;
	matrix = {
			{1.0, 0.0, 0.0, 0.0},
			{0.0, 1.0, 0.0, 0.0},
			{0.0, 0.0, 1.0, 0.0},
			{0.0, 0.0, 0.0, 1.0}
	};
}

double Shape::x(int n)
{
	if (v <= n) return -1;
	return vertices[n].x * vertices[n].scale;
}

double Shape::y(int n)
{
	if (n >= v) return -1;
	return vertices[n].y * vertices[n].scale;
}

double Shape::z(int n)
{
	if (n >= v) return -1;
	return vertices[n].z * vertices[n].scale;
}

double Shape::getV()
{
	return v;
}

Vertex::Vertex(double mX, double mY) : x(mX), y(mY), z(0), scale(1) {}
Vertex::Vertex(double mX, double mY, double mZ) : x(mX), y(mY), z(mZ), scale(1) {}

Vertex::Vertex() {}

bool Vertex::operator ==(const Vertex &b) const 
{
	if (x == b.x && y == b.y && z == b.z) {
		return true;
	}
	else {
		return false;
	}
}

Vertex Vertex::operator-(const Vertex &vertex) const {
	return Vertex(x - vertex.x, y - vertex.y, z - vertex.z);
}

Vertex Vertex::operator+(const Vertex &vertex) const {
	return Vertex(x + vertex.x, y + vertex.y, z + vertex.z);
}

Vertex Vertex::operator^(const Vertex &vertex) const {
	return Vertex(
		y * vertex.z - z * vertex.y,
		z * vertex.x - x * vertex.z,
		x * vertex.y - y * vertex.x
	);
}

double Vertex::operator*(const Vertex &vertex) const {
	return x * vertex.x + y * vertex.y + z * vertex.z;
}


Vertex Vertex::operator*(const double &v) const {
	return Vertex(x*v, y*v, z*v);
}

Vertex Vertex::operator%(const Vertex &v) const {
	return Vertex(
		y * v.z - v.y * z,
		x * v.z - v.x * z,
		x * v.y - v.x * y
	);
}

double& Vertex::operator[](int i) {
	switch (i)
	{
	case 0: return x;
	case 1: return y;
	case 2: return z;
	default: return scale;
	}
}

double Vertex::modulo() {
	return sqrt(*this * *this);
}

std::istream & operator>>(std::istream & is, Vertex & vertex)
{
	is >> vertex.x >> vertex.y;
	return is;
}

std::ostream & operator<<(std::ostream & is, Vertex & vertex)
{
	is << vertex.x << " " << vertex.y << " " << vertex.z << " " << vertex.scale;
	return is;
}

void Shape::addVertex(double x, double y)
{
	if (vXmin < 0 || x < this->x(vXmin)) vXmin = v;
	if (vYmin < 0 || y < this->y(vYmin)) vYmin = v;
	if (vXMAX < 0 || x > this->x(vXMAX)) vXMAX = v;
	if (vYMAX < 0 || y > this->y(vYMAX)) vYMAX = v;
	addVertex(Vertex(x, y));
}

void Shape::addVertex(Vertex v)
{
	if (!this->v) position = v;
	v.x = v.x;
	//v.x = v.x - position.x;
	v.y = v.y;
	//v.y = v.y - position.y;
	v.z = v.z;
	//v.z = v.z - position.z;

	verticesSalvo.push_back(v);
	vertices.push_back(v);
	this->v++;
}

void Shape::addEdge(int v1, int v2)
{
	if (v <= v1 || v <= v2) return;
	edges.push_back(pair<int, int>(v1, v2));
	e++;
}

void Shape::addEdgeCurva(int v1, int v2, Vertex v3, Vertex v4) {
	if (v <= v1 || v <= v2) return;
	addEdge(v1, v2);
	curvas[e - 1] = pair<Vertex, Vertex>(v3, v4);
	curvasSalvo[e - 1] = pair<Vertex, Vertex>(v3, v4);
}

Face::Face() : color(0), isVisible(true), isLighted(false) {}

void Shape::addFace(vector<int> lados, int color)
{
	Face fa;
	fa.edges = lados;
	fa.color = color;
	f++;
	faces.push_back(fa);
}

void Shape::setSRU(double x, double y)
{
	height = y;
	width = x;
}

void Shape::scale(double n)
{
	for (int i = 0; i < v; i++) {
		vertices[i].scale = n;
	}
}

/*
void Shape::tranform()
{
	vector<Vertex> matrixAux = vertices;
	//cout << "matrix" << endl;
	for (auto it = matrix.begin(); it != matrix.end(); it++) {
		for (auto it2 = it->begin(); it2 != it->end(); it2++) {
			//cout << *it2 << " ";
		}
		//cout << endl;
	}
	//cout << ";" << endl;
	//cout << "verticesPrint" << endl;
	for (auto it = verticesPrint.begin(); it != verticesPrint.end(); it++) {
		//cout << *it << endl;
	}
	//cout << ";" << endl;
	double sum;
	int i, j, k;
	for (j = 0; j < this->v; j++) {
		for (i = 0; i < 4; i++) {
			sum = 0;
			for (k = 0; k < 4; k++) {
				if (k < 3) {
					sum += matrixAux[j][k] * matrix[k][i];
				}
				else {
					sum += (matrixAux[j][k] / matrixAux[j][k]) * matrix[k][i];
				}
			}
			Vertex &v = verticesPrint[j];
			v[i] = sum;
		}
	}
	//cout << "verticesPrint" << endl;
	for (auto it = verticesPrint.begin(); it != verticesPrint.end(); it++) {
		//cout << *it << endl;
	}
	//cout << ";" << endl;
}
//*/

//*
void Shape::translade(double x, double y)
{
	position.x += x;
	position.y += y;
}

double Shape::minVertex() {
	double minimo = -1;
	Vertex v1;
	for (int i = 0; i < edges.size(); i++) {
		v1 = vertices[edges[i].second] - vertices[edges[i].first];
		if (v1.modulo() < minimo || minimo == -1) {
			minimo = v1.modulo();
		}
	}
	return minimo;
}



void Shape::rotate(double theta)
{
	theta = theta / 180.0 * M_PI;
	//cout << theta << endl;
	vector<vector<double>> rotate = {
		{cos(theta), sin(theta), 0, 0 },
		{-sin(theta), cos(theta), 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}
	};
	//cout << "sqrt(theta): " << sqrt(theta) << endl;
	Vertex aux, aux2;
	double sum, duba;
	int i, j, k;
	for (j = 0; j < this->v; j++) {
		for (i = 0; i < 4; i++) {
			sum = 0;
			for (k = 0; k < 4; k++) {
				sum += verticesSalvo[j][k] * rotate[k][i];
			}
			aux[i] = sum;
		}
		vertices[j] = aux + position;
		//cout << aux << endl;
	}

	for (auto w = curvas.begin(); w != curvas.end(); w++) {
		for (i = 0; i < 4; i++) {
			sum = 0;
			duba = 0;
			for (k = 0; k < 4; k++) {
				sum += w->second.first[k] * rotate[k][i];
				duba += w->second.second[k] * rotate[k][i];
			}
			aux[i] = sum;
			aux[2] = duba;
		}
		w->second.first = aux;
		w->second.first = aux2;
	}
}

void Shape::reseta() {
	vertices = verticesSalvo;
	curvas = curvasSalvo;
}

Vertex quartenionRotate(double s, Vertex &v, Vertex &r) {
	return (
		(r * (s*s))
		- (r * (v*v))
		+ (v*(v*r * 2.0)
		+ ((v ^ r) * s * 2.0)
	));
}

void Shape::rotateQ(double theta, int x, int y, int z) {
	theta = (theta / 180.0) * M_PI;
	pair<double, Vertex> q, q_;
	double s;
	Vertex ve(x, y, z);
	s = cos(theta / 2);
	ve = ve * sin(theta / 2);

	int i = 0;
	for (i = 0; i < v; i++) {
		vertices[i] = quartenionRotate(s, ve, vertices[i]) + position;
	}
	for (auto it = curvasSalvo.begin(); it != curvasSalvo.end(); it++) {
		curvas[it->first] = pair<Vertex, Vertex>(
			quartenionRotate(s, ve, it->second.first) + position,
			quartenionRotate(s, ve, it->second.second) + position
		);
	}
}

void Shape::rotate(double thetaX, double thetaY, double thetaZ){
	vector<vector<double>> rotate_;
	if (thetaX != 0) {
		rotate(thetaX);
	}
	if (thetaY != 0){
		thetaZ = thetaZ / 180.0 * M_PI;
		//cout << theta << endl;
		rotate_ = {
			{cos(thetaY), 0, sin(thetaY), 0},
			{0, 1, 0, 0 },
			{-sin(thetaY), 0, cos(thetaY), 0},
			{0, 0, 0, 1}
		};
		//cout << "sqrt(theta): " << sqrt(theta) << endl;
		Vertex aux;
		double sum;
		int i, j, k;
		for (j = 0; j < this->v; j++) {
			for (i = 0; i < 4; i++) {
				sum = 0;
				for (k = 0; k < 4; k++) {
					sum += vertices[j][k] * rotate_[k][i];
				}
				aux[i] = sum;
			}
			vertices[j] = aux;
			//cout << aux << endl;
		}
	}
	if (thetaZ != 0) {
		thetaZ = thetaZ / 180.0 * M_PI;
		//cout << theta << endl;
		rotate_ = {
			{1, 0, 0 ,0},
			{0, cos(thetaZ), -sin(thetaZ), 0 },
			{0, sin(thetaZ), cos(thetaZ), 0},
			{0, 0, 0, 1}
		};
		//cout << "sqrt(theta): " << sqrt(theta) << endl;
		Vertex aux;
		double sum;
		int i, j, k;
		for (j = 0; j < this->v; j++) {
			for (i = 0; i < 4; i++) {
				sum = 0;
				for (k = 0; k < 4; k++) {
					sum += vertices[j][k] * rotate_[k][i];
				}
				aux[i] = sum;
			}
			vertices[j] = aux;
			//cout << aux << endl;
		}
	}
}
//*/

void Shape::setPosition(double x, double y)
{
	position.x = x;
	position.y = y;
}

void Shape::hider(double x, double y, double z, bool islight)
{
	double p0x, p0y, p0z, p1x, p1y, p1z, px, py, pz, V;
	int v0, v1, v2;
	Vertex p0, p1;
	Vertex N, p;
	Vertex coor(x, y, z);
	for (int i = 0; i < f; i++) {
		//cout << "a[0] " << vertices[edges[faces[i].edges[0]].first] << endl;
		//cout << "a[1] " << vertices[edges[faces[i].edges[0]].second] << endl;
		//cout << "b[0] " << vertices[edges[faces[i].edges[1]].first] << endl;
		//cout << "b[1] " << vertices[edges[faces[i].edges[1]].second] << endl;
		if (edges[faces[i].edges[0]].first == edges[faces[i].edges[1]].first) {
			//cout << "Triangulo " << i << " " << 1 << endl;
			p0 = vertices[edges[faces[i].edges[0]].first] - vertices[edges[faces[i].edges[0]].second];
			p1 = vertices[edges[faces[i].edges[1]].second] - vertices[edges[faces[i].edges[0]].second];
			v0 = edges[faces[i].edges[0]].first;
		}
		else if (edges[faces[i].edges[0]].first == edges[faces[i].edges[1]].second) {
			//cout << "Triangulo " << i << " " << 2 << endl;
			p0 = vertices[edges[faces[i].edges[0]].first] - vertices[edges[faces[i].edges[0]].second];
			p1 = vertices[edges[faces[i].edges[1]].first] - vertices[edges[faces[i].edges[0]].second];
			
			v0 = edges[faces[i].edges[0]].first;

		}
		else if (edges[faces[i].edges[0]].second == edges[faces[i].edges[1]].first) {
			//cout << "Triangulo " << i << " " << 3 << endl;
			
			p0 = vertices[edges[faces[i].edges[0]].second] - vertices[edges[faces[i].edges[0]].first];
			p1 = vertices[edges[faces[i].edges[1]].second] - vertices[edges[faces[i].edges[0]].first];
			v0 = edges[faces[i].edges[0]].first;

		}
		else {
			//cout << "Triangulo " << i << " " << 4 << endl;
			p0 = vertices[edges[faces[i].edges[0]].second] - vertices[edges[faces[i].edges[0]].first];
			p1 = vertices[edges[faces[i].edges[1]].first] - vertices[edges[faces[i].edges[0]].first];
			v0 = edges[faces[i].edges[0]].first;
		}
		//cout << "v0 " << vertices[v0] << endl;
		//cout << "p0 " << p0 << endl;
		//cout << "p1 " << p1 << endl;
		N = p0 ^ p1;
		//p = (vertices[v0] + position -coor);
		p = (vertices[v0] -coor);
		//cout << p << endl;
		p = vertices[v0] - coor;
		//p = vertices[v0] + position - coor;
		V = N * p;
		if (!islight) {
			if (V >= 0) {
				faces[i].isVisible = true;
			}
			else {
				faces[i].isVisible = false;
			}
		}
		else {
			double parX, parY, parZ, seno, raPar, raN, raLuz;
			Vertex par = N ^ coor;
			//raPar = sqrt((parX * parX) + (parY * parY) + (parZ * parZ));
			raPar = N * coor;
			raN = N.modulo();
			raLuz = coor.modulo();
			seno = raPar / (raN * raLuz);
			if (raN == 0) {
				seno = 1;
			}

		
			//cout << "face: " << i << " lighted: " << (V > 0) << endl;
			//cout << "face: " << i << " seno: " << seno << endl;
			//cout << "face: " << i << " N: " << N << endl;
			//cout << "face: " << i << " coor: " << coor << endl;
			//cout << "face: " << i << " par: " << par << endl;
			//cout << "face: " << i << " coordenada: " << coor << endl;
			//cout << "face: " << i << " parX * parX: " << par.x * par.x << endl;
			//cout << "face: " << i << " parY * parY: " << par.y * par.y << endl;
			//cout << "face: " << i << " parZ * parZ: " << par.z * par.z << endl;
			//cout << "face: " << i << " N.x * N.x: " << N.x * N.x << endl;
			//cout << "face: " << i << " N.y * N.y: " << N.y * N.y << endl;
			//cout << "face: " << i << " N.z * N.z: " << N.z * N.z << endl;
			//cout << "face: " << i << " x * x: " << x * x << endl;
			//cout << "face: " << i << " y * y: " << y * y << endl;
			//cout << "face: " << i << " z * z: " << z * z << endl;
			//cout << "face: " << i << " sqrt(par): " << raPar << endl;
			//cout << "face: " << i << " sqrt(N): " << raN << endl;
			//cout << "face: " << i << " sqrt(luz): " << raLuz << endl;
			//cout << "face: " << i << "sqrt(N) * sqrt(luz): " << raN * raLuz << endl;
			seno = abs(seno);
			if (V >= 0) {
				faces[i].isLighted = true;
				faces[i].seno = seno;
			}
			else {
				faces[i].isLighted = false;
				faces[i].seno = seno;
			}
		}
		
	}
}

void Shape::slide(double tam)
{
	int i, j, ver = v, fac = f, edg = e;
	for (i = 0; i < ver; i++) {
		addVertex(vertices[i].x, vertices[i].y);
		//addVertex(vertices[i].x + position.x , vertices[i].y + position.y);

		vertices[i + ver].scale = vertices[i].scale;
		vertices[i + ver].z = tam;
		//vertices[i + ver].z = tam - position.z;
		verticesSalvo[i + ver].z = tam;
		//verticesSalvo[i + ver].z = tam - position.z;
	}

	for (i = 0; i < edg; i++) {
		addEdge(edges[i].first + 8, edges[i].second + 8);
		if (curvasSalvo.find(i) != curvasSalvo.end()) {
			curvasSalvo[e - 1] = curvasSalvo[i];
			curvas[e - 1] = curvas[i];

			curvasSalvo[e - 1].first = curvasSalvo[e - 1].second;
			curvasSalvo[e - 1].second = curvasSalvo[i].first;
			curvas[e - 1].first = curvas[e - 1].second;
			curvas[e - 1].second = curvas[i].first;
			
			curvasSalvo[e - 1].first.z = tam - position.z;
			curvasSalvo[e - 1].second.z = tam - position.z;
			curvas[e - 1].first.z = tam - position.z;
			curvas[e - 1].second.z = tam - position.z;
		}
	}


	vector<int> arestas;
	//cout << faces.size() << endl;
	for (i = 0; i < fac; i++) {
		arestas.clear();
		for (auto it2 = faces[i].edges.begin(); it2 != faces[i].edges.end(); ++it2) {
			arestas.push_back(*it2 + 8);
		}
		std::reverse(arestas.begin(), arestas.end());
		addFace(arestas, 900);
	}
	

	
	int c = 16;
	addEdge(0, 0 + 8);
	for (i = 1; i < ver; ++i) {
		//addEdge(i - 1, i + 10);
		addEdge(i, i + 8);
		
		arestas.clear();
		if (i == ver - 1) {
			arestas.push_back(i + 7);
			arestas.push_back(c + 1);
			arestas.push_back(i - 1);
			arestas.push_back(c);
		}
		else {
			arestas.push_back(c);
			arestas.push_back(i - 1);
			arestas.push_back(c + 1);
			arestas.push_back(i + 7);
		}
	
		addFace(arestas, 900);
		
		c += 1;
		/*
		arestas.clear();
		arestas.push_back(c);
		arestas.push_back(i - 1);
		arestas.push_back(c + 1);
		addFace(arestas, 900);

		c += 1;
		*/
	}
	arestas = { 7, 16, 15, 23 };
	addFace(arestas, 900);
	/*
	addEdge(9, 10);
	arestas = { 53, 26, 52 };
	addFace(arestas, 900);

	arestas = { 34, 53, 9 };
	addFace(arestas, 900);
	
	*/
//arestas = { 9, 10, 29, 19 };
	//addFace(arestas, 50);
	//arestas.clear();
	//for (i = 20; i < 30; i++) { arestas.push_back(i); }
	//addFace(arestas, 9);
	iter_swap(faces.begin(), faces.end()-1);
}

void projectionMult(vector<vector<double>> &projc, Vertex &v) {
	vector<double> temp = { 0, 0, 0, 0 };
	for (int i = 0; i < 4; i++) {
		double sum = 0;
		for (int k = 0; k < 4; k++) {
			//cout << k << " " << j << " " << projc[k][i] << " projection" << endl;
			sum += v[k] * projc[k][i];
		}
		temp[i] = sum;
	}
	v.x = temp[0];
	v.y = temp[1];
	v.z = temp[2];
	v.scale = temp[3];
}

void Shape::projection(double theta)
{
	degree = M_PI * (theta / 180.0);
	//*
	double cosT = cos(degree), sinT = sin(degree), sum = 0;
	int i, j, k;
	vector<vector<double>> projc = { {1, 0, 0, 0},
									 {0, 1, 0, 0},
									 {cosT, sinT, 0, 0},
									 {0, 0, 0, 1} };
	for (j = 0; j < v; j++) {
		projectionMult(projc, vertices[j]);
	}
	for (auto it = curvas.begin(); it != curvas.end(); it++) {
		projectionMult(projc, it->second.first);
		projectionMult(projc, it->second.second);
	}
	//*/
}

/*
void Shape::printShape(SDL_Renderer* renderer, int width, int height)
{
	int w, h;
	int Xdi, Xdf, Ydi, Ydf;
	SDL_GetRendererOutputSize(renderer, &w, &h);
	//cout << "printShape" << endl;
	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		Xdi = (int) (((this->x(edge->first) + position.x * position.scale) * width) / this->width);
		Ydi = (int) (((this->y(edge->first) + position.y * position.scale) * -height) / this->height + height);
		Xdf = (int) (((this->x(edge->second) + position.x * position.scale) * width) / this->width);
		Ydf = (int) (((this->y(edge->second) + position.y * position.scale) * -height) / this->height + height);
		//cout << Xdi << " " << Ydi << " " << Xdf << " " << Ydf << endl;
		SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE); // ... (renderer, r, g, b, alpha)
		SDL_RenderDrawLine(
			renderer,
			Xdi, Ydi,
			Xdf, Ydf
		);
	}
}
//*/



int _gfxPrimitivesCompareInt(const void *a, const void *b)
{
	// Code from SDL gfx: https://github.com/keera-studios/SDL2_gfx
	return (*(const int *)a) - (*(const int *)b);
}

int hline(SDL_Renderer * renderer, Sint16 x1, Sint16 x2, Sint16 y)
{
	// Code from SDL gfx: https://github.com/keera-studios/SDL2_gfx
	return SDL_RenderDrawLine(renderer, x1, y, x2, y);;
}

int filledPolygonRGBAMT(SDL_Renderer * renderer, const Sint16 * vx, const Sint16 * vy, int n, Uint8 r, Uint8 g, Uint8 b, Uint8 a, int **polyInts, int *polyAllocated)
{
	// Code from SDL gfx: https://github.com/keera-studios/SDL2_gfx
	int result;
	int i;
	int y, xa, xb;
	int miny, maxy;
	int x1, y1;
	int x2, y2;
	int ind1, ind2;
	int ints;
	int *gfxPrimitivesPolyInts = NULL;
	int *gfxPrimitivesPolyIntsNew = NULL;
	int gfxPrimitivesPolyAllocated = 0;

	/*
	* Vertex array NULL check
	*/
	if (vx == NULL) {
		return (-1);
	}
	if (vy == NULL) {
		return (-1);
	}

	/*
	* Sanity check number of edges
	*/
	if (n < 3) {
		return -1;
	}


	static int *gfxPrimitivesPolyIntsGlobal = NULL;
	static int gfxPrimitivesPolyAllocatedGlobal = 0;
	/*
	* Map polygon cache
	*/
	if ((polyInts == NULL) || (polyAllocated == NULL)) {
		/* Use global cache */
		gfxPrimitivesPolyInts = gfxPrimitivesPolyIntsGlobal;
		gfxPrimitivesPolyAllocated = gfxPrimitivesPolyAllocatedGlobal;
	}
	else {
		/* Use local cache */
		gfxPrimitivesPolyInts = *polyInts;
		gfxPrimitivesPolyAllocated = *polyAllocated;
	}

	/*
	* Allocate temp array, only grow array
	*/
	if (!gfxPrimitivesPolyAllocated) {
		gfxPrimitivesPolyInts = (int *)malloc(sizeof(int) * n);
		gfxPrimitivesPolyAllocated = n;
	}
	else {
		if (gfxPrimitivesPolyAllocated < n) {
			gfxPrimitivesPolyIntsNew = (int *)realloc(gfxPrimitivesPolyInts, sizeof(int) * n);
			if (!gfxPrimitivesPolyIntsNew) {
				if (!gfxPrimitivesPolyInts) {
					free(gfxPrimitivesPolyInts);
					gfxPrimitivesPolyInts = NULL;
				}
				gfxPrimitivesPolyAllocated = 0;
			}
			else {
				gfxPrimitivesPolyInts = gfxPrimitivesPolyIntsNew;
				gfxPrimitivesPolyAllocated = n;
			}
		}
	}

	/*
	* Check temp array
	*/
	if (gfxPrimitivesPolyInts == NULL) {
		gfxPrimitivesPolyAllocated = 0;
	}

	/*
	* Update cache variables
	*/
	if ((polyInts == NULL) || (polyAllocated == NULL)) {
		gfxPrimitivesPolyIntsGlobal = gfxPrimitivesPolyInts;
		gfxPrimitivesPolyAllocatedGlobal = gfxPrimitivesPolyAllocated;
	}
	else {
		*polyInts = gfxPrimitivesPolyInts;
		*polyAllocated = gfxPrimitivesPolyAllocated;
	}

	/*
	* Check temp array again
	*/
	if (gfxPrimitivesPolyInts == NULL) {
		return(-1);
	}

	/*
	* Determine Y maxima
	*/
	miny = vy[0];
	maxy = vy[0];
	for (i = 1; (i < n); i++) {
		if (vy[i] < miny) {
			miny = vy[i];
		}
		else if (vy[i] > maxy) {
			maxy = vy[i];
		}
	}

	/*
	* Draw, scanning y
	*/
	result = 0;
	for (y = miny; (y <= maxy); y++) {
		ints = 0;
		for (i = 0; (i < n); i++) {
			if (!i) {
				ind1 = n - 1;
				ind2 = 0;
			}
			else {
				ind1 = i - 1;
				ind2 = i;
			}
			y1 = vy[ind1];
			y2 = vy[ind2];
			if (y1 < y2) {
				x1 = vx[ind1];
				x2 = vx[ind2];
			}
			else if (y1 > y2) {
				y2 = vy[ind1];
				y1 = vy[ind2];
				x2 = vx[ind1];
				x1 = vx[ind2];
			}
			else {
				continue;
			}
			if (((y >= y1) && (y < y2)) || ((y == maxy) && (y > y1) && (y <= y2))) {
				gfxPrimitivesPolyInts[ints++] = ((65536 * (y - y1)) / (y2 - y1)) * (x2 - x1) + (65536 * x1);
			}
		}

		qsort(gfxPrimitivesPolyInts, ints, sizeof(int), _gfxPrimitivesCompareInt);

		/*
		* Set color
		*/
		result = 0;
		result |= SDL_SetRenderDrawBlendMode(renderer, (a == 255) ? SDL_BLENDMODE_NONE : SDL_BLENDMODE_BLEND);
		result |= SDL_SetRenderDrawColor(renderer, r, g, b, a);

		for (i = 0; (i < ints); i += 2) {
			xa = gfxPrimitivesPolyInts[i] + 1;
			xa = (xa >> 16) + ((xa & 32768) >> 15);
			xb = gfxPrimitivesPolyInts[i + 1] - 1;
			xb = (xb >> 16) + ((xb & 32768) >> 15);
			result |= hline(renderer, xa, xb, y);
		}
	}

	return (result);
}


void Shape::printShape(SDL_Renderer* renderer, int width, int height, int mode) // Print usando faces
{
	int w, h, fa;
	int Xdi, Xdf, Ydi, Ydf;
	SDL_GetRendererOutputSize(renderer, &w, &h);
	//cout << "printShape" << endl;
	//cout << "e ae" << faces.size() << f << endl;
	//faces[11].isVisible = false;
	//faces[7].isVisible = false;

	//cout << "Entrei" << endl;
	for (fa = 0; fa < f; fa++) {
		//cout << fa << endl;
		if (faces[fa].isVisible) {
			//cout << "entrou" << endl;
			vector<Vertex> novo = sortFace(faces[fa]);

			Vertex ultimo = novo[novo.size() - 1];

			int tam = (int)novo.size();
			Sint16 *vx = new Sint16[tam];
			Sint16 *vy = new Sint16[tam];


			for (int i = 0; i < novo.size(); i++) {
				Vertex &atual = novo[i];
				//*
				//Xdi = (int)((((ultimo.x + ultimo.z * cos(degree)) + position.x * position.scale) * width) / this->width);
				//Ydi = (int)((((ultimo.y + ultimo.z * sin(degree)) + position.y * position.scale) * -height) / this->height + height);
				//Xdf = (int)((((atual.x + atual.z * cos(degree)) + position.x * position.scale) * width) / this->width);
				//Ydf = (int)((((atual.y + atual.z * sin(degree)) + position.y * position.scale) * -height) / this->height + height);
				Xdi = (int)((((ultimo.x) * ultimo.scale) * width) / this->width);
				Ydi = (int)((((ultimo.y) * ultimo.scale) * -height) / this->height + height);
				Xdf = (int)((((atual.x) * atual.scale) * width) / this->width);
				Ydf = (int)((((atual.y) * atual.scale) * -height) / this->height + height);


				//*/
				/*
				Xdi = (int)(((ultimo.x + position.x * position.scale) * width) / this->width);
				Ydi = (int)(((ultimo.y + position.y * position.scale) * -height) / this->height + height);
				Xdf = (int)(((atual.x + position.x * position.scale) * width) / this->width);
				Ydf = (int)(((atual.y + position.y * position.scale) * -height) / this->height + height);
				//*/
				if (mode != WIRE_FRAME)
				{
					vx[i] = (Sint16)Xdf;
					vy[i] = (Sint16)Ydf;
				}
				if (mode == WIRE_FRAME)
				{
					SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE); // ... (renderer, r, g, b, alpha)
					SDL_RenderDrawLine(
						renderer,
						Xdi, Ydi,
						Xdf, Ydf
					);
				}
				ultimo = atual;
			}

			if (mode != WIRE_FRAME) {
				int r = 0, g = 0, b = 0;
				if (faces[fa].isLighted && mode == SOLID_WITH_LIGHT) {
					//cout << "Iluminado " << fa << endl;
					r = (int)(255 * faces[fa].seno);
					g = (int)(127 * faces[fa].seno);
					b = (int)0;
					
					//SDL_SetRenderDrawColor(renderer,, , 0, SDL_ALPHA_OPAQUE);//cosseno disfarcado
					//cout <<"face: " << fa << " islight: " << faces[fa].seno << endl;
				}
				else {
					//cout << "Apagado " << fa << endl;
					r = (int)(255 * r / 9);
					g = (int)(255 * g / 9);
					b = (int)(255 * b / 9);
					//SDL_SetRenderDrawColor(renderer,, , , SDL_ALPHA_OPAQUE);
				}
				Uint32 color = 0;
				Uint8 *c = (Uint8 *)&color;
				filledPolygonRGBAMT(renderer, vx, vy, tam, r, g, b, SDL_ALPHA_OPAQUE, NULL, NULL);

			}
		}

	}
}

void Shape::metaShape()
{
	
	Vertex Xdi, Ydi;
	for (int i = 0; i < f; i++) {
		cout << i;
		for (int j = 0; j < faces[i].edges.size(); j++) {
			cout << " " << faces[i].edges[j];
		}
		//Xdi = vertices[i] + position;
		cout << endl;

	}
}

void Shape::readShape(string name)
{
	vertices.clear();
	edges.clear();
	v = 0;
	e = 0;
	scale(1);
	ifstream ifs(name, ifstream::in);
	if (!ifs.bad())
	{
		int v_aux;
		ifs >> v_aux;
		for (int i = 0; i < v_aux; i++) {
			Vertex vertex;
			ifs >> vertex;
			//cout << vertex << endl;
			addVertex(vertex);
		}
		int e_aux;
		ifs >> e_aux;
		for (int i = 0; i < e_aux; i++) {
			int v1, v2;
			ifs >> v1 >> v2;
			addEdge(v1, v2);
		}
		int s;
		ifs >> s;
		scale(s);
		ifs.close();
	}
}

vector<vector<double>> Shape::matrixMult(vector<vector<double>> mat1, vector<vector<double>> mat2)
{
	vector<vector<double>> matrix = mat1;
	double sum;
	int i, j, k;
	for (j = 0; j < 4; j++) {
		for (i = 0; i < 4; i++) {
			sum = 0;
			for (k = 0; k < 4; k++) {
				sum += mat1[j][k] * mat2[k][i];
			}
			matrix[j][i] = sum;
		}
	}
	return matrix;
}

vector<Vertex> Shape::sortFace(Face poly) {
	vector<Vertex> vertics;
	int steps = 30;
	int vertic, primeiro;
	if (edges[poly.edges[0]].first == edges[poly.edges[1]].first) {
		vertic = edges[poly.edges[0]].second;
	}
	else if (edges[poly.edges[0]].first == edges[poly.edges[1]].second) {
		vertic = edges[poly.edges[0]].second;
	}
	else if (edges[poly.edges[0]].second == edges[poly.edges[1]].first) {
		vertic = edges[poly.edges[0]].first;
	}
	else {
		vertic = edges[poly.edges[0]].first;
	}
	Vertex p1 = vertices[vertic];
	primeiro = vertic;
	vertics.push_back(vertices[primeiro]);
	for (int i = 0; i < poly.edges.size(); i++) {
		Vertex p2;

		if (vertic == edges[poly.edges[i]].first) {
			p2 = vertices[edges[poly.edges[i]].second];
			vertic = edges[poly.edges[i]].second;
		}
		else {
			p2 = vertices[edges[poly.edges[i]].first];
			vertic = edges[poly.edges[i]].first;
		}
		if (curvas.find(poly.edges[i]) != curvas.end()) {
			//Vertex t1 = curvas[poly.edges[i]].first + position;
			//Vertex t2 = curvas[poly.edges[i]].second + position;
			Vertex t1 = curvas[poly.edges[i]].first;
			Vertex t2 = curvas[poly.edges[i]].second;
			for (int t = 0; t <= steps; t++) {
				double s = (double)t / (double)steps;
				double s1 = 1 - s;
				/*

				double h1 = - s*s*s + 3 * s*s -3 *s + 1;
				double h2 = 3 * s*s*s - 6 * s*s +3 *s;
				double h3 = -3 * s*s*s + 3 * s*s;
				double h4 = 3 * s*s*s;
				double h1 = 2 * s*s*s - 3 * s*s + 1;
				double h2 = -2 * s*s*s + 3 * s*s;
				double h3 = s * s*s - 2 * s*s + s;
				double h4 = s * s*s - s * s;
				
				//
				Vertex p = p1 * h1 + p2 * h2 + t2 * h3 + t1 * h4;
				*/
				//
				Vertex p = p1 * s1*s1*s1 + t1 * 3 * s1*s1*s + t2 * 3 * s1*s*s + p2 * s*s*s;
				//p.z = p1.z;
				vertics.push_back(p);
				//cout << poly.edges[i] << " interpolado " << p << endl;
			}
		}
		else {
			vertics.push_back(p2);
			//cout << "padrao " << p2 << endl;
		}
		p1 = p2;

	}
	vertics.push_back(vertices[primeiro]);
	return vertics;
}
