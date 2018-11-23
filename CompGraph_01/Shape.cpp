#include "stdafx.h"
#include "Shape.h"

#include <cmath>

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
	is << vertex.x << " " << vertex.y << " " << vertex.z;
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
	v.x = v.x - position.x;
	v.y = v.y - position.y;
	v.z = v.z - position.z;

	vertices.push_back(v);
	this->v++;
}

void Shape::addEdge(int v1, int v2)
{
	if (v <= v1 || v <= v2) return;
	edges.push_back(pair<int, int>(v1, v2));
	e++;
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
	Vertex aux;
	double sum;
	int i, j, k;
	for (j = 0; j < this->v; j++) {
		for (i = 0; i < 4; i++) {
			sum = 0;
			for (k = 0; k < 4; k++) {
				sum += vertices[j][k] * rotate[k][i];
			}
			aux[i] = sum;
		}
		vertices[j] = aux;
		//cout << aux << endl;
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
	int v0, v1, v2, p0x, p0y, p0z, p1x, p1y, p1z, px, py, pz, V;
	Vertex p0, p1;
	Vertex N, p;
	Vertex coor(x, y, z);
	for (int i = 0; i < f; i++) {
		cout << "a[0] " << vertices[edges[faces[i].edges[0]].first] << endl;
		cout << "a[1] " << vertices[edges[faces[i].edges[0]].second] << endl;
		cout << "b[0] " << vertices[edges[faces[i].edges[1]].first] << endl;
		cout << "b[1] " << vertices[edges[faces[i].edges[1]].second] << endl;
		if (edges[faces[i].edges[0]].first == edges[faces[i].edges[1]].first) {
			cout << "Triangulo " << i << " " << 1 << endl;
			p0 = vertices[edges[faces[i].edges[0]].first] - vertices[edges[faces[i].edges[0]].second];
			p1 = vertices[edges[faces[i].edges[1]].second] - vertices[edges[faces[i].edges[0]].second];
			v0 = edges[faces[i].edges[0]].first;
		}
		else if (edges[faces[i].edges[0]].first == edges[faces[i].edges[1]].second) {
			cout << "Triangulo " << i << " " << 2 << endl;
			p0 = vertices[edges[faces[i].edges[0]].first] - vertices[edges[faces[i].edges[0]].second];
			p1 = vertices[edges[faces[i].edges[1]].first] - vertices[edges[faces[i].edges[0]].second];
			
			v0 = edges[faces[i].edges[0]].first;

		}
		else if (edges[faces[i].edges[0]].second == edges[faces[i].edges[1]].first) {
			cout << "Triangulo " << i << " " << 3 << endl;
			
			p0 = vertices[edges[faces[i].edges[0]].second] - vertices[edges[faces[i].edges[0]].first];
			p1 = vertices[edges[faces[i].edges[1]].second] - vertices[edges[faces[i].edges[0]].first];
			v0 = edges[faces[i].edges[0]].first;

		}
		else {
			cout << "Triangulo " << i << " " << 4 << endl;
			p0 = vertices[edges[faces[i].edges[0]].second] - vertices[edges[faces[i].edges[0]].first];
			p1 = vertices[edges[faces[i].edges[1]].first] - vertices[edges[faces[i].edges[0]].first];
			v0 = edges[faces[i].edges[0]].first;
		}
		cout << "v0 " << vertices[v0] << endl;
		cout << "p0 " << p0 << endl;
		cout << "p1 " << p1 << endl;
		N = p0 ^ p1;
		p = (vertices[v0] + position -coor);
		cout << p << endl;
		p = vertices[v0] + position - coor;
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

		
			cout << "face: " << i << " lighted: " << (V > 0) << endl;
			cout << "face: " << i << " seno: " << seno << endl;
			cout << "face: " << i << " N: " << N << endl;
			cout << "face: " << i << " coor: " << coor << endl;
			cout << "face: " << i << " par: " << par << endl;
			cout << "face: " << i << " coordenada: " << coor << endl;
			cout << "face: " << i << " parX * parX: " << par.x * par.x << endl;
			cout << "face: " << i << " parY * parY: " << par.y * par.y << endl;
			cout << "face: " << i << " parZ * parZ: " << par.z * par.z << endl;
			cout << "face: " << i << " N.x * N.x: " << N.x * N.x << endl;
			cout << "face: " << i << " N.y * N.y: " << N.y * N.y << endl;
			cout << "face: " << i << " N.z * N.z: " << N.z * N.z << endl;
			cout << "face: " << i << " x * x: " << x * x << endl;
			cout << "face: " << i << " y * y: " << y * y << endl;
			cout << "face: " << i << " z * z: " << z * z << endl;
			cout << "face: " << i << " sqrt(par): " << raPar << endl;
			cout << "face: " << i << " sqrt(N): " << raN << endl;
			cout << "face: " << i << " sqrt(luz): " << raLuz << endl;
			cout << "face: " << i << "sqrt(N) * sqrt(luz): " << raN * raLuz << endl;
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
	int i, ver = v;
	for (i = 0; i < ver; i++) {
		addVertex(vertices[i].x + position.x + 1, vertices[i].y + position.y);
		addEdge(i, i + ver);
		vertices[i + ver].scale = vertices[i].scale;
		vertices[i + ver].z = tam - position.z;
	}
	for (i = ver; i < v - 1; i++) addEdge(i+ 1, i);
	addEdge(v - 1, ver);
	vector<int> arestas;
	for (i = 0; i < 9; i++) {
		arestas.clear();
		arestas.push_back(i);
		arestas.push_back(i + 11);
		arestas.push_back(i + 20);
		arestas.push_back(i + 10);
		addFace(arestas, 900);
	}
	arestas.clear();
	arestas = { 9, 10, 29, 19 };
	addFace(arestas, 50);
	arestas.clear();
	for (i = 20; i < 30; i++) { arestas.push_back(i); }
	addFace(arestas, 9);
	iter_swap(faces.begin(), faces.end()-1);
}

void Shape::projection(double theta)
{
	double cosT = cos(theta), sinT = sin(theta), sum = 0;
	int i, j, k;
	vector<vector<double>> projc = { {1, 0, 0, 0},
									 {0, 1, 0, 0},
									 {cosT, sinT, 0, 0},
									 {0, 0, 0, 1} };
	vector<double> temp = { 0, 0, 0, 0 };
	for (j = 0; j < v; j++) {
		for (i = 0; i < 4; i++) {
			sum = 0;
			for (k = 0; k < 4; k++) {
				//cout << k << " " << j << " " << projc[k][i] << " projection" << endl;
				sum += vertices[j][k] * projc[k][i];
			}
			temp[i] = sum;
		}
		vertices[j].x = temp[0];
		vertices[j].y = temp[1];
		vertices[j].z = temp[2];
		vertices[j].scale = temp[3];
	}
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

struct EdgeBucket
{
	int yMax; // Maximum y position of the edge
	int yMin; // Minimum y position of the edge
	int x; /* The current x position along the scan line,
	initially starting at the same point as the yMin of the edge */
	int sign; // The sign of the edge’s slope ( either -1 or 1)
	int dx; // The absolute delta x between the edge’s vertex points
	int dy; // The absolute delta y between the edge’s vertex points
	int sum; /*Initiated to zero.
	Used as the scan lines are being filled to x to the next position */

	int key;
	std::string stringValue;
	EdgeBucket(int k, const std::string& s) : key(k), stringValue(s) {}
	bool operator < (const EdgeBucket& str) const
	{
		return (key < str.key);
	}
	EdgeBucket();
};

EdgeBucket::EdgeBucket() {}

EdgeBucket addToBucket( int Xdi, int Ydi, int Xdf, int Ydf) {
	EdgeBucket nEB;
	nEB.dy = Ydi - Ydf;
	if (nEB.dy != 0)
	{
		nEB.dx = Xdi - Xdf;
		if (nEB.dx < 0)
		{
			nEB.dx = -nEB.dx;
		}
		if (nEB.dy > 0)
		{
			nEB.yMax = Ydi;
			nEB.yMin = Ydf;
			nEB.x = Xdf;
			nEB.sign = -1;
		}
		else
		{
			nEB.yMax = Ydf;
			nEB.yMin = Ydi;
			nEB.x = Xdf;
			nEB.sign = 1;
			nEB.dy = -nEB.dy;
		}
		nEB.sum = 0;
		nEB.key = nEB.yMin;
		return nEB;
	}
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
			vector<EdgeBucket> ET; // Edge Table
			vector<EdgeBucket> AL; // Active List
			for (auto edge = faces[fa].edges.begin(); edge != faces[fa].edges.end(); ++edge) {
				//cout << (*edge) << endl;
				Xdi = (int)(((this->x(edges[*edge].first) + position.x * position.scale) * width) / this->width);
				Ydi = (int)(((this->y(edges[*edge].first) + position.y * position.scale) * -height) / this->height + height);
				Xdf = (int)(((this->x(edges[*edge].second) + position.x * position.scale) * width) / this->width);
				Ydf = (int)(((this->y(edges[*edge].second) + position.y * position.scale) * -height) / this->height + height);
				//cout << Xdi << " " << Ydi << " " << Xdf << " " << Ydf << endl;
				if(mode != WIRE_FRAME)
				{
					if (Ydi != Ydf) {
						ET.push_back(addToBucket(Xdi, Ydi, Xdf, Ydf));
					}
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
			}
			std::sort(ET.begin(), ET.end());
			if (mode != WIRE_FRAME) {
				//cout << "ET empty: " << ET.empty() << endl;
				int k = 0;
				vector<vector<EdgeBucket>> divs;
				/*
				if (fa == 0 || fa == faces.size() - 1)
				{
					int i = edges[faces[fa].edges[0]].first, j = i + 9;
					k = 14;
					ET.clear();
					while (k> 0) {
						Xdi = (int)(((this->x(i) + position.x * position.scale) * width) / this->width);
						Ydi = (int)(((this->y(i) + position.y * position.scale) * -height) / this->height + height);
						Xdf = (int)(((this->x(j) + position.x * position.scale) * width) / this->width);
						Ydf = (int)(((this->y(j) + position.y * position.scale) * -height) / this->height + height);
						ET.push_back(addToBucket(Xdi, Ydi, Xdf, Ydf));
						if (k == 14) { j = i + 2; i += 9; }
						if (k == 13) { j = i; i--; }
						if (k == 12) { j = i; i--; }
						if (k == 11) { 
							divs.push_back(ET);
							ET.clear();
							j += 9; i ++; 
						}
						if (k == 10) { i = j; j--; }
						if (k == 9) { i--; j--; }
						if (k == 8) { i--; j--; }
						if (k == 7) { i = j; j -= 3; }
						if (k == 6) { i = j; j--; }
						if (k == 5) { 
							divs.push_back(ET);
							ET.clear();
							j += 4; }
						if (k == 4) { i = j; j--; }
						if (k == 3) { i--; j--; }
						if (k == 2) { i--; j--; }
						if (k == 1) {
							divs.push_back(ET);
							ET.clear();
						}
						k--;
					}
					k = 3;
				}
				//*/
				if (faces[fa].isLighted && mode == SOLID_WITH_LIGHT) {
					//cout << "Iluminado " << fa << endl;
					SDL_SetRenderDrawColor(renderer, 255 * faces[fa].seno, 165 * faces[fa].seno, 0 * faces[fa].seno, SDL_ALPHA_OPAQUE);
					//cout <<"face: " << fa << " islight: " << faces[fa].seno << endl;
				}
				else {
					//cout << "Apagado " << fa << endl;
					int r, g, b;
					r = faces[fa].color % 100;
					if (r > 9) { r = 9; }
					g = faces[fa].color % 10 - r * 10;
					b = faces[fa].color - g * 10;
					SDL_SetRenderDrawColor(renderer, 255 * r / 9, 255 * g / 9, 255 * b / 9, SDL_ALPHA_OPAQUE);
					//SDL_SetRenderDrawColor(renderer, 128, 128, 128, SDL_ALPHA_OPAQUE);
				}
				while (k >= 0) {
					if (k > 0) { ET.swap(divs[k-1]); }

					int y = ET[0].yMin;
					while (!ET.empty()) {
						if (!AL.empty()) {
							int i = 0;
							int j = AL.size();
							while (i < j) {
								if (AL[i].yMax == y) {
									AL.erase(AL.begin() + i);
									j--;
								}
								else { i++; }
							}
							i = 0;
							j = ET.size();
							while (i < j) {
								if (ET[i].yMax == y) {
									ET.erase(ET.begin() + i);
									j--;
								}
								else { i++; }
							}
						}
						for (int i = 0; i < ET.size(); ++i) {
							if (ET[i].yMin == y)
							{
								AL.push_back(ET[i]);
								AL[AL.size() - 1].key = AL[AL.size() - 1].x;
							}
						}
						std::sort(AL.begin(), AL.end());
						//cout << AL.size() << endl;
						int tam = AL.size();
						for (int i = 0; i < tam - 1; i++) {
							//cout << AL[i].x << endl;
							SDL_RenderDrawLine(
								renderer,
								AL[i].x, y,
								AL[i + 1].x, y
							);
						}
						for (int i = 0; i < AL.size(); ++i) {
							if (AL[i].dx != 0) {
								AL[i].sum += AL[i].dx;
							}
							while (AL[i].sum >= AL[i].dy) {
								AL[i].x += AL[i].sign;
								AL[i].sum -= AL[i].dy;
							}
						}
						y++;
					}
					if (k == 1) { k--; }
					k--;
				}
			}
		}
	}
}

void Shape::metaShape(SDL_Renderer* renderer, int width, int height)
{
	
	int Xdi, Ydi;
	for (int i = 0; i < v; i++) {
		Xdi = (int)(((this->x(i) + position.x * position.scale) * width) / this->width);
		Ydi = (int)(((this->y(i) + position.y * position.scale) * -height) / this->height + height);
		SDL_SetRenderDrawColor(renderer, 255/v * i, 255 / v * i, 255 / v * i, SDL_ALPHA_OPAQUE); // ... (renderer, r, g, b, alpha)
		SDL_Rect rect = { Xdi - 2, Ydi - 2, 4, 4 };
		SDL_RenderFillRect(renderer, &rect);
		//cout << i << "x:" << vertices[i].x << " y:" << vertices[i].y << " s:" << vertices[i].scale << endl;
		//cout << i << "xd:" << Xdi << " yd:" << Ydi << endl;

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
