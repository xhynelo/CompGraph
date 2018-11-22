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

double& Vertex::operator[](int i) {
	switch (i)
	{
	case 0: return x;
	case 1: return y;
	case 2: return z;
	default: return scale;
	}
}

std::istream & operator>>(std::istream & is, Vertex & vertex)
{
	is >> vertex.x >> vertex.y;
	return is;
}

std::ostream & operator<<(std::ostream & is, Vertex & vertex)
{
	is << vertex.x << " " << vertex.y;
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

void Shape::addFace(vector<int> lados)
{
	Face fa;
	fa.edges = lados;
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

void Shape::hider()
{
	int v0, v1, v2, p0x, p0y, p0z, p1x, p1y, p1z, px, py, pz, V;
	Vertex N;
	for (int i = 0; i < f; i++) {
		if (edges[faces[i].edges[0]].first == edges[faces[i].edges[1]].first) {
			v0 = edges[faces[i].edges[0]].second;
			v1 = edges[faces[i].edges[0]].first;
			v2 = edges[faces[i].edges[1]].second;
		}
		else if (edges[faces[i].edges[0]].first == edges[faces[i].edges[1]].second) {
			v0 = edges[faces[i].edges[0]].second;
			v1 = edges[faces[i].edges[0]].first;
			v2 = edges[faces[i].edges[1]].first;
		}
		else if (edges[faces[i].edges[0]].second == edges[faces[i].edges[1]].first) {
			v0 = edges[faces[i].edges[0]].first;
			v1 = edges[faces[i].edges[0]].second;
			v2 = edges[faces[i].edges[1]].second;
		}
		else {
			v0 = edges[faces[i].edges[0]].first;
			v1 = edges[faces[i].edges[0]].second;
			v2 = edges[faces[i].edges[1]].first;
		}
		p0x = (vertices[v0].x + position.x) - (vertices[v1].x + position.x);
		p0y = (vertices[v0].y + position.y) - (vertices[v1].y + position.y);
		p0z = (vertices[v0].z + position.z) - (vertices[v1].z + position.z);
		p1x = (vertices[v2].x + position.x) - (vertices[v1].x + position.x);
		p1y = (vertices[v2].y + position.y) - (vertices[v1].y + position.y);
		p1z = (vertices[v2].z + position.z) - (vertices[v1].z + position.z);
		N.x = p0y * p1z - p0z * p1y;
		N.y = p0x * p1z - p0z * p1x;
		N.z = p0x * p1y - p0y * p1x;
		px = (vertices[v0].x + position.x) - 0;
		py = (vertices[v0].y + position.x) - 0;
		pz = (vertices[v0].z + position.x) - 100;
		V = px * N.x + py * N.y + pz * N.z;
		if (V >= 0) {
			faces[i].isVisible = true;
		}
		else {
			faces[i].isVisible = false;
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
	for (i = ver; i < v - 1; i++) addEdge(i, i + 1);
	addEdge(v - 1, ver);
	vector<int> arestas;
	for (i = 0; i < 9; i++) {
		arestas.clear();
		arestas.push_back(i);
		arestas.push_back(i + 11);
		arestas.push_back(i + 20);
		arestas.push_back(i + 10);
		addFace(arestas);
	}
	arestas.clear();
	arestas = { 9, 10, 29, 19 };
	addFace(arestas);
	arestas.clear();
	arestas = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
	addFace(arestas);
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
void Shape::printShape(SDL_Renderer* renderer, int width, int height, int mode) // Print usando faces
{
	int w, h, fa;
	int Xdi, Xdf, Ydi, Ydf;
	SDL_GetRendererOutputSize(renderer, &w, &h);
	//cout << "printShape" << endl;
	//cout << "e ae" << faces.size() << f << endl;
	//faces[11].isVisible = false;
	//faces[7].isVisible = false;

	
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
						ET.push_back(nEB);
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
				int y = ET[0].yMin;
				while (!ET.empty()) {
					if (!AL.empty()) {
						int i = 0;
						while (i < AL.size()) {
							if (AL[i].yMax == y) {
								AL.erase(AL.begin() + i);
							}
							else { i++; }
						}
						i = 0;
						while (i < ET.size()) {
							if (ET[i].yMax == y) {
								ET.erase(AL.begin() + i);
							}
							else { i++; }
						}
					}
					for (int i = 0; i < ET.size(); ++i) {
						if (ET[i].yMin == y)
						{
							AL.push_back(ET[i]);
							AL[-1].key = AL[-1].x;
						}
					}
					std::sort(AL.begin(), AL.end());
					// NingenAki
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
