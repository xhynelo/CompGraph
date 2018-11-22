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

Face::Face() : color(0), isVisible(true) {}

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
//*/

void Shape::setPosition(double x, double y)
{
	position.x = x;
	position.y = y;
}

void Shape::hider()
{
	cout << "entrei" << endl;
	int v0, v1, v2, p0x, p0y, p0z, p1x, p1y, p1z, px, py, pz, V;
	Vertex N, teste1, teste2;
	teste1 = Vertex(2, 9);
	teste2 = Vertex(2, 9);
	if (teste1 == teste2) { cout << "teste1 e teste2 sao iguai" << endl; }
	for (int i = 0; i < f; i++) {
		cout << "face: " << i << endl;
		for (int j = 0; j < faces[i].edges.size(); j++) {
			cout << "aresta " << j << ": " << faces[i].edges[j] << " -> ";
		}
		cout << endl;
		if (edges[faces[i].edges[0]].first == edges[faces[i].edges[1]].first) {
			v0 = edges[faces[i].edges[0]].second;
			cout << "face: " << i << " vertice: " << v0 << " v0.x " << vertices[v0].x << " v0.y: " << vertices[v0].y << " v0.y: " << vertices[v0].z << endl;
			v1 = edges[faces[i].edges[0]].first;
			cout << "face: " << i << " vertice: " << v1 << " v0.x " << vertices[v1].x << " v1.y: " << vertices[v1].y << " v1.y: " << vertices[v1].z << endl;
			v2 = edges[faces[i].edges[1]].second;
			cout << "face: " << i << " vertice: " << v2 << " v0.x " << vertices[v2].x << " v2.y: " << vertices[v2].y << " v2.y: " << vertices[v2].z << endl;
			cout << "face: " << i << " igual" << endl;
		}
		else {
			v0 = edges[faces[i].edges[0]].first;
			cout << "face: " << i << " vertice: " << v0 << " v0.x " << vertices[v0].x << " v0.y: " << vertices[v0].y << " v0.y: " << vertices[v0].z << endl;
			v1 = edges[faces[i].edges[0]].second;
			cout << "face: " << i << " vertice: " << v1 << " v0.x " << vertices[v1].x << " v1.y: " << vertices[v1].y << " v1.y: " << vertices[v1].z << endl;
			v2 = edges[faces[i].edges[1]].second;
			cout << "face: " << i << " vertice: " << v2 << " v0.x " << vertices[v2].x << " v2.y: " << vertices[v2].y << " v2.y: " << vertices[v2].z << endl;
			cout << "face: " << i << " vertice a first: " << edges[faces[i].edges[0]].first << " vertice a second: " << edges[faces[i].edges[0]].second << " vertice b first:: " << edges[faces[i].edges[1]].first << " vertice b second: " << edges[faces[i].edges[1]].second << endl;
			cout << "face: " << i << " diferente" << endl;
		}
		p0x = (vertices[v0].x + position.x) - (vertices[v1].x + position.x);
		p0y = (vertices[v0].y + position.y) - (vertices[v1].y + position.y);
		p0z = (vertices[v0].z + position.z) - (vertices[v1].z + position.z);
		p1x = (vertices[v2].x + position.x) - (vertices[v1].x + position.x);
		p1y = (vertices[v2].y + position.y) - (vertices[v1].y + position.y);
		p1z = (vertices[v2].z + position.z) - (vertices[v1].z + position.z);
		cout << "face: " << i << " vertice " << v0 << " " << v1 << " p0x: " << p0x << " p0y: " << p0y << " p0z: " << p0z << endl;
		cout << "face: " << i << " vertice " << v2 << " " << v1 << " p1x: " << p1x << " p1y: " << p1y << " p1z: " << p1z << endl;
		N.x = p0y * p1z - p0z * p1y;
		N.y = p0x * p1z - p0z * p1x;
		N.z = p0x * p1y - p0y * p1x;
		cout << "face: " << i << " Nx: " << N.x << " Ny: " << N.y << " Nz: " << N.z << endl;
		px = (vertices[v0].x + position.x) + 100;
		py = (vertices[v0].y + position.x) + 50;
		pz = (vertices[v0].z + position.x) + 100;
		cout << "face: " << i << " px: " << px << " py: " << py << " pz: " << pz << endl;
		V = px * N.x + py * N.y + pz * N.z;
		cout << "face: " << i << " Normal: " << V << endl;
		if (V >= 0) {
			faces[i].isVisible = true;
			cout << "face: " << i << " visivel" << endl;
		}
		else {
			faces[i].isVisible = false;
			cout << "face: " << i << " nao_visivel" << endl;
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
	arestas = { 29, 28, 27, 26, 25, 24, 23, 22, 21, 20 };
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
void Shape::printShape(SDL_Renderer* renderer, int width, int height) // Print usando faces
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
			for (auto edge = faces[fa].edges.begin(); edge != faces[fa].edges.end(); ++edge) {
				//cout << (*edge) << endl;
				Xdi = (int)(((this->x(edges[*edge].first) + position.x * position.scale) * width) / this->width);
				Ydi = (int)(((this->y(edges[*edge].first) + position.y * position.scale) * -height) / this->height + height);
				Xdf = (int)(((this->x(edges[*edge].second) + position.x * position.scale) * width) / this->width);
				Ydf = (int)(((this->y(edges[*edge].second) + position.y * position.scale) * -height) / this->height + height);
				//cout << Xdi << " " << Ydi << " " << Xdf << " " << Ydf << endl;
				SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE); // ... (renderer, r, g, b, alpha)
				SDL_RenderDrawLine(
					renderer,
					Xdi, Ydi,
					Xdf, Ydf
				);
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
