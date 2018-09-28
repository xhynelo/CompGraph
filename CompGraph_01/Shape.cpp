#include "stdafx.h"
#include "Shape.h"

#define XuMin 0
#define YuMin 0
#define XuMax 100
#define YuMax 100

Vertex::Vertex(int mX, int mY) : x(mX), y(mY), z(0), scale(1) {}

Vertex::Vertex() {}

int Vertex::operator[](int i) {
	switch (i)
	{
	case 0: return x;
	case 1: return y;
	case 2: return z;
	case 3: return scale;
	default:
		return -1;
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

int Shape::x(int n)
{
	if (v <= n) return -1;
	return vertices[n].x * vertices[n].scale;
}

int Shape::y(int n)
{
	if (n >= v) return -1;
	return vertices[n].y * vertices[n].scale;
}

int Shape::getV()
{
	return v;
}

void Shape::scale(int n)
{
	for (int i = 0; i < v; i++) {
		vertices[i].scale = n;
	}
}

void Shape::addVertex(int x, int y)
{
	if (vXmin < 0 || x < this->x(vXmin)) vXmin = v;
	if (vYmin < 0 || y < this->y(vYmin)) vYmin = v;
	if (vXMAX < 0 || x > this->x(vXMAX)) vXMAX = v;
	if (vYMAX < 0 || y > this->y(vYMAX)) vYMAX = v;
	addVertex(Vertex(x, y));
}

void Shape::addVertex(Vertex v)
{
	vertices.push_back(v);
	this->v++;
}


void Shape::move(int x, int y)
{
	for (int i = 0; i < v; i++) {
		vertices[i].x += x;
		vertices[i].y += y;
	}
}

void Shape::printShape(HDC hdc)
{
	/*MoveToEx(hdc, this->x(0), this->y(0), NULL);	
	for (int i = 1; i < this->getV(); i++) {
		LineTo(hdc, this->x(i), this->y(i));
	}
	LineTo(hdc, this->x(0), this->y(0));*/
	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		MoveToEx(hdc, this->x(edge->first), this->y(edge->first), NULL);
		LineTo(hdc, this->x(edge->second), this->y(edge->second));
	}
}

void Shape::metaShape()
{
	cout << "oi, como vai?" << endl;
	printf("Hello, world\n");
}

void Shape::addEdge(int v1, int v2)
{
	if (v <= v1 || v <= v2) return;
	edges.push_back(pair<int, int>(v1, v2));
	e++;
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
			cout << vertex << endl;
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

void Shape::remaping(HWND hWnd, HDC hdc)
{
	int Xumin = 0, Yumin = 0, XuMAX = 100, YuMAX = 100;
	//int XDmin = this->x(vXmin), YDmin = this->y(vYmin), XDMAX = this->x(vXMAX), YDMAX = this->y(vYMAX);
	int XDmin = 0, YDmin = 0, XDMAX, YDMAX;
	int Xu, Yu, scalexy, newScale;
	int width, height;
	RECT rect;
	if (GetWindowRect(hWnd, &rect))
	{
		width = rect.right - rect.left;
		height = rect.bottom - rect.top;
	}
	//*
	XDMAX = width;
	YDMAX = height;
	//*/
	/*
	(3,1)->(9,9)
	XD = [(Xu - Xumin)(XDMAX - XDmin) / (XuMAX - Xumin)] + XDmin)
	YD = [(Yu - Yumin)(YDMAX - YDmin) / (YuMAX - Yumin)] + YDmin)
	*/
	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		//scalexy = vertices[0].scale;
		//newScale = (scalexy * 100) / XDMAX;
		//scale(newScale);
		Xu = this->x(edge->first);
		Yu = this->y(edge->first);
		MoveToEx(hdc, (Xu - Xumin)*(XDMAX - XDmin) / (XuMAX - Xumin) + XDmin, (Yu - Yumin)*(YDMAX - YDmin) / (YuMAX - Yumin) + YDmin, NULL);
		Xu = this->x(edge->second);
		Yu = this->y(edge->second);
		LineTo(hdc, (Xu - Xumin)*(XDMAX - XDmin) / (XuMAX - Xumin) + XDmin, (Yu - Yumin)*(YDMAX - YDmin) / (YuMAX - Yumin) + YDmin);
	}
}

void Shape::matrixMult(vector<vector<float>> mat)
{
	vector<vector<float>> matrixAux = this->matrix;
	float sum;
	int i, j, k;
	for (i=0; i < 4; i++) {
		for (j=0; j < 4; j++) {
			sum = 0;
			for (k = 0; k < 4; k++) {
				sum += matrixAux[j][k] * mat[k][j];
			}
			matrix[i][j] = sum;
		}
	}
}
