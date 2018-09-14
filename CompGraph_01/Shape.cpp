#include "stdafx.h"
#include "Shape.h"


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