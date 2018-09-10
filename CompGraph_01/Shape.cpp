#include "stdafx.h"
#include "Shape.h"


Shape::Shape()
{
	v = 0;
	e = 0;
	f = 0;
	delete[] vertices;
	vertices = NULL;
}

int Shape::x(int n)
{
	if (v <= n) return -1;
	return vertices[n * 4] * vertices[n * 4 + 3];
}

int Shape::y(int n)
{
	if (n >= v) return -1;
	return vertices[n * 4 + 1] * vertices[n * 4 + 3];
}

int Shape::getV()
{
	return v;
}

void Shape::scale(int n)
{
	for (int i = 0; i < v; i++) {
		vertices[i * 4 + 3] = n;
	}
}

void Shape::addVertex(int x, int y)
{
	if (!vertices) {
		vertices = new int[4]{ x, y, 0, 1 };
		v = 1;
	}
	else {
		int *aux = vertices;
		vertices = new int[4 * (v + 1)];
		for (int i = 0; i < v * 4; i++) {
			vertices[i] = aux[i];
		}
		delete[] aux;
		vertices[v * 4] = x;
		vertices[v * 4 + 1] = y;
		vertices[v * 4 + 2] = 0;
		vertices[v * 4 + 3] = 1;
		v++;
	}
}

void Shape::move(int x, int y)
{
	for (int i = 0; i < v; i++) {
		vertices[i * 4] += x;
		vertices[i * 4 + 1] += y;
	}
}
