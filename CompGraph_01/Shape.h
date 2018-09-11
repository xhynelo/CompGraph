#pragma once
#include <windef.h>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

struct Vertex
{
	int x = 0;
	int y = 0;
	int z = 0;
	int scale = 1;

	Vertex(int mX, int mY);
	Vertex();
	int operator[](int i);
	friend std::istream& operator>>(std::istream &is, Vertex &vertex);
	friend std::ostream& operator<<(std::ostream &is, Vertex &vertex);
};

class Shape
{
	int *faces = NULL;
	vector<Vertex> vertices;
	vector<pair<int, int>> edges;
	//number of elements on the above arrays
	int v = 0, e = 0, f = 0;

	public:
		Shape();

		int x(int n);
		int y(int n);
		int getV();

		void scale(int n);
		void addVertex(int x, int y);
		void addVertex(Vertex v);
		void move(int x, int y);
		void printShape(HDC hdc);
		void addEdge(int v1, int v2);
		void readShape(string name);
};