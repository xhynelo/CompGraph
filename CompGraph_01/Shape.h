#pragma once
#include <windef.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <SDL.h>

using namespace std;


struct Face
{
	double color;
	vector<int> edges;
	bool isVisible;
	
	Face();
};

struct Vertex
{
	double x = 0;
	double y = 0;
	double z = 0;
	double scale = 1;

	

	Vertex(double mX, double mY);
	Vertex();
	double& operator[](int i);
	friend std::istream& operator>>(std::istream &is, Vertex &vertex);
	friend std::ostream& operator<<(std::ostream &is, Vertex &vertex);
	bool operator ==(const Vertex &b) const;
};

class Shape
{
	vector<Face> faces;
	vector<Vertex> vertices;
	vector<pair<int, int>> edges;
	//number of elements on the above arrays
	int v = 0, e = 0, f = 0;
	//int XDmin = -1, YDmin = -1, XDMAX = 0, YDMAX = 0;
	int vXmin = -1, vYmin = -1, vXMAX = -1, vYMAX = -1;

	Vertex position;
	int pivot = 0;
	double height, width;
	int invertX = 0, invertY = 0;

	public:
		Shape();

		double x(int n);
		double y(int n);
		double z(int n);
		double getV();
		vector<vector<double>> matrix;
		
		void addVertex(double x, double y);
		void addVertex(Vertex v);
		void addEdge(int v1, int v2);
		void addFace(vector<int>);
		void setSRU(double x, double y);
		void scale(double n);
		void tranform();
		void translade(double x, double y);
		void rotate(double theta);
		void rotate(double thetaX, double thetaY, double thetaZ);
		void setPosition(double x, double y);
		void hider();
		void slide(double tam);
		void projection(double theta);
		void printShape(SDL_Renderer* renderer, int width, int height);
		void metaShape(SDL_Renderer* renderer, int width, int height);
		void readShape(string name);

		vector<vector<double>> matrixMult(vector<vector<double>> mat1, vector<vector<double>> mat2);
};