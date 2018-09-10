#pragma once
class Shape
{
	int *vertices = NULL, *edges = NULL, *faces = NULL;
	//number of elements on the above arrays
	int v = 0, e = 0, f = 0;

	public:
		Shape();

		int x(int n);
		int y(int n);
		int getV();

		void scale(int n);
		void addVertex(int x, int y);
		void move(int x, int y);
};