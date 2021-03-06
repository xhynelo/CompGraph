// CompGraph_01.cpp : Defines the entry point for the application.
//
#include "stdafx.h"
#include "CompGraph_01.h"
#include <iostream>
#include "Shape.h"
#include <fstream>
#include <cmath>
#include <ctime>

#include <SDL.h>

#define WIDTH 1000
#define HEIGHT 1000

using namespace std;

int main(int argc, char* argv[])
{
	Shape u;
	int i = 0;
	vector<vector<float>> translacao;

	u.setSRU(100, 100);
	
	
	
	u.addVertex(1, 8);
	u.addVertex(3, 8);
	u.addVertex(3, 3);
	u.addVertex(5, 3);
	u.addVertex(5, 8);
	u.addVertex(7, 8);
	u.addVertex(7, 1);
	//u.addVertex(6, 0);
	//u.addVertex(2, 0);
	u.addVertex(1, 1);
	u.setPosition(100,50);
	for (int i = 0; i < 7; i++) {
		if (i == 6) {
			u.addEdgeCurva(i + 1, i, Vertex(1, -1, 0), Vertex(7, -1, 0));
		}
		else {
			u.addEdge(i + 1, i);
		}
	}
	//u.addEdge(7, 6);
	u.addEdge(7, 0);
	//u.addEdge(9, 8);
	//u.addEdge(9, 0);

	/*
	u.addEdge(2, 0); // 10
	u.addEdge(9, 2); // 11
	u.addEdge(2, 8); // 12
	u.addEdge(8, 3); // 13
	u.addEdge(3, 7); // 14
	u.addEdge(6, 3); // 15
	u.addEdge(3, 5); // 16
	*/

	u.scale(1);
	vector<int> edges;
	//edges = { 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	edges = { 7, 6, 5, 4, 3, 2, 1, 0 };
	//edges = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	/*
	edges = { 10, 1, 0 };
	//edges = { 0, 1, 10 };
	u.addFace(edges, 900);

	edges = { 11, 10, 9 };
	//edges = { 9, 10, 11 };
	u.addFace(edges, 900);

	edges = { 8, 12, 11 };
	//edges = { 11, 12, 8 };
	u.addFace(edges, 900);

	edges = { 13, 2, 12 };
	//edges = { 12, 2, 13 };
	u.addFace(edges, 900);

	edges = { 7, 14,13 };
	//edges = { 13, 14, 7 };
	u.addFace(edges, 900);

	edges = { 6, 15, 14 };
	//edges = { 14, 15, 6 };
	u.addFace(edges, 900);

	edges = { 5, 16, 15 };
	//edges = { 15, 16, 5 };
	u.addFace(edges, 900);

	edges = { 16, 4, 3 };
	//edges = { 3, 4, 16 };
	*/
	u.addFace(edges, 900);
	//u.metaShape();
	//u.readShape("testeGraph.txt");
	//u.move(3, 1);
	
	
	//*
	translacao = {
			{1.0, 0.0, 0.0, 0.0},
			{0.0, 1.0, 0.0, 0.0},
			{0.0, 0.0, 1.0, 0.0},
			{2.0, 2.0, 0.0, 1.0}
	};

	int m;
	//u.matrixMult(translacao);
	//*/

	/*
	cout << width << "  " << height << endl;
	cout << "-------&&&-------" << endl;

	while (i < 20) {


		u.matrixMult(translacao);
		u.tranform();
		i++;
		cout << u.x(0) << "  " << u.y(0) << endl;
		//cout << i << endl;

		cout << "-----__-----" << endl;
		ShowWindow(hWnd, nCmdShow);
		UpdateWindow(hWnd);
		mySleep(2);
	}*/
	clock_t msec = 100;
	
	if (SDL_Init(SDL_INIT_VIDEO) == 0) {
		SDL_Window* window = NULL;
		SDL_Renderer* renderer = NULL;

		if (SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, 0, &window, &renderer) == 0) {
			SDL_bool done = SDL_FALSE;
			//cout << "antes" << endl;
			//u.metaShape(renderer, WIDTH, HEIGHT);
			clock_t start_time = clock();
			//u.remaping(.0, 100.0);
			//u.setPosition(50, 50);
			//u.remaping(WIDTH, HEIGHT);
			//u.DrawFilledPolygon(u.faces[i], renderer);
			m = u.minVertex();
			//cout << m << endl;
			u.slide(m);
			//u.metaShape();

			//u.rotate(0, -30, 0);
			//u.hider(60, 60, -100, 0);
			//u.hider(200, 200, -100, 1);
			while (!done) {
				SDL_Event event;

				SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
				SDL_RenderClear(renderer);
				
				clock_t time = clock();

				u.reseta();
				//pair<double, Vertex> quartenion = u.multQuaternion(180, Vertex(1, 0, 0), 90, Vertex(0, 0, 1));
				//cout << "angulo " << quartenion.first << " vetor " << quartenion.second.x << " " << quartenion.second.y << " " << quartenion.second.z << endl;
				//u.rotateQ1(quartenion.first, quartenion.second.x, quartenion.second.y, quartenion.second.z);
				u.rotateQ((i + 1) * 3, 0, 0, 1);

				//*
				if ((i < 30) && (time - start_time >= msec)) {
					u.translade(-(100 / 30), (100 / 60));

					//u.translade(-(100/30.0), (100 /60.0));
					//u.rotate((i+1)*3);
					//u.rotateQ(45, 0, 0, 1);
					

					//u.rotate(pi/2);
					//cout << "antes" << endl;

					
					i++;
					start_time = start_time + msec;
				}
				/*/
				//u.metaShape(renderer, WIDTH, HEIGHT);
				SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE); // ... (renderer, r, g, b, alpha)
				SDL_RenderDrawLine(
					renderer,
					70, 10,
					70, 70
				);
				/*
				SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE); // ... (renderer, r, g, b, alpha)
				SDL_RenderDrawLine(
					renderer,
					70, 10,
					70, 70
				);
				SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE); // ... (renderer, r, g, b, alpha)
				SDL_RenderDrawLine(
					renderer,
					70, 10,
					70, 70
				);
				*/

				//u.hider(30, 60, -100, 0);

				u.projection(120);

				//u.printShape(renderer, WIDTH, HEIGHT, SOLID_WITH_LIGHT);
				u.printShape(renderer, WIDTH, HEIGHT, WIRE_FRAME);

				SDL_RenderPresent(renderer);

				while (SDL_PollEvent(&event)) {
					if (event.type == SDL_QUIT) {
						done = SDL_TRUE;
					}
				}
			}
		}

		if (renderer) {
			SDL_DestroyRenderer(renderer);
		}
		if (window) {
			SDL_DestroyWindow(window);
		}
	}
	SDL_Quit();
	return 0;
}