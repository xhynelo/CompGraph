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
	u.addVertex(6, 0);
	u.addVertex(2, 0);
	u.addVertex(1, 1);
	u.setPosition(50,50);
	for (int i = 0; i < 9; i++) u.addEdge(i, i + 1);
	u.addEdge(9, 0);
	u.scale(1);
	vector<int> arestas;
	for (int i = 0; i < 10; i++) {
		arestas.push_back(i);
	}
	u.addFace(arestas);

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
			u.slide(0.2);
			while (!done) {
				SDL_Event event;

				SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
				 SDL_RenderClear(renderer);
				
				clock_t time = clock();
				

				//*
				if ((i < 0) && (time - start_time >= msec)) {

					u.translade(-(100/30.0), (100 /60.0));
					u.rotate(3);
					//u.rotate(pi/2);
					//cout << "antes" << endl;

					
					i++;
					start_time = start_time + msec;
				}
				//*/
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

				u.projection(120);
				u.printShape(renderer, WIDTH, HEIGHT);
				
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