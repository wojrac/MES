#pragma once
#include"Node.h"
#include"Element.h"

struct Grid
{
	int nHG , nLG;
	//Grid(int, int);
	Grid();
	Node** nodes;
	Grid(int, int,int);
	Element* elements;

	~Grid()
	{
		delete[] nodes;
		delete[] elements;
	}
	//Element* elements = new Element[100000];
	
};


