#include"Grid.h"

/*Grid::Grid(Node* arraynodes, int noN)
{
	this->numberOfNodes = noN;
	arraynodes = new Node[noN];
}
Grid::Grid(Element* arrayelements, int noE)
{
	this->numberOfElements = noE;
	arrayelements = new Element[noE];
}*/
/*Grid::Grid(int noN, int noE)
{
	this->numberOfElements = noE;
	this->numberOfNodes = noN;
}*/
Grid::Grid()
{

}
Grid::Grid(int n, int m , int ne )
{
	this->nLG = n;
	this->nHG = m;
	nodes = new Node*[n];
	for (int i = 0; i < n; i++)
	{
		nodes[i] = new Node[m];
	}
	elements = new Element[ne];
}