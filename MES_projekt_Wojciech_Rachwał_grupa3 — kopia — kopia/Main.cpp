#include"Node.h"
#include"Element.h"
#include"Grid.h"
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<string>
#include"functions.h"



int main()
{
	fstream plik;
	plik.open("Properties.txt" , ios::in);
	if (plik.good() == false)
	{
		cout << "plik nie istnieje " << endl;
		exit(0);
	}
	int ne = 0;
	int nodesId = 0;
	double H;
	double L;
	int nH ;
	int nL ;
	double tab1[4][4];                          //z zeszytu jest to tabela dN[i] /d ksi
	double tab2[4][4];							//z zeszytu jest to tabela dN[i] /d eta
	double funkcjeKsztaltu[4][4];	
	double funkcjeKsztaltu_dS[8][4];						//tabelka podobna do powyzszych ale bez pochodnych
	double przedjakobianDlaPierwszegoPktCalk[2][2];    //jakobian w tym wypadu rozumiem jako ta macierz pierwotna przed odwroceniem
	double przedjakobianDlaDrugiegoPktCalk[2][2];
	double przedjakobianDlaTrzeciegoPktCalk[2][2];
	double przedjakobianDlaCzwartegoPktCalk[2][2];
	double det1punkt = 0;
	double det2punkt = 0; 
	double det3punkt = 0;
	double det4punkt = 0;
	double Jakobian1pkt[2][2];
	double Jakobian2pkt[2][2];
	double Jakobian3pkt[2][2];
	double Jakobian4pkt[2][2];
	double TabDNPoDx[4][4];
	double TabDNPoDy[4][4];
	double bokA;   //1 i 3 powierzchnia
	double bokB;	// 2 i 4 powierzchnia
	double tab1dla1pktdnpodx[4];
	double tab2dla2pktdnpodx[4];
	double tab3dla3pktdnpodx[4];
	double tab4dla4pktdnpodx[4];
	double tab5dla1pktdnpody[4];
	double tab6dla2pktdnpody[4];
	double tab7dla3pktdnpody[4];
	double tab8dla4pktdnpody[4];
	double tab2DpoTraspozycji1pktpodx[4][4];
	double tab2DpoTraspozycji2pktpodx[4][4];
	double tab2DpoTraspozycji3pktpodx[4][4];
	double tab2DpoTraspozycji4pktpodx[4][4];
	double tab2DpoTraspozycji1pktpody[4][4];
	double tab2DpoTraspozycji2pktpody[4][4];
	double tab2DpoTraspozycji3pktpody[4][4];
	double tab2DpoTraspozycji4pktpody[4][4];
	double tabPoDodaniuDyPlusDxrazkpkt1[4][4];
	double tabPoDodaniuDyPlusDxrazkpkt2[4][4];
	double tabPoDodaniuDyPlusDxrazkpkt3[4][4];
	double tabPoDodaniuDyPlusDxrazkpkt4[4][4];
	double Ndla1punktudoC[4];
	double Ndla2punktudoC[4]; 
	double Ndla3punktudoC[4];
	double Ndla4punktudoC[4];
	double NpoTranspondla1punktudoC[4][4];
	double NpoTranspondla2punktudoC[4][4];
	double NpoTranspondla3punktudoC[4][4];
	double NpoTranspondla4punktudoC[4][4];
	double deltaTau = 1;
	double initialTemperature = 100;
	double simulation_time =20;
	int number_iteration = (int)(simulation_time / deltaTau);
	Node  pkt_calkowania_dS[8];
	double wspolrzedna_calkowania = 1 / (sqrt(3));
	string linia;
	int line_counter = 1;
	string::size_type st;
	while (getline(plik, linia))
	{
		switch(line_counter)
		{
			case 1:
				H = stod(linia,&st); break;
			case 2:
				L = stod(linia, &st); break;
			case 3:
				nH = atoi(linia.c_str()); break;
			case 4:
				nL = atoi(linia.c_str()); break;
		}
		line_counter++;
	}
	/*cout << "H = " << H << endl;
	cout << "L = " << L << endl;
	cout << "nH = " << nH << endl;
	cout << "nL = " << nL << endl;*/
	double** MatrixHGlobal = new double*[nL*nH];
	for (int i = 0; i < nL*nH; i++)
	{
		MatrixHGlobal[i] = new double[nL*nH];
	}
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0 ;j < (nH * nL); j++)
		{
			MatrixHGlobal[i][j] = 0;
		}
	}
	double** MatrixCGlobal = new double*[nL*nH];
	for (int i = 0; i < nL*nH; i++)
	{
		MatrixCGlobal[i] = new double[nL*nH];
	}
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH * nL); j++)
		{
			MatrixCGlobal[i][j] = 0;
		}
	}
	double** MatrixCGlobalprzezdTau = new double*[nL*nH];
	for (int i = 0; i < nL*nH; i++)
	{
		MatrixCGlobalprzezdTau[i] = new double[nL*nH];
	}
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH * nL); j++)
		{
			MatrixCGlobalprzezdTau[i][j] = 0;
		}
	}
	double** MatrixHzDaszkiem = new double*[nL*nH];
	for (int i = 0; i < nL*nH; i++)
	{
		MatrixHzDaszkiem[i] = new double[nL*nH];
	}
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH * nL); j++)
		{
			MatrixHzDaszkiem[i][j] = 0;
		}
	}
	double * VectorPGlobal = new double[nL*nH];
	for (int i = 0; i < nH * nL; i++)
	{
		VectorPGlobal[i]  =0;
	}
	double * VectorPzDaszkiem = new double[nL*nH];
	for (int i = 0; i < nH * nL; i++)
	{
		VectorPzDaszkiem[i] = 0;
	}
	int* nodesIds = new int[nH*nL];
	double jumpX = (double)L / (double)(nL-1);
	double jumpY = (double)H / (double)(nH-1);
	Node* arrayNodesMain = new Node[nH*nL];
	for (int i = 0; i < nH*nL; i++)
	{
		arrayNodesMain[i].x = 0;
		arrayNodesMain[i].y = 0;
	}
	bokA = jumpX;
	bokB = jumpY;
	ne = (nH - 1)*(nL - 1); 
	Grid grid(nL, nH, ne);
	Element* element = new Element[ne];
	for (int i = 0; i < ne; i++) //do kazdego elementu musze przypisac wezly 
	{
		// jest typ element  grid.elements[i] = element.nodes_ids; 4 idiki wezlow dla kazdego elementu po prawej stronie musze miec typ element

	}
	for (int i = 0; i < nL; i++)
	{
		for (int j = 0; j < nH; j++)
		{
			nodesId++;
			grid.nodes[i][j].x = arrayNodesMain[j].x + jumpX *i;
			grid.nodes[i][j].y = arrayNodesMain[j].y + jumpY *j;
			nodesIds[j + nH*i] = nodesId;
		}
	}
	for (int i = 0; i < nL; i++)
	{
		for (int j = 0; j < nH; j++)
		{
			arrayNodesMain[j + nH*i] = grid.nodes[i][j];
		}
	}
	for (int i = 0; i < nL - 1; i++)
	{
		for (int j = 0; j < nH - 1; j++)
		{

			for (int k = 1; k <= 4; k++)
			{
				int temp = i + j + (nH - 1)*i;
				if (k == 1)
				{
					temp = i + j + (nH - 1)*i;
				}
				if (k == 2)
				{
					temp = temp + nH;
				}
				if (k == 3)
				{
					temp = temp + nH + 1;
				}
				if (k == 4)
				{
					temp++;
				}
				element[j + (nH - 1)*i].nodes_ids[k - 1] = nodesIds[temp];
			}
			grid.elements[j + (nH - 1)*i] = element[j + (nH - 1)*i];
		}
	}
	/*for (int i = 0; i < nL - 1; i++)    // no i TU TEZ
	{
		for (int j = 0; j < nH - 1; j++)
		{
			cout << "pierwszy wezel elementu " << j + (nH - 1)*i + 1 << ": " << grid.elements[j + (nH - 1)*i].nodes_ids[0] << endl;
			cout << "drugi wezel elementu " << j + (nH - 1)*i + 1 << ": " << grid.elements[j + (nH - 1)*i].nodes_ids[1] << endl;
			cout << "trzeci wezel elementu " << j + (nH - 1)*i + 1 << ": " << grid.elements[j + (nH - 1)*i].nodes_ids[2] << endl;
			cout << "czwarty wezel elementu " << j + (nH - 1)*i + 1 << ": " << grid.elements[j + (nH - 1)*i].nodes_ids[3] << endl;
			

			cout << endl;
		}
		cout << endl;
		cout << endl;
	}  */                                                        //to do odkomentowania
	/*for (int j = 0; j < ne; j++)    //ladnie , czysciej i przejrzysciej od tego na gorze
	{
		cout << "element nr: " << j+1 << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << grid.elements[j].nodes_ids[i] <<" "<< endl;
		}
		cout << endl;
	}*
	/*for (int i = 0; i < nH*nL; i++)
	{
		cout << nodesIds[i] << endl;
	}
	
	//arrayNodesMain[grid.elements[0].nodes_ids[0]].x 
	for (int i = 0; i < nL; i++)
	{
		for (int j = 0; j < nH; j++)
		{
			cout << "x" << i << ": " << grid.nodes[i][j].x << endl;
			cout << "y" << j << ": " << grid.nodes[i][j].y << endl;
			cout << "id dla tego wezla: " << nodesIds[j + nH *i] << endl;
		}
	}		*/																						//2D->1D
	/*for (int i = 0; i < nH*nL; i++)									//TU TEz
	{
		cout << "x: " << arrayNodesMain[i].x << endl;
		cout << "y: " << arrayNodesMain[i].y << endl;
		cout << "id dla tego wezla: " << nodesIds[i] << endl;
	}*/
	//cout << "Skrajne elementy mojej siatki to maja wezly i wspolrzedne przedstawione ponizej: " << endl;
	/*for (int i = 0; i < 4; i++) 
	{
		cout<<arrayNodesMain[grid.elements[0].nodes_ids[i]-1].x << endl;
		cout << arrayNodesMain[grid.elements[0].nodes_ids[i]-1].y << endl;

	}*/
	//przypisanie wezlom
	pkt_calkowania_dS[0].ksi = (-1)*wspolrzedna_calkowania;   //pow 1 pkt
	pkt_calkowania_dS[0].eta = -1;			
	pkt_calkowania_dS[1].ksi = wspolrzedna_calkowania;//pow pkt 2
	pkt_calkowania_dS[1].eta = -1;
	pkt_calkowania_dS[2].ksi = 1;			//pow 2
	pkt_calkowania_dS[2].eta = (-1)*wspolrzedna_calkowania;
	pkt_calkowania_dS[3].ksi = 1;
	pkt_calkowania_dS[3].eta = wspolrzedna_calkowania;
	pkt_calkowania_dS[4].ksi = wspolrzedna_calkowania;    // pow 3
	pkt_calkowania_dS[4].eta = 1;
	pkt_calkowania_dS[5].ksi = (-1)* wspolrzedna_calkowania;
	pkt_calkowania_dS[5].eta = 1;
	pkt_calkowania_dS[6].ksi = -1;				//pow 4
	pkt_calkowania_dS[6].eta = wspolrzedna_calkowania;
	pkt_calkowania_dS[7].ksi = -1;
	pkt_calkowania_dS[7].eta = (-1)* wspolrzedna_calkowania;
	for (int i = 0; i < 8; i++)
	{
		funkcjeKsztaltu_dS[i][0] = 0.25*(1 - pkt_calkowania_dS[i].ksi)*(1 - pkt_calkowania_dS[i].eta);  //dla kazdego pkt calkowania liczymy fk
		funkcjeKsztaltu_dS[i][1] = 0.25*(1 + pkt_calkowania_dS[i].ksi)*(1 - pkt_calkowania_dS[i].eta);
		funkcjeKsztaltu_dS[i][2] = 0.25*(1 + pkt_calkowania_dS[i].ksi)*(1 + pkt_calkowania_dS[i].eta);
		funkcjeKsztaltu_dS[i][3] = 0.25*(1 - pkt_calkowania_dS[i].ksi)*(1 + pkt_calkowania_dS[i].eta);
	}
	//funkcje ksztaltu:
	/*for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << funkcjeKsztaltu_dS[i][j] << endl;
		}
		cout << endl;
	}*/
	//cout << "o wspolrzednych: " < endl;
	/*for (int i = 0; i < 4; i++)
	{
		cout << "dla " << grid.elements[0].nodes_ids[i] << " wezla " << endl;
		cout << "wspolrzedne wynosza " << endl;
		cout << arrayNodesMain[grid.elements[0].nodes_ids[i] - 1].x << endl;
		cout << arrayNodesMain[grid.elements[0].nodes_ids[i] - 1].y << endl;
	}
	for (int i = 0; i < 4; i++)
	{
		cout << "dla " << grid.elements[ne-1].nodes_ids[i] << " wezla " << endl;
		cout << "wspolrzedne wynosza " << endl;
		cout << arrayNodesMain[grid.elements[ne-1].nodes_ids[i] - 1].x << endl;
		cout << arrayNodesMain[grid.elements[ne-1].nodes_ids[i] - 1].y << endl;
	}*/
	/////////////////////////////////////////***************************

	for (int e = 0; e < ne; e++)
	{
		grid.elements[e].ambient_temperature = 1200;
		grid.elements[e].alfa = 300;
		grid.elements[e].specific_heat = 700;
		grid.elements[e].density = 7800;
		grid.elements[e].k = 25;										//przypisuje kazdemu elementowi conductivity
		grid.elements[e].pkt_calkowania[0].ksi = (-1)* wspolrzedna_calkowania;  //kazdy element ma 4 pkt calkowania (takie same)
		grid.elements[e].pkt_calkowania[0].eta = (-1)* wspolrzedna_calkowania;	// tu sa one w ukladzie lokalnym
		grid.elements[e].pkt_calkowania[1].ksi = wspolrzedna_calkowania;
		grid.elements[e].pkt_calkowania[1].eta = (-1)* wspolrzedna_calkowania;
		grid.elements[e].pkt_calkowania[2].ksi = wspolrzedna_calkowania;
		grid.elements[e].pkt_calkowania[2].eta = wspolrzedna_calkowania;
		grid.elements[e].pkt_calkowania[3].ksi = (-1)* wspolrzedna_calkowania;
		grid.elements[e].pkt_calkowania[3].eta = wspolrzedna_calkowania;

		grid.elements[e].pkt_calkowania[0].x = arrayNodesMain[grid.elements[e].nodes_ids[0] - 1].x;  //tu sa w ukladzie globalnym
		grid.elements[e].pkt_calkowania[0].y = arrayNodesMain[grid.elements[e].nodes_ids[0] - 1].y;
		grid.elements[e].pkt_calkowania[1].x = arrayNodesMain[grid.elements[e].nodes_ids[1] - 1].x;
		grid.elements[e].pkt_calkowania[1].y = arrayNodesMain[grid.elements[e].nodes_ids[1] - 1].y;
		grid.elements[e].pkt_calkowania[2].x = arrayNodesMain[grid.elements[e].nodes_ids[2] - 1].x;
		grid.elements[e].pkt_calkowania[2].y = arrayNodesMain[grid.elements[e].nodes_ids[2] - 1].y;
		grid.elements[e].pkt_calkowania[3].x = arrayNodesMain[grid.elements[e].nodes_ids[3] - 1].x;
		grid.elements[e].pkt_calkowania[3].y = arrayNodesMain[grid.elements[e].nodes_ids[3] - 1].y;
		/*for (int i = 0; i < 4; i++)
		{
			cout << "ksi " << grid.elements[e].pkt_calkowania[i].ksi << endl;
			cout << "eta " << grid.elements[e].pkt_calkowania[i].eta << endl;
		}*/
		for (int i = 0; i < 4; i++)
		{
			tab1[i][0] = 0.25*(grid.elements[e].pkt_calkowania[i].eta - 1);			//z1  dla i =0  
			tab1[i][1] = 0.25*(1 - grid.elements[e].pkt_calkowania[i].eta);			//z2 dla i =0 
			tab1[i][2] = 0.25*(1 + grid.elements[e].pkt_calkowania[i].eta);			//z3 dla i =0
			tab1[i][3] = 0.25*((-1)* (grid.elements[e].pkt_calkowania[i].eta) - 1);			//z4 dla i =0
		}
		for (int i = 0; i < 4; i++)
		{
			tab2[i][0] = 0.25*(grid.elements[e].pkt_calkowania[i].ksi - 1);
			tab2[i][1] = 0.25*((-1) *(grid.elements[e].pkt_calkowania[i].ksi) - 1);
			tab2[i][2] = 0.25*(grid.elements[e].pkt_calkowania[i].ksi + 1);
			tab2[i][3] = 0.25*((-1) *(grid.elements[e].pkt_calkowania[i].ksi) + 1);

		}
		for (int i = 0; i< 4; i++)
		{
			funkcjeKsztaltu[i][0] = 0.25*(1 - grid.elements[e].pkt_calkowania[i].ksi)*(1 - grid.elements[e].pkt_calkowania[i].eta);//dla kazdego pc liczymy fk
			funkcjeKsztaltu[i][1] = 0.25*(1 + grid.elements[e].pkt_calkowania[i].ksi)*(1 - grid.elements[e].pkt_calkowania[i].eta);
			funkcjeKsztaltu[i][2] = 0.25*(1 + grid.elements[e].pkt_calkowania[i].ksi)*(1 + grid.elements[e].pkt_calkowania[i].eta);
			funkcjeKsztaltu[i][3] = 0.25*(1 - grid.elements[e].pkt_calkowania[i].ksi)*(1 + grid.elements[e].pkt_calkowania[i].eta);
		}
		/*cout << "TAB1" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << "w " << i + 1 << " punkcie calkowania , dla " << j + 1 << " funkcji kszaltu  :" << tab1[i][j] << endl;
			}
		}
		cout << "TAB2" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << "w " << i + 1 << " punkcie calkowania , dla  " << j + 1 << " funkcji ksztaltu  :" << tab2[i][j] << endl;
			}
		}*/
		// pierwszy punkt calkowania
		przedjakobianDlaPierwszegoPktCalk[0][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaPierwszegoPktCalk[0][0] = przedjakobianDlaPierwszegoPktCalk[0][0] + tab1[0][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaPierwszegoPktCalk[0][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaPierwszegoPktCalk[0][1] = przedjakobianDlaPierwszegoPktCalk[0][1] + tab2[0][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaPierwszegoPktCalk[1][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaPierwszegoPktCalk[1][0] = przedjakobianDlaPierwszegoPktCalk[1][0] + tab1[0][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		przedjakobianDlaPierwszegoPktCalk[1][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaPierwszegoPktCalk[1][1] = przedjakobianDlaPierwszegoPktCalk[1][1] + tab2[0][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		//drugi pkt calkowania
		przedjakobianDlaDrugiegoPktCalk[0][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaDrugiegoPktCalk[0][0] = przedjakobianDlaDrugiegoPktCalk[0][0] + tab1[1][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaDrugiegoPktCalk[0][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaDrugiegoPktCalk[0][1] = przedjakobianDlaDrugiegoPktCalk[0][1] + tab2[1][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaDrugiegoPktCalk[1][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaDrugiegoPktCalk[1][0] = przedjakobianDlaDrugiegoPktCalk[1][0] + tab1[1][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		przedjakobianDlaDrugiegoPktCalk[1][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaDrugiegoPktCalk[1][1] = przedjakobianDlaDrugiegoPktCalk[1][1] + tab2[1][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		//trzeci pkt calkowania
		przedjakobianDlaTrzeciegoPktCalk[0][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaTrzeciegoPktCalk[0][0] = przedjakobianDlaTrzeciegoPktCalk[0][0] + tab1[2][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaTrzeciegoPktCalk[0][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaTrzeciegoPktCalk[0][1] = przedjakobianDlaTrzeciegoPktCalk[0][1] + tab2[2][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaTrzeciegoPktCalk[1][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaTrzeciegoPktCalk[1][0] = przedjakobianDlaTrzeciegoPktCalk[1][0] + tab1[2][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		przedjakobianDlaTrzeciegoPktCalk[1][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaTrzeciegoPktCalk[1][1] = przedjakobianDlaTrzeciegoPktCalk[1][1] + tab2[2][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		//czwarty pkt calkowania
		przedjakobianDlaCzwartegoPktCalk[0][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaCzwartegoPktCalk[0][0] = przedjakobianDlaCzwartegoPktCalk[0][0] + tab1[3][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaCzwartegoPktCalk[0][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaCzwartegoPktCalk[0][1] = przedjakobianDlaCzwartegoPktCalk[0][1] + tab2[3][i] * grid.elements[0].pkt_calkowania[i].x;
		}
		przedjakobianDlaCzwartegoPktCalk[1][0] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaCzwartegoPktCalk[1][0] = przedjakobianDlaCzwartegoPktCalk[1][0] + tab1[3][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		przedjakobianDlaCzwartegoPktCalk[1][1] = 0;
		for (int i = 0; i < 4; i++)
		{
			przedjakobianDlaCzwartegoPktCalk[1][1] = przedjakobianDlaCzwartegoPktCalk[1][1] + tab2[3][i] * grid.elements[0].pkt_calkowania[i].y;
		}
		det1punkt = policzWyznacznik(przedjakobianDlaPierwszegoPktCalk);
		det2punkt = policzWyznacznik(przedjakobianDlaDrugiegoPktCalk);
		det3punkt = policzWyznacznik(przedjakobianDlaTrzeciegoPktCalk);
		det4punkt = policzWyznacznik(przedjakobianDlaTrzeciegoPktCalk);
		//cout << "det 1 :" << det1punkt << endl;
	//TU MAM TE WYZNACZNIKI KTORYCH UZYJE TEZ DO MACIERZY H ,WIEC TO JEST MOJ PUNKT WYJSCIA	 , funkcje ksztaltu tez mam wyliczone
			// no to chyba trzeba porozdzielac na 8 wektorow i te wektory stransponowac 
		for (int i = 0; i < 4; i++)
		{
			Ndla1punktudoC[i] = funkcjeKsztaltu[0][i];
			Ndla2punktudoC[i] = funkcjeKsztaltu[1][i];
			Ndla3punktudoC[i] = funkcjeKsztaltu[2][i];
			Ndla4punktudoC[i] = funkcjeKsztaltu[3][i];
		} //ok
		//N po transpon
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				NpoTranspondla1punktudoC[i][j] = transponujMnoz(Ndla1punktudoC)[i][j] * det1punkt * grid.elements[e].density * grid.elements[e].specific_heat;		//dN/dx * (dN/x) T
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				NpoTranspondla2punktudoC[i][j] = transponujMnoz(Ndla2punktudoC)[i][j] * det2punkt * grid.elements[e].density * grid.elements[e].specific_heat;		//dN/dx * (dN/x) T
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				NpoTranspondla3punktudoC[i][j] = transponujMnoz(Ndla3punktudoC)[i][j] * det3punkt * grid.elements[e].density * grid.elements[e].specific_heat;		//dN/dx * (dN/x) T
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				NpoTranspondla4punktudoC[i][j] = transponujMnoz(Ndla4punktudoC)[i][j] * det4punkt * grid.elements[e].density * grid.elements[e].specific_heat;		//dN/dx * (dN/x) T
			}
		}//tu wszystko ok, teraz tylko dodac te macierze i bedzie macierz C
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixC[i][j] = NpoTranspondla1punktudoC[i][j] + NpoTranspondla2punktudoC[i][j] + NpoTranspondla3punktudoC[i][j] + NpoTranspondla4punktudoC[i][j];
			}
		}//tu bylo wtracone liczenie macierzy C teraz wracam do macierzy H po dV
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				Jakobian1pkt[i][j] = (1 / det1punkt)* odwrocPrzedJakobian(przedjakobianDlaPierwszegoPktCalk)[i][j];
				Jakobian2pkt[i][j] = (1 / det2punkt)*odwrocPrzedJakobian(przedjakobianDlaDrugiegoPktCalk)[i][j];
				Jakobian3pkt[i][j] = (1 / det3punkt)*odwrocPrzedJakobian(przedjakobianDlaTrzeciegoPktCalk)[i][j];
				Jakobian4pkt[i][j] = (1 / det4punkt)*odwrocPrzedJakobian(przedjakobianDlaCzwartegoPktCalk)[i][j];
			}
		}
		//x 1pc
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDx[0][i] = Jakobian1pkt[0][0] * tab1[0][i] - Jakobian1pkt[0][1] * tab1[0][i];
		}
		//2pc
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDx[1][i] = Jakobian2pkt[0][0] * tab1[1][i] - Jakobian2pkt[0][1] * tab1[1][i];
		}
		//3pc
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDx[2][i] = Jakobian3pkt[0][0] * tab1[2][i] - Jakobian3pkt[0][1] * tab1[2][i];
		}
		//4pc
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDx[3][i] = Jakobian4pkt[0][0] * tab1[3][i] - Jakobian4pkt[0][1] * tab1[3][i];
		}
		//y
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDy[0][i] = (-1)*Jakobian1pkt[1][0] * tab2[0][i] + Jakobian1pkt[1][1] * tab2[0][i];
		}
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDy[1][i] = (-1)*Jakobian2pkt[1][0] * tab2[1][i] + Jakobian1pkt[1][1] * tab2[1][i];
		}
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDy[2][i] = (-1)*Jakobian3pkt[1][0] * tab2[2][i] + Jakobian1pkt[1][1] * tab2[2][i];
		}
		for (int i = 0; i < 4; i++)
		{
			TabDNPoDy[3][i] = (-1)*Jakobian4pkt[1][0] * tab2[3][i] + Jakobian1pkt[1][1] * tab2[3][i];
		}
		/*cout << "DN po DX" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << TabDNPoDx[i][j] << endl;
			}
		}
		cout << "Dn po DY" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << TabDNPoDy[i][j] << endl;
			}
		}*/
		// dn po dx i pody liczy swietnie 
		// na razie rozbijam moje dwie tablice dn po dx i dn po dy na 8 tablic 1d zeby przetransponowac i sprawzdzic wyniki
		for (int i = 0; i < 4; i++)
		{
			tab1dla1pktdnpodx[i] = TabDNPoDx[0][i];
		}
		for (int i = 0; i < 4; i++)
		{
			tab2dla2pktdnpodx[i] = TabDNPoDx[1][i];
		}
		for (int i = 0; i < 4; i++)
		{
			tab3dla3pktdnpodx[i] = TabDNPoDx[2][i];
		}
		for (int i = 0; i < 4; i++)
		{
			tab4dla4pktdnpodx[i] = TabDNPoDx[3][i];
		}

		for (int i = 0; i < 4; i++)
		{
			tab5dla1pktdnpody[i] = TabDNPoDy[0][i];
		}
		for (int i = 0; i < 4; i++)
		{
			tab6dla2pktdnpody[i] = TabDNPoDy[1][i];
		}
		for (int i = 0; i < 4; i++)
		{
			tab7dla3pktdnpody[i] = TabDNPoDy[2][i];
		}
		for (int i = 0; i < 4; i++)
		{
			tab8dla4pktdnpody[i] = TabDNPoDy[3][i];
		}
		
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji1pktpodx[i][j] = transponujMnoz(tab1dla1pktdnpodx)[i][j] * det1punkt;		//dN/dx * (dN/x) T  1pc
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji2pktpodx[i][j] = transponujMnoz(tab2dla2pktdnpodx)[i][j] * det2punkt;		//dN/dx * (dN/x) T  2pc
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji3pktpodx[i][j] = transponujMnoz(tab3dla3pktdnpodx)[i][j] * det3punkt;		//dN/dx * (dN/x) T  3pc
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji4pktpodx[i][j] = transponujMnoz(tab4dla4pktdnpodx)[i][j] * det4punkt;		//dN/dx * (dN/x) T   4pc
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji1pktpody[i][j] = transponujMnoz(tab5dla1pktdnpody)[i][j] * det1punkt;		//dN/dy * (dN/y) T
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji2pktpody[i][j] = transponujMnoz(tab6dla2pktdnpody)[i][j] * det2punkt;		//dN/dy * (dN/y) T
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji3pktpody[i][j] = transponujMnoz(tab7dla3pktdnpody)[i][j] * det3punkt;		//dN/dy * (dN/y) T
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tab2DpoTraspozycji4pktpody[i][j] = transponujMnoz(tab8dla4pktdnpody)[i][j] * det4punkt;		//dN/dy * (dN/y) T
			}
		}
		/*cout << "test dla 1 pkt calkowania po transpozycji " << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << tab2DpoTraspozycji4pktpody[i][j] << endl;
			}
			cout << endl;
		}  */   //czyli te moje tablice po transpozycji dzialaja  i pomnozeniu przez wyznacznik 
		 //teraz by sie przydalo dodac do siebie odpowiednio po dx i po dy dla poszczegolnych punktow calkowania i pomnozyc przez k
		
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tabPoDodaniuDyPlusDxrazkpkt1[i][j] = (tab2DpoTraspozycji1pktpodx[i][j] + tab2DpoTraspozycji1pktpody[i][j])*grid.elements[0].k;
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tabPoDodaniuDyPlusDxrazkpkt2[i][j] = (tab2DpoTraspozycji2pktpodx[i][j] + tab2DpoTraspozycji2pktpody[i][j])*grid.elements[0].k;
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tabPoDodaniuDyPlusDxrazkpkt3[i][j] = (tab2DpoTraspozycji3pktpodx[i][j] + tab2DpoTraspozycji3pktpody[i][j])*grid.elements[0].k;
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				tabPoDodaniuDyPlusDxrazkpkt4[i][j] = (tab2DpoTraspozycji4pktpodx[i][j] + tab2DpoTraspozycji4pktpody[i][j])*grid.elements[0].k;
			}
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixH[i][j] = tabPoDodaniuDyPlusDxrazkpkt1[i][j] + tabPoDodaniuDyPlusDxrazkpkt2[i][j] + tabPoDodaniuDyPlusDxrazkpkt3[i][j] + tabPoDodaniuDyPlusDxrazkpkt4[i][j];
			}
		}
		/*cout << "MACIERZ H DLA " <<e+1<< " elementu" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << grid.elements[e].MatrixH[i][j] << " \t";
			}
			cout << endl;
		}*/
	}
	/*for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHpow1[i][j] = 0;
				grid.elements[e].MatrixHpow2[i][j] = 0;
				grid.elements[e].MatrixHpow3[i][j] = 0;
				grid.elements[e].MatrixHpow4[i][j] = 0;
			}
		}
	}*/
	/*cout << "Macierz H bez ds " << endl;
	for (int e = 0; e < ne; e++)
	{
		cout << "e :" << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << grid.elements[e].MatrixH[i][j] << " \t";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}*/
	for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				 grid.elements[e].MatrixHDVplusDS[i][j] =0;
			}
		}
	}
	for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorP[i] = 0;
		}
	}
	//cout << "TERAZ MACIERZ H   i  wektor P pO DS DLA PIERWSZEGO SKRAJNEGO ELEMENTU " << endl;
	//pow 1
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[0].MatrixHpow1[i][j] = liczHpow1(funkcjeKsztaltu_dS, grid.elements[0].alfa, bokA)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[0].VectorPpow1[i] = liczPpow1(funkcjeKsztaltu_dS, grid.elements[0].ambient_temperature,grid.elements[0].alfa, bokA)[i];
	}
	//pow 4
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[0].MatrixHpow4[i][j] = liczHpow4(funkcjeKsztaltu_dS, grid.elements[0].alfa, bokB)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[0].VectorPpow4[i] = liczPpow4(funkcjeKsztaltu_dS, grid.elements[0].ambient_temperature,grid.elements[0].alfa, bokB)[i];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[0].MAtrixSumForFirstExtremeElement[i][j] = grid.elements[0].MatrixHpow1[i][j] + grid.elements[0].MatrixHpow4[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[0].VectorPSumForFirstExtremeElement[i] = grid.elements[0].VectorPpow1[i] + grid.elements[0].VectorPpow4[i];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[0].MatrixHDVplusDS[i][j] = grid.elements[0].MatrixH[i][j] + grid.elements[0].MAtrixSumForFirstExtremeElement[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[0].VectorP[i] += grid.elements[0].VectorPSumForFirstExtremeElement[i];
	}
	/*cout << "matrix  h z pierwszego skrajnego elementu" << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << grid.elements[0].MAtrixSumForFirstExtremeElement[i][j] << endl;        //Dziala !!
		}
		cout << endl;
	}*/
	/*cout << "wektor P dla pierwszego skrajnego elementu " << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << grid.elements[0].VectorPSumForFirstExtremeElement[i] << endl;
	}*/
	//TERAZ TO SAMO DLA 3 POZOSTALYCH SKRAJNYCH POWIERZCHNII
	//cout << "TERAZ MACIERZ H i wektor P pO DS DLA Drugiego SKRAJNEGO ELEMENTU " << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[nH-2].MatrixHpow4[i][j] = liczHpow4(funkcjeKsztaltu_dS, grid.elements[nH-2].alfa, bokB)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[nH-2].MatrixHpow3[i][j] = liczHpow3(funkcjeKsztaltu_dS, grid.elements[nH-2].alfa, bokA)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[nH-2].MAtrixSumForSecondExtremeElement[i][j] = grid.elements[nH-2].MatrixHpow3[i][j] + grid.elements[nH-2].MatrixHpow4[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[nH-2].VectorPpow3[i] = liczPpow3(funkcjeKsztaltu_dS, grid.elements[nH-2].ambient_temperature,grid.elements[nH-2].alfa, bokA)[i];
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[nH-2].VectorPpow4[i] = liczPpow4(funkcjeKsztaltu_dS, grid.elements[nH-2].ambient_temperature,grid.elements[nH-2].alfa, bokB)[i];
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[nH-2].VectorPSumForSecondExtremeElement[i] = grid.elements[nH-2].VectorPpow3[i] + grid.elements[nH-2].VectorPpow4[i];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[nH-2].MatrixHDVplusDS[i][j]= grid.elements[nH-2].MatrixH[i][j] + grid.elements[nH-2].MAtrixSumForSecondExtremeElement[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[nH-2].VectorP[i] += grid.elements[nH-2].VectorPSumForSecondExtremeElement[i];
	}
	//cout << "TERAZ MACIERZ H i wektor P pO DS DLA Trzciego SKRAJNEGO ELEMENTU " << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne-1].MatrixHpow2[i][j] = liczHpow2(funkcjeKsztaltu_dS, grid.elements[ne-1].alfa, bokB)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne-1].MatrixHpow3[i][j] = liczHpow3(funkcjeKsztaltu_dS, grid.elements[ne-1].alfa, bokA)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne-1].MAtrixSumForThirdExtremeElement[i][j] = grid.elements[ne-1].MatrixHpow2[i][j] + grid.elements[ne-1].MatrixHpow3[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne-1].VectorPpow3[i] = liczPpow3(funkcjeKsztaltu_dS, grid.elements[ne-1].ambient_temperature,grid.elements[ne-1].alfa, bokA)[i];
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne-1].VectorPpow2[i] = liczPpow2(funkcjeKsztaltu_dS, grid.elements[ne-1].ambient_temperature, grid.elements[ne - 1].alfa, bokB)[i];
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne-1].VectorPSumForThirdExtremeElement[i] = grid.elements[ne-1].VectorPpow3[i] + grid.elements[ne-1].VectorPpow2[i];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne -  1].MatrixHDVplusDS[i][j] = grid.elements[ne-1].MatrixH[i][j]+ grid.elements[ne-1].MAtrixSumForThirdExtremeElement[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne-1].VectorP[i] += grid.elements[ne-1].VectorPSumForThirdExtremeElement[i];
	}
	
	//cout << "TERAZ MACIERZ H i wektor P pO DS DLA czwaretgo SKRAJNEGO ELEMENTU " << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne - (nH-1)].MatrixHpow2[i][j] = liczHpow2(funkcjeKsztaltu_dS, grid.elements[ne- (nH-1)].alfa, bokB)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne - (nH-1)].MatrixHpow1[i][j] = liczHpow1(funkcjeKsztaltu_dS, grid.elements[ne - (nH-1)].alfa, bokA)[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne-(nH-1)].MAtrixSumForFourthExtremeElement[i][j] = grid.elements[ne - (nH-1)].MatrixHpow1[i][j] + grid.elements[ne- (nH-1)].MatrixHpow2[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne - (nH-1)].VectorPpow1[i] = liczPpow1(funkcjeKsztaltu_dS, grid.elements[ne -(nH- 1)].ambient_temperature,grid.elements[ne-(nH-1)].alfa, bokA)[i];
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne- (nH-1)].VectorPpow2[i] = liczPpow2(funkcjeKsztaltu_dS, grid.elements[ne - (nH-1)].ambient_temperature, grid.elements[ne - (nH - 1)].alfa, bokB)[i];
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne - (nH-1)].VectorPSumForFourthExtremeElement[i] = grid.elements[ne -(nH-1)].VectorPpow1[i] + grid.elements[ne -(nH-1)].VectorPpow2[i];
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			grid.elements[ne-(nH-1)].MatrixHDVplusDS[i][j] = grid.elements[ne-(nH-1)].MatrixH[i][j]+ grid.elements[ne-(nH-1)].MAtrixSumForFourthExtremeElement[i][j];
		}
	}
	for (int i = 0; i < 4; i++)
	{
		grid.elements[ne -(nH-1)].VectorP[i] += grid.elements[ne - (nH-1)].VectorPSumForFourthExtremeElement[i];
	}
	
	//teraz juz musze wziac te skrajne dla ktorych jest po 1 powierzchnii i od razu moglbym dodac to do macierzy H po dV
	//bok 1 (miedzy 1 i 2 skrajny)
	
	for (int e = 1; e < (nH - 2); e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHpow4[i][j] = liczHpow4(funkcjeKsztaltu_dS, grid.elements[e].alfa, bokB)[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorPpow4[i] = liczPpow4(funkcjeKsztaltu_dS, grid.elements[e].ambient_temperature,grid.elements[e].alfa, bokB)[i];
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHDVplusDS[i][j]= grid.elements[e].MatrixH[i][j]+ grid.elements[e].MatrixHpow4[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorP[i] += grid.elements[e].VectorPpow4[i];
		}
	} //ok
	//bok 2 (miedzy 2 i 3 skrajnym)
	for (int counter = 0, e = ((nH - 2) + (nH - 1)); counter < (nL - 3); e += (nH - 1) ,counter++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHpow3[i][j] = liczHpow3(funkcjeKsztaltu_dS, grid.elements[e].alfa, bokA)[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorPpow3[i] = liczPpow3(funkcjeKsztaltu_dS, grid.elements[e].ambient_temperature, grid.elements[e].alfa, bokA)[i];
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHDVplusDS[i][j] = grid.elements[e].MatrixH[i][j]+ grid.elements[e].MatrixHpow3[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorP[i] += grid.elements[e].VectorPpow3[i];
		}
	}  //git
	//bok 3 (miedzy 3 i 4 skrajnym)
	for (int e = ne - (nH-2); e < (ne - 1); e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHpow2[i][j] = liczHpow2(funkcjeKsztaltu_dS, grid.elements[e].alfa, bokB)[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorPpow2[i] = liczPpow2(funkcjeKsztaltu_dS, grid.elements[e].ambient_temperature, grid.elements[e].alfa, bokB)[i];
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHDVplusDS[i][j] = grid.elements[e].MatrixH[i][j]+ grid.elements[e].MatrixHpow2[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorP[i] += grid.elements[e].VectorPpow2[i];
		}
	} //git
	  //bok 4 (miedzy 4 i 1 skrajnym)
	for (int e = nH-1; e < (ne - (nH - 1)); e+=(nH-1))
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHpow1[i][j] = liczHpow1(funkcjeKsztaltu_dS, grid.elements[e].alfa, bokA)[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorPpow1[i] = liczPpow1(funkcjeKsztaltu_dS, grid.elements[e].ambient_temperature, grid.elements[e].alfa, bokA)[i];
		}
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				grid.elements[e].MatrixHDVplusDS[i][j] = grid.elements[e].MatrixH[i][j]+ grid.elements[e].MatrixHpow1[i][j];
			}
		}
		for (int i = 0; i < 4; i++)
		{
			grid.elements[e].VectorP[i] += grid.elements[e].VectorPpow1[i];
		}
	}
	/*cout << "MAcierz H " << endl;
	for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if (grid.elements[e].MatrixHDVplusDS[i][j] == 0)
					grid.elements[e].MatrixHDVplusDS[i][j] = grid.elements[e].MatrixH[i][j];
			}
			
		}
		/*cout << "e : " << e << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << grid.elements[e].MatrixHDVplusDS[i][j] << " \t";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}*/
	for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if (grid.elements[e].MatrixHDVplusDS[i][j] == 0)
					grid.elements[e].MatrixHDVplusDS[i][j] = grid.elements[e].MatrixH[i][j];
			}

		
		}
	}
	
	/*cout << "MACIERZe  C lokalne " << endl;
	for (int e = 0; e < ne; e++)
	{
		cout << "e : " << e << endl;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				cout << grid.elements[e].MatrixC[i][j] << " \t";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}*/

	/*cout << "Wektor P " << endl;
	for (int e = 0; e < ne; e++)
	{
		cout << "e : " << e << endl;
		for (int i = 0; i < 4; i++)
		{
			cout << grid.elements[e].VectorP[i] << endl;
		}
		cout << endl;
		cout << endl;
	}*/
	// I TERAZ MAM CALA MACIERZ [H]  i WEKTOR {P} i  MACIERZ [C]!!
	//robimy macierz  z daszkiem
	// proba Hl -> Hg
	for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				MatrixHGlobal[grid.elements[e].nodes_ids[i] -1][grid.elements[e].nodes_ids[j]-1 ] += grid.elements[e].MatrixHDVplusDS[i][j];
			}
		}
	}
	//proba CL->CG
	for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				MatrixCGlobal[grid.elements[e].nodes_ids[i] - 1][grid.elements[e].nodes_ids[j] - 1] += grid.elements[e].MatrixC[i][j];
			}
		}
	}
	//[H^]
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH*nL); j++)
		{
			MatrixCGlobalprzezdTau[i][j] =MatrixCGlobal[i][j] / deltaTau;
		}
	}
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH*nL); j++)
		{
			MatrixHzDaszkiem[i][j] = MatrixHGlobal[i][j] +MatrixCGlobalprzezdTau[i][j];
		}
	}
	
	// proba Pl -> Pg 
	for (int e = 0; e < ne; e++)
	{
		for (int i = 0; i < 4; i++)
		{
			VectorPGlobal[grid.elements[e].nodes_ids[i] - 1] += grid.elements[e].VectorP[i];
		}
	}
	//cout << "MAcierz H GLOBALNA :" << endl;

	/*for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH*nL); j++)
		{
			cout << MatrixHGlobal[i][j] << " ";
		}
		cout << endl;
	}*//**/  //Macierz H byla bez dS w wynikach , byc moze bedzie potrzebna
	/*cout << "MAcierz C GLOBALNA :" << endl;
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH*nL); j++)
		{
			cout << MatrixCGlobal[i][j] << " ";
		}
		cout << endl;
	} *///ok
	/*cout << "MAcierz H z daszkiem GLOBALNA :" << endl;
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH*nL); j++)
		{
			cout << MatrixHzDaszkiem[i][j] << " ";
		}
		cout << endl;
	} */ //// ok
	/*cout << "MAcierz C przez dTau GLOBALNA :" << endl;
	for (int i = 0; i < (nH*nL); i++)
	{
		for (int j = 0; j < (nH*nL); j++)
		{
			cout << MatrixCGlobalprzezdTau[i][j] << " ";
		}
		cout << endl;
	}*/
	/*cout << "Wektor P globalny " << endl;
	for (int i = 0; i < nH * nL; i++)
	{
		cout << VectorPGlobal[i] << endl;
	}*/
	//definuje wektor {t0} 
	double * wektorT0 = new double[nH*nL];
	for (int i = 0; i < nH* nL; i++)
	{
		wektorT0[i] = initialTemperature;
	}
	for (int i = 0; i < nH*nL; i++)
	{
		VectorPzDaszkiem[i] = mnozMacierziWektor(MatrixCGlobalprzezdTau, wektorT0, nH*nL)[i] - VectorPGlobal[i] ; //!!!!!!! mnozylem razy -1 zeby wektor P mial dodatnie wartosci
	}
	//sprawdzam czy wektor P z daszkiem dziala...
	//cout << "Wektor P z daszkiem " << endl;
	/*for (int i = 0; i < nH*nL; i++)
	{
		cout << VectorPzDaszkiem[i] << endl;    
		
	}*/
	//musze ten wektor zrobic dla kazdej iteracji czasu <- trza poprobowac i oczywiscie wtedy  WEKTORT0 BEDZIE BEZUTYCZNY BOZ TEGO WSZYSTKO BEDZIE BRANE	
	double** wektorT1 = new double*[number_iteration];
	for (int i = 0; i < number_iteration; i++)
	{
		wektorT1[i] = new double[nH*nL];
	}
	// deklaruje ze pierwsze wartosci dla 0 iteracji sa t0 = initial temperature
	for (int i = 0; i < nH*nL; i++)
	{
		wektorT1[0][i] = initialTemperature;
	}
	// TESTY KODU Z EDUINF TROCHE PRZEROBIONEGO
	//musze zrobic nowa macierz w ktorej obok macierzy H bedzie wektor wyrazow wolnych {P}
	double** MacierzHiPwjednym = new double*[nH * nL ];    //u nas 17x16
	for (int i = 0; i < nH * nL +1; i++)
	{
		MacierzHiPwjednym[i] = new double[nH*nL +1];
	}
	// juz wszystkie iteracje
	for (int ic = 0; ic < number_iteration; ic++)
	{
		for (int i = 0; i < nH*nL; i++)
		{
			VectorPzDaszkiem[i] = mnozMacierziWektor(MatrixCGlobalprzezdTau, wektorT1[ic], nH*nL)[i] - VectorPGlobal[i]; 
		}
		for (int i = 0; i < nH*nL; i++)
		{
			for (int j = 0; j < nH*nL + 1; j++)
			{
				MacierzHiPwjednym[i][j - 1] = MatrixHzDaszkiem[i][j - 1];

			}
			MacierzHiPwjednym[i][nH*nL] = VectorPzDaszkiem[i];
		}
		
		/*cout << "H  i P w jednym " << endl;
		for (int i = 0; i < nH*nL; i++)
		{
			for (int j = 0; j < nH*nL + 1; j++)
			{
				cout << MacierzHiPwjednym[i][j] << " ";   // dziala  , teraz  juz dac to do funkcji i powoli tego gaussa robic :) 

			}
			cout << endl;
		}*/
		//cout << "wektor temperatur  t1 dla "<<ic+1<<" iteracji" << endl;
		cout << "ic " << ic +1<< " :";
		if (gauss(nH*nL, MacierzHiPwjednym, wektorT1[ic]))
		{
			for (int i = 0; i< nH* nL; i++)
			{
				if (ic<number_iteration - 1)
					wektorT1[ic + 1][i] = wektorT1[ic][i];
				VectorPzDaszkiem[i] = 0;
			}
			/*for (int i = 0; i < nH * nL; i++)
			{
				cout<< wektorT1[ic][i] << endl;     // to teraz musze to zrobic dla kazdej iteracji
			}*/
			
		}

		
		
	}  cout << endl;
		double maxValue;
		double minValue;
		for (int i = 0; i < number_iteration; i++)
		{
		maxValue = liczMax(wektorT1, i, nH*nL);
		minValue = liczMin(wektorT1, i, nH*nL);
		cout << "t " << i + 1 << " min: " << minValue << endl;
		cout << "t " << i + 1 << " max: " << maxValue << endl;
		}
			   // trace dostep do poprzednich iteracji ale w wektorze P smialo moge tracic   
	//cout << "macierz H  dla 12 elementu " << grid.elements[11].MatrixH[0][0] << endl;  <- moge sie odniesc do macierzy w poszczgolnym elemencie 
	cout << "wspolrzedne skrajnego (lewego dolnego) wezla: " << endl;
	cout << "x: " << grid.nodes[0][0].x << endl;
	cout<< "y: " << grid.nodes[0][0].y << endl;
	cout << endl;
	cout << "wspolrzedne skrajnego (prawego gornego) wezla: " << endl;
	cout << "x: " << grid.nodes[nL-1][nH-1].x << endl;
	cout << "y: " << grid.nodes[nL-1][nH-1].y << endl;
	cout << "wezly dla ostatniego elementu: " << endl;
	plik.close();
	delete[] nodesIds;
	delete[] arrayNodesMain;
	delete[] element;
	delete[] MatrixCGlobal;
	delete[] MatrixHGlobal;
	delete[] MatrixHzDaszkiem;
	delete[] VectorPGlobal;
	delete[] MatrixCGlobalprzezdTau;
	delete[] MacierzHiPwjednym;
	delete[] VectorPzDaszkiem;
	delete[] wektorT0;
	delete[] wektorT1;
	system("pause");
	
}
