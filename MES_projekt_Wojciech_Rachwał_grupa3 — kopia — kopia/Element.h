#pragma once
#include "Node.h"

struct Element
{
	double k;
	double alfa;
	double ambient_temperature;
	double specific_heat;
	double density;
	int nodes_ids[4];
	double MatrixH[4][4];
	double MatrixHDVplusDS[4][4];
	double MatrixC[4][4];
	double MatrixHpow1[4][4];
	double MatrixHpow2[4][4];
	double MatrixHpow3[4][4];
	double MatrixHpow4[4][4];
	double MAtrixSumForFirstExtremeElement[4][4];
	double MAtrixSumForSecondExtremeElement[4][4];
	double MAtrixSumForThirdExtremeElement[4][4];
	double MAtrixSumForFourthExtremeElement[4][4];
	double VectorP[4];
	double VectorPpow1[4];
	double VectorPpow2[4];
	double VectorPpow3[4];
	double VectorPpow4[4];
	double VectorPSumForFirstExtremeElement[4];
	double VectorPSumForSecondExtremeElement[4];
	double VectorPSumForThirdExtremeElement[4];
	double VectorPSumForFourthExtremeElement[4];
	
	Node  pkt_calkowania[4];
	//Node  pkt_calkowania_dS[8];
	Element();
};

