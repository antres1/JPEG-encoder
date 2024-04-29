#ifndef NXNDCT_H_
#define NXNDCT_H_
#include <QDebug>
		

void GenerateDCTmatrix(double* DCTKernel, int order);

void DCT(const char input[], int16_t output[], int N, double* DCTKernel);

void IDCT(const short input[], char output[], int N, double* DCTKernel);
   
void extendBorders(char* input, int xSize, int ySize, int N, char** output, int* newXSize, int* newYSize);

void cropImage(char* input, int xSize, int ySize, char* output, int newXSize, int newYSize);

#endif // NXNDCT_H_
