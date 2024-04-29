#define _USE_MATH_DEFINES
#include <math.h>
#include "NxNDCT.h"

#define PI 3.14159265359

static char quantizationMatrix[64] =
{
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

void GenerateDCTmatrix(double* DCTKernel, int order)
{
    int i, j;
    int N;
    N = order;
    double alpha;
    double denominator;
    for (j = 0; j <= N - 1; j++)
    {
        DCTKernel[0, j] = sqrt(1 / (double)N);
    }
    alpha = sqrt(2 / (double)N);
    denominator = (double)2 * N;
    for (j = 0; j <= N - 1; j++)
        for (i = 1; i <= N - 1; i++)
        {
            DCTKernel[i*N + j] = alpha * cos(((2 * j + 1) * i * PI) / denominator);
        }
}

void DCT(const char input[], int16_t output[], int N, double* DCTKernel)
{
    double* temp = new double[N*N];
    double* DCTCoefficients = new double[N*N];

    double sum;
    for (int i = 0; i <= N - 1; i++)
    {
        for (int j = 0; j <= N - 1; j++)
        {
            sum = 0;
            for (int k = 0; k <= N - 1; k++)
            {
                sum = sum + DCTKernel[i*N+k] * (input[k*N+j]-128.0);
            }
            temp[i*N + j] = sum;
        }
    }

    for (int i = 0; i <= N - 1; i++)
    {
        for (int j = 0; j <= N - 1; j++)
        {
            sum = 0;
            for (int k = 0; k <= N - 1; k++)
            {
                sum = sum + temp[i*N+k] * DCTKernel[j*N+k];
            }
            DCTCoefficients[i*N+j] = sum;
        }
    }

    for(int i = 0; i < N*N; i++)
    {
        output[i] = floor(DCTCoefficients[i]+0.5);
    }

    delete[] temp;
    delete[] DCTCoefficients;

    return;
}



void IDCT(const short input[], char output[], int N, double* DCTKernel)
{
    /* TO DO */
}

void extendBorders(char* input, int xSize, int ySize, int N, char** p_output, int* newXSize, int* newYSize)
{
    int deltaX = N-xSize%N;
    int deltaY = N-ySize%N;

    *newXSize = xSize+deltaX;
    *newYSize = ySize+deltaY;

    char* output = new char[(xSize+deltaX)*(ySize+deltaY)];

    for (int i=0; i<ySize; i++)
    {
        memcpy(&output[i*(xSize+deltaX)], &input[i*(xSize)], xSize);
        if(deltaX != 0)
        {
            memset(&output[i*(xSize+deltaX)+xSize], output[i*(xSize+deltaX)+xSize - 1], deltaX);
        }
    }

    for (int i=0; i<deltaY; i++)
    {
        memcpy(&output[(i + ySize)*(xSize+deltaX)], &output[(ySize)*(xSize+deltaX)], xSize+deltaX);
    }

    *p_output = output;
}

void cropImage(char* input, int xSize, int ySize, char* output, int newXSize, int newYSize)
{
    for (int i=0; i<newYSize; i++)
    {
        memcpy(&output[i*(newXSize)], &input[i*(xSize)], newXSize);
    }
}

