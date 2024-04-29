#include "JPEG.h"
#include "NxNDCT.h"
#include <math.h>

#include "JPEGBitStreamWriter.h"


#define DEBUG(x) do{ qDebug() << #x << " = " << x;}while(0)



// quantization tables from JPEG Standard, Annex K
uint8_t QuantLuminance[8*8] =
    { 16, 11, 10, 16, 24, 40, 51, 61,
      12, 12, 14, 19, 26, 58, 60, 55,
      14, 13, 16, 24, 40, 57, 69, 56,
      14, 17, 22, 29, 51, 87, 80, 62,
      18, 22, 37, 56, 68,109,103, 77,
      24, 35, 55, 64, 81,104,113, 92,
      49, 64, 78, 87,103,121,120,101,
      72, 92, 95, 98,112,100,103, 99 };
uint8_t QuantChrominance[8*8] =
    { 17, 18, 24, 47, 99, 99, 99, 99,
      18, 21, 26, 66, 99, 99, 99, 99,
      24, 26, 56, 99, 99, 99, 99, 99,
      47, 66, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99,
      99, 99, 99, 99, 99, 99, 99, 99 };

void DCTUandV(const char input[], int16_t output[], int N, double* DCTKernel)
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
                sum = sum + DCTKernel[i*N+k] * (input[k*N+j]);
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

uint8_t quantQuality(uint8_t quant, uint8_t quality) {
    // Convert to an internal JPEG quality factor, formula taken from libjpeg
    int16_t q = quality <= 100 ? 5000 / quality : 200 - quality * 2;
    return clamp((quant * q + 50) / 100, 1, 255);
}

static void doZigZag1(int16_t block[], int N)
{
    /* TO DO */
    int currDiagonalWidth = 1;
    int i;
    int row;
    int col;
    int k = 0;

    int16_t* out = new int16_t[N*N];

    while (currDiagonalWidth<N)
    {
        for (i = 0; i<currDiagonalWidth; i++)
        {
            if (currDiagonalWidth % 2 == 1)
            {
                row = currDiagonalWidth - i - 1;
                col = i;
            }
            else
            {
                row = i;
                col = currDiagonalWidth - i - 1;
            }

            out[k] = block[row * N + col];
            k++;
        }
        currDiagonalWidth++;
    }
    while (currDiagonalWidth> 0)
    {
        for (i = currDiagonalWidth; i> 0; i--)
        {
            if (currDiagonalWidth % 2 == 1){
                row = N - currDiagonalWidth + i - 1;
                col = N - i;
            }
            else
            {
                row = N - i;
                col = N - currDiagonalWidth + i - 1;
            }
            out[k] = block[row * N + col];
            k++;
        }

        currDiagonalWidth--;
    }

    for(int i = 0; i < N*N; i++)
    {
        block[i] = out[i];
    }
    delete[] out;
}

static void doZigZag2(uint8_t block[],int N)
{
    int currDiagonalWidth = 1;
    int i;
    int row;
    int col;
    int k = 0;

    short* out = new short[N*N];

    while (currDiagonalWidth<N)
    {
        for (i = 0; i<currDiagonalWidth; i++)
        {
            if (currDiagonalWidth % 2 == 1)
            {
                row = currDiagonalWidth - i - 1;
                col = i;
            }
            else
            {
                row = i;
                col = currDiagonalWidth - i - 1;
            }

            out[k] = block[row * N + col];
            k++;
        }
        currDiagonalWidth++;
    }
    while (currDiagonalWidth> 0)
    {
        for (i = currDiagonalWidth; i> 0; i--)
        {
            if (currDiagonalWidth % 2 == 1){
                row = N - currDiagonalWidth + i - 1;
                col = N - i;
            }
            else
            {
                row = N - i;
                col = N - currDiagonalWidth + i - 1;
            }
            out[k] = block[row * N + col];
            k++;
        }

        currDiagonalWidth--;
    }

    for(int i = 0; i < N*N; i++)
    {
        block[i] = out[i];
    }
    delete[] out;
}

/* perform DCT */
void performDCT(char input[], int xSize, int ySize, int N, uint8_t* quantType)
{
	// TO DO
    // Create NxN buffer for one input block
    char* inBlock = new char[N*N];

    // Generate DCT kernel
    double* DCTKernel = new double[N*N];
    GenerateDCTmatrix(DCTKernel, N);

    // Create buffer for DCT coefficients
    int16_t* dctCoeffs = new int16_t[N*N];

    // Extend image borders
    int x_size2, y_size2;
    char* input2;
    extendBorders(input, xSize , ySize, N, &input2, &x_size2, &y_size2);

    for(int y = 0; y < y_size2/N; y++){
        for(int x = 0; x < x_size2/N; x++){

            // Fill input block buffer
            for(int j=0; j<N; j++){
                for(int i=0; i<N; i++){
                    inBlock[j*N+i] = input2[(N*y+j)*(x_size2)+N*x+i];
                }
            }

            // Invoke DCT
            DCTUandV(inBlock, dctCoeffs, N, DCTKernel);


            for(int j = 0; j < N; j++){
                for(int i = 0; i < N; i++){
                    dctCoeffs[j*N + i] = round(dctCoeffs[j*N+i]/quantType[j*N + i]);
                }
            }

            doZigZag1(dctCoeffs, 8);

            // Write output values to output image
            for(int j=0; j<N; j++){
                for(int i=0; i<N; i++){
                    input2[(N*y+j)*(x_size2)+N*x+i] = dctCoeffs[j*N+i];
                }
            }
        }
    }

    cropImage(input2, x_size2, y_size2, input, xSize, ySize);

    // Delete used memory buffers coefficients
    delete[] dctCoeffs;
    delete[] inBlock;
    delete[] DCTKernel;
    delete[] input2;
}

//JPEGBitStreamWriter streamer("example.jpg");
void performJPEGEncoding(uchar Y_buff[], char U_buff[], char V_buff[], int xSize, int ySize, int quality)
{
	DEBUG(quality);
	
	
    auto s = new JPEGBitStreamWriter("example.jpg");
	// TO DO
    int N = 8;
    uint8_t* luminance = new uint8_t[N*N];
    uint8_t* chrominance = new uint8_t[N*N];

    char Y_char[xSize*ySize];

    for(int x = 0; x < xSize; x++)
    {
        for(int y = 0; y < ySize; y++)
        {
            Y_char[y*xSize + x] = Y_buff[y*xSize + x] - 128;
        }
    }


    for(int i = 0;i < N; i++)
    {
        for(int k = 0; k < N; k++)
        {
            luminance[i*N + k] = quantQuality(QuantLuminance[i*N + k], quality);
            chrominance[i*N + k] = quantQuality(QuantChrominance[i*N + k], quality);
        }
    }


    performDCT(Y_char, xSize, ySize, 8, luminance);
    performDCT(U_buff, xSize/2, ySize/2, 8, chrominance);
    performDCT(V_buff, xSize/2, ySize/2, 8, chrominance);

    doZigZag2(luminance, N);
    doZigZag2(chrominance, N);



    s->writeHeader();
    s->writeQuantizationTables(luminance, chrominance);
    s->writeImageInfo(xSize, ySize);
    s->writeHuffmanTables();

    int16_t* inBlock1 = new int16_t[16*16];
    int16_t* inBlock2 = new int16_t[8*8];
    for(int y = 0; y < ySize/16; y++){
        for(int x = 0; x < xSize/16; x++){
            for (int j=0; j<16; j++)
            {
                for (int i=0; i<16; i++)
                {
                    inBlock1[j*16+i] = Y_char[(16*y+j)*xSize+16*x+i];
                }
            }
            for(int y2 = 0; y2 < 2; y2++){
                for(int x2 = 0; x2 < 2; x2++){
                    for(int j = 0; j < 8; j++){
                        for(int i = 0; i < 8; i++){
                            inBlock2[j*8 + i] = inBlock1[(8*y2+j)*16+8*x2+i];
                        }
                    }
                    s->writeBlockY(inBlock2);
                }
            }
            for(int j = 0; j < 8; j++){
                for(int i = 0; i < 8; i++){
                    inBlock2[j*8 + i] = U_buff[(8*y+j)*(xSize/2)+8*x+i];
                }
            }
            s->writeBlockU(inBlock2);
            for(int j = 0; j < 8; j++){
                for(int i = 0; i < 8; i++){
                    inBlock2[j*8 + i] = V_buff[(8*y+j)*(xSize/2)+8*x+i];
                }
            }
            s->writeBlockV(inBlock2);
        }
    }
    s->finishStream();

    delete[] inBlock1;
    delete[] inBlock2;
    delete s;
    delete[] luminance;
    delete[] chrominance;
}
