#ifndef JPEG_H_
#define JPEG_H_

#include <QDebug>


void performJPEGEncoding(uchar Y_buff[], char U_buff[], char V_buff[], int xSize, int ySize, int compressionRate);


#endif // JPEG_H_

