
#include "ImageProcessing.h"
#include "ColorSpaces.h"
#include "JPEG.h"

#include <cmath>

#include <QDebug>
#include <QString>
#include <QImage>

void imageProcessingFun(const QString& progName, QImage& outImgs, const QImage& inImgs, const QVector<double>& params)
{
	// TO DO

    /* Create buffers for YUV image */
    uchar* Y_buff = new uchar[inImgs.width()*inImgs.height()];
    char* U_buff = new char[inImgs.width()*inImgs.height()/4];
    char* V_buff = new char[inImgs.width()*inImgs.height()/4];

	/* Create empty output image */
	outImgs = QImage(inImgs.width(), inImgs.height(), inImgs.format());

	/* Convert input image to YUV420 image */
    RGBtoYUV420(inImgs.bits(), inImgs.width(), inImgs.height(), Y_buff, U_buff, V_buff);

    if(progName == QString("JPEG Encoder"))
	{	
		/* Perform NxN DCT */
        performJPEGEncoding(Y_buff, U_buff, V_buff, inImgs.width(), inImgs.height(), params[0]);
	}

    /* Convert YUV image back to RGB */
    YUV420toRGB(Y_buff, U_buff, V_buff, inImgs.width(), inImgs.height(), outImgs.bits());

    outImgs = QImage("example.jpg");

	/* Delete used memory buffers */
    delete[] Y_buff;
    delete[] U_buff;
    delete[] V_buff;
}

