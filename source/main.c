#include <3ds.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "font.h"

#define SCREENWIDTH 400
#define SCREENHEIGHT 240
#define INITIALITERATIONS 256
#define INITIALTARGETX -0.5
#define INITIALTARGETY 0.0
#define INITIALSIZE 3.1
#define COLOURTABLESIZE 64
#define TIMELIMIT (1000000 * 88/100 / 60)
#define SPLITS 8

double targetX = INITIALTARGETX, targetY = INITIALTARGETY, size = INITIALSIZE;
int iterations = INITIALITERATIONS;

bool quit = false;
bool drawfull = true;
bool fixRowNext = true;
u64 startTime;

u8 fontData[fontTilesLen];
u16 colourTable[COLOURTABLESIZE];
double colXValues[SCREENWIDTH];
double rowYValues[SCREENHEIGHT];
double oldColXValues[SCREENWIDTH];
double oldRowYValues[SCREENHEIGHT];
double colErrors[SCREENWIDTH];
double rowErrors[SCREENHEIGHT];

void InitColourTable()
{
    for (int i = 0; i < COLOURTABLESIZE; ++i)
    {
        double u = (double)i / COLOURTABLESIZE;
        const u32 r = 255 * fmax(fmin(fabs(0.5-u)*6-1,1),0);
        const u32 g = 255 * fmax(fmin(fabs(0.5-fmod(u+2./3,1))*6-1,1),0);
        const u32 b = 255 * fmax(fmin(fabs(0.5-fmod(u+1./3,1))*6-1,1),0);
        colourTable[i] = ((r>>3)<<11) | ((g>>2)<<5) | ((b>>3)<<0);
    }
}

void InitConsole()
{
    consoleInit(GFX_BOTTOM, NULL);

    for (int i = 0; i < fontTilesLen; ++i)
    {
        u8 row = ((u8*)fontTiles)[i];
        fontData[i] = ((row &  1) << 7) | ((row &  2) << 5) | ((row &  4) << 3) | ((row &   8) << 1) | 
                      ((row & 16) >> 1) | ((row & 32) >> 3) | ((row & 64) >> 5) | ((row & 128) >> 7);
    }

    ConsoleFont font;
    font.gfx = fontData;
    font.asciiOffset = 32;
    font.numChars = 256-font.asciiOffset;
    consoleSetFont(NULL, &font);
}

void ColourText(int x1, int y1, int x2, int y2, int xo, int yo)
{
    u16* buffer = (u16*)gfxGetFramebuffer(GFX_BOTTOM, GFX_LEFT, NULL, NULL);
    for (int y = y1 + yo; y < y2; ++y)
    {
        for (int x = x1 + xo; x < x2; ++x)
        {
            double u = (double)(x-x1+y-y1) / (x2-x1);
            const u32 r = 255 * fmax(fmin(fabs(0.5-u)*6-1,1),0);
            const u32 g = 255 * fmax(fmin(fabs(0.5-fmod(u+2./3,1))*6-1,1),0);
            const u32 b = 255 * fmax(fmin(fabs(0.5-fmod(u+1./3,1))*6-1,1),0);
            u32 col = ((r>>3)<<11) | ((g>>2)<<5) | ((b>>3)<<0);

            if (buffer[x*240+239-y] > 0)
            {
                buffer[x*240+239-y] = col;
            }
        }
    }
}

void Instructions()
{
    consoleGetDefault()->cursorX = 0;
    consoleGetDefault()->cursorY = 0;

    for (int i = 128; i < 224; ++i)
    {
        if (!(i%32))
        {
            iprintf("\n    ");
        }
    	iprintf("%c",i);
    }
    ColourText(32, 8+3, 320-32, 32-3, 0, 0);
	printf(
        "\n\n"
        "      Circle-Pad : Move around\n"
        "             A/B : Zoom in/out\n"
        "             X/Y : Inc/dec iterations\n"
        "          Select : Reset view\n"
        "           Start : Exit\n"
        "\n"
        " Iteration count : \x1b[s                    \n"
        "        Real pos :                     \n"
        "        Imag pos :                     \n"
        "            Zoom :                     \n"
        "\n"
        "           Frame :                     \n"
        "           Vsync :                     \n"
        "\n"
        "            By Movie Vertigo\n"
        "        youtube.com/movievertigo\n"
        "        twitter.com/movievertigo"
    );
}

void Controls()
{
    circlePosition circlePos;
    scanKeys();

    int pressed = hidKeysHeld();

    hidCircleRead(&circlePos);
    targetX += (circlePos.dx/16)*size/512;
    targetY += (circlePos.dy/16)*size/512;
    if (pressed & KEY_A) size *= 0.99;
    if (pressed & KEY_B) size /= 0.99;

    pressed = hidKeysDown();
    if (pressed & KEY_SELECT) { size = INITIALSIZE; targetX = INITIALTARGETX; targetY = INITIALTARGETY; drawfull = true; }
    if (pressed & KEY_START) { quit = true; }
}

u64 getusec()
{
    return svcGetSystemTick() / CPU_TICKS_PER_USEC;
}

inline int mandelbrotpixel(double pX, double pY)
{
    double x = pX;
    double y = pY;
    double x2 = x*x;
    double y2 = y*y;

    int i = 0;
    while (++i < iterations && x2 + y2 < 4.0)
    {
        y = pY + 2 * x * y;
        x = pX + x2 - y2;
        x2 = x*x;
        y2 = y*y;
    }
    return i;
}

inline bool mandelbrotsplitcol(u16* ptr, double pX)
{
    for (int split = 0; split < SPLITS; ++split)
    {
        if (getusec()-startTime > TIMELIMIT)
        {
            return false;
        }

        int start = split * (SCREENHEIGHT/SPLITS);
        int end = (split+1) * (SCREENHEIGHT/SPLITS);

        for (int row = start; row < end; ++row)
        {
            int i = mandelbrotpixel(pX, rowYValues[row]);
            *ptr = i == iterations ? 0 : colourTable[i%COLOURTABLESIZE];
            ++ptr;
        }
    }

    return true;
}

inline bool mandelbrotsplitrow(u16* ptr, double pY)
{
    for (int split = 0; split < SPLITS; ++split)
    {
        if (getusec()-startTime > TIMELIMIT)
        {
            return false;
        }

        int start = split * (SCREENWIDTH/SPLITS);
        int end = (split+1) * (SCREENWIDTH/SPLITS);
    
        for (int col = start; col < end; ++col)
        {
            int i = mandelbrotpixel(colXValues[col], pY);
            *ptr = i == iterations ? 0 : colourTable[i%COLOURTABLESIZE];
            ptr += SCREENHEIGHT;
        }
    }

    return true;
}

inline void mandelbrotrow(u16* ptr, double pY)
{
    for (int col = 0; col < SCREENWIDTH; ++col)
    {
        int i = mandelbrotpixel(colXValues[col], pY);
        *ptr = i == iterations ? 0 : colourTable[i%COLOURTABLESIZE];
        ptr += SCREENHEIGHT;
    }
}

void fullmandelbrot(u16* buffer)
{
    double sizeX = (SCREENWIDTH * size) / SCREENHEIGHT;
    double minX = targetX - sizeX/2;
    double minY = targetY - size/2;
    u16* ptr = buffer;

    for (int col = 0; col < SCREENWIDTH; ++col)
    {
        colXValues[col] = minX + (col * sizeX) / SCREENWIDTH;
        colErrors[col] = 0.0;
    }

    for (int row = 0; row < SCREENHEIGHT; ++row)
    {
        double pY = minY + (row * size) / SCREENHEIGHT;
        rowYValues[row] = pY;
        rowErrors[row] = 0.0;
        mandelbrotrow(ptr, pY);
        ++ptr;
    }
}

void scale(u16* srcBuffer, u16* dstBuffer)
{
    double sizeX = (SCREENWIDTH * size) / SCREENHEIGHT;
    double minX = targetX - sizeX/2;
    double minY = targetY - size/2;

    memcpy(oldColXValues, colXValues, sizeof(double) * SCREENWIDTH);
    memcpy(oldRowYValues, rowYValues, sizeof(double) * SCREENHEIGHT);

    int srcCol = 0;
    int srcCols[SCREENWIDTH];
    for (int dstCol = 0; dstCol < SCREENWIDTH; ++dstCol)
    {
        double pX = minX + (dstCol * sizeX) / SCREENWIDTH;
        while (oldColXValues[srcCol] < pX && srcCol < SCREENWIDTH - 1)
        {
            ++srcCol;
        }
        if (srcCol > 0 && fabs(oldColXValues[srcCol-1]-pX) < fabs(oldColXValues[srcCol]-pX))
        {
            --srcCol;
        }
        srcCols[dstCol] = srcCol;
        colXValues[dstCol] = oldColXValues[srcCol];
        colErrors[dstCol] = fabs(colXValues[dstCol]-pX);
    }

    int srcRow = 0;
    int srcRows[SCREENHEIGHT];
    for (int dstRow = 0; dstRow < SCREENHEIGHT; ++dstRow)
    {
        double pY = minY + (dstRow * size) / SCREENHEIGHT;
        while (oldRowYValues[srcRow] < pY && srcRow < SCREENHEIGHT - 1)
        {
            ++srcRow;
        }
        if (srcRow > 0 && fabs(oldRowYValues[srcRow-1]-pY) < fabs(oldRowYValues[srcRow]-pY))
        {
            --srcRow;
        }
        srcRows[dstRow] = srcRow;
        rowYValues[dstRow] = oldRowYValues[srcRow];
        rowErrors[dstRow] = fabs(rowYValues[dstRow]-pY);
    }

    for (int dstCol = 0; dstCol < SCREENWIDTH; ++dstCol)
    {
        for (int dstRow = 0; dstRow < SCREENHEIGHT; ++dstRow)
        {
            dstBuffer[dstCol * SCREENHEIGHT + dstRow] = srcBuffer[srcCols[dstCol] * SCREENHEIGHT + srcRows[dstRow]];
        }
    }
}

bool fixworstcol(u16* buffer)
{
    double worstError = 0.0;
    int worstCol = 0;
    for (int col = 0; col < SCREENWIDTH; ++col)
    {
        if (colErrors[col] > worstError)
        {
            worstError = colErrors[col];
            worstCol = col;
        }
    }

    double sizeX = (SCREENWIDTH * size) / SCREENHEIGHT;
    double minX = targetX - sizeX/2;
    double pX = minX + (worstCol * sizeX) / SCREENWIDTH;
    if (mandelbrotsplitcol(buffer+worstCol*SCREENHEIGHT, pX))
    {
        colXValues[worstCol] = pX;
        colErrors[worstCol] = 0.0;
        return true;
    }
    else
    {
        return false;
    }
}

bool fixworstrow(u16* buffer)
{
    double worstError = 0.0;
    int worstRow = 0;
    for (int row = 0; row < SCREENHEIGHT; ++row)
    {
        if (rowErrors[row] > worstError)
        {
            worstError = rowErrors[row];
            worstRow = row;
        }
    }

    double minY = targetY - size/2;
    double pY = minY + (worstRow * size) / SCREENHEIGHT;
    if (mandelbrotsplitrow(buffer+worstRow, pY))
    {
        rowYValues[worstRow] = pY;
        rowErrors[worstRow] = 0.0;
        return true;
    }
    else
    {
        return false;
    }
}


int main(void)
{
    osSetSpeedupEnable(true);

    gfxInit(GSP_RGB565_OES,GSP_RGB565_OES,false);
    gfxSetDoubleBuffering(GFX_TOP, true);
    gfxSetDoubleBuffering(GFX_BOTTOM, false);
    gfxSet3D(false);

    InitConsole();
    InitColourTable();
    Instructions();

    u16* oldBuffer;
    u16* buffer = (u16*)gfxGetFramebuffer(GFX_TOP, GFX_LEFT, NULL, NULL);

    startTime = getusec();
	while (aptMainLoop() && !quit)
	{
        oldBuffer = buffer;
        buffer = (u16*)gfxGetFramebuffer(GFX_TOP, GFX_LEFT, NULL, NULL);

        if (drawfull)
        {
            fullmandelbrot(buffer);
            drawfull = false;
        }
        else
        {
            scale(oldBuffer, buffer);
            bool timeleft = true;
            while (timeleft)
            {
                if (fixRowNext)
                {
                    timeleft = fixworstrow(buffer);
                }
                else
                {
                    timeleft = fixworstcol(buffer);
                }
                if (timeleft)
                {
                    fixRowNext = !fixRowNext;
                }
            }
        }

        u64 usec = getusec()-startTime;
		gfxFlushBuffers();
		gfxSwapBuffers();
		gspWaitForVBlank();
        u64 usecvsync = getusec()-startTime;
        startTime = getusec();
        printf("\x1b[u%d  \x1b[u\x1b[1B%.16f  \x1b[u\x1b[2B%.16f  \x1b[u\x1b[3B%lld  \x1b[u\x1b[5B%lldms %lldfps  \x1b[u\x1b[6B%lldms %lldfps  ", iterations, targetX, targetY, (u64)(1/size), usec/1000, (2000000/usec+1)/2, usecvsync/1000, (2000000/usecvsync+1)/2);

        Controls();
	}

	gfxExit();
	return 0;
}
