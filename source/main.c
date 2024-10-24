#include <3ds.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "font.h"

#define SCREENWIDTH 400
#define SCREENHEIGHT 240
#define ITERATIONSNEW3DS 256
#define ITERATIONSOLD3DS 128
#define INITIALTARGETX -0.5
#define INITIALTARGETY 0.0
#define INITIALSIZE 3.1
#define COLOURTABLESIZE 64
#define TIMELIMIT (1000000 * 88/100 / 60)
#define THREADSTACKSIZE (4 * 1024)

double targetX = INITIALTARGETX, targetY = INITIALTARGETY, size = INITIALSIZE;
volatile int iterations;

bool new3DS = true;
bool quit = false;
bool drawfull = true;
bool fixRowNext = true;
int colCycle = 0;
int rowCycle = 0;
u64 startTime;

#define MAXCORECOUNT 4
Thread threads[MAXCORECOUNT];
volatile int threadJob[MAXCORECOUNT];
u16* volatile threadBuffer;
volatile double threadPX;
volatile double threadPY;
int coreCount;
volatile bool runThreads = false;
volatile bool timeLeft = false;


u8 fontData[fontTilesLen];
u16 colourTable[COLOURTABLESIZE];
double colXValues[SCREENWIDTH];
double rowYValues[SCREENHEIGHT];
double oldColXValues[SCREENWIDTH];
double oldRowYValues[SCREENHEIGHT];
double colErrors[SCREENWIDTH];
double rowErrors[SCREENHEIGHT];
int colRand[SCREENWIDTH];
int rowRand[SCREENHEIGHT];

void InitRand(int* table, int count)
{
    for (int i = 0; i < count; ++i)
    {
        table[i] = i;
    }
    
    for (int i = count-1; i > 0; --i)
    {
        int j = rand() % (i + 1);
        
        int temp = table[i];
        table[i] = table[j];
        table[j] = temp;
    }
}

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

void Instructions(bool logo)
{
    iprintf("\x1b[1;1H");

    if (logo)
    {
        for (int i = 128; i < 224; ++i)
        {
            if (!(i%32))
            {
                iprintf("\n    ");
            }
            iprintf("%c",i);
        }
        ColourText(32, 8+3, 320-32, 32-3, 0, 0);
    }
    else
    {
        iprintf("\n\n\n");
    }

	iprintf(
        "\n\n\n\n"
        "      Circle-Pad : Move around\n"
        "             A/B : Zoom in/out\n"
        "             X/Y : Inc/dec iterations\n"
        "          Select : Reset view\n"
        "           Start : Exit\n"
        "\n\n"
        " Iteration count : \x1b[s                    \n"
        "        Real pos :                     \n"
        "        Imag pos :                     \n"
        "            Zoom :                     \n"
        "\n"
        "          Strips :                     \n"
        "            Time :                     \n"
        "\n\n\n\n"
        "            By Movie Vertigo\n"
        "        youtube.com/movievertigo\n"
        "        twitter.com/movievertigo"
    );
}

void Controls()
{
    circlePosition circlePos;
    scanKeys();

    int pressed = hidKeysDownRepeat();
    if (pressed & (KEY_X|KEY_Y))
    {
        if (pressed & KEY_X) ++iterations;
        if (pressed & KEY_Y) iterations = iterations > 1 ? iterations-1 : iterations;
        return;
    }

    pressed = hidKeysHeld();
    hidCircleRead(&circlePos);
    targetX += (circlePos.dx/16)*size/512;
    targetY += (circlePos.dy/16)*size/512;
    if (pressed & KEY_A) size *= 0.99;
    if (pressed & KEY_B) size /= 0.99;

    pressed = hidKeysDown();
    if (pressed & KEY_SELECT)
    {
        drawfull = true;
        size = INITIALSIZE; targetX = INITIALTARGETX; targetY = INITIALTARGETY;
        iterations = new3DS ? ITERATIONSNEW3DS : ITERATIONSOLD3DS;
    }
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
    double q = (x-0.25)*(x-0.25)+y2;
    if (q*(q+x-0.25) <= 0.25*y2 || (x+1.0)*(x+1.0)+y2 <= 0.0625) return iterations;

    int i = -1;
    while (++i < iterations && x2 + y2 < 4.0)
    {
        y = pY + 2 * x * y;
        x = pX + x2 - y2;
        x2 = x*x;
        y2 = y*y;
    }
    return i;
}

inline void threadcol(u16* ptr, double pX, int rowStart, int rowStep)
{
    for (int row = rowStart; timeLeft && row < SCREENHEIGHT; row += rowStep)
    {
        int i = mandelbrotpixel(pX, rowYValues[row]);
        *ptr = (i == iterations ? 0 : colourTable[i%COLOURTABLESIZE]);
        ptr += rowStep;
    }
}

inline void threadrow(u16* ptr, double pY, int colStart, int colStep)
{
    for (int col = colStart; timeLeft && col < SCREENWIDTH; col += colStep)
    {
        int i = mandelbrotpixel(colXValues[col], pY);
        *ptr = (i == iterations ? 0 : colourTable[i%COLOURTABLESIZE]);
        ptr += colStep * SCREENHEIGHT;
    }
}

void mandelbrotcol(u16* ptr, double pX)
{
    threadBuffer = ptr;
    threadPX = pX;
    __dsb();
    for (int core = 1; core < coreCount; ++core)
    {
        threadJob[core] = 1;
    }

    for (int row = 0; timeLeft && row < SCREENHEIGHT; row += coreCount)
    {
        int i = mandelbrotpixel(pX, rowYValues[row]);
        *ptr = (i == iterations ? 0 : colourTable[i%COLOURTABLESIZE]);
        ptr += coreCount;
    }

    for (int core = 1; core < coreCount; ++core)
    {
        while (threadJob[core]);
    }
}


void mandelbrotrow(u16* ptr, double pY)
{
    threadBuffer = ptr;
    threadPY = pY;
    __dsb();
    for (int core = 1; core < coreCount; ++core)
    {
        threadJob[core] = 2;
    }

    for (int col = 0; timeLeft && col < SCREENWIDTH; col += coreCount)
    {
        int i = mandelbrotpixel(colXValues[col], pY);
        *ptr = (i == iterations ? 0 : colourTable[i%COLOURTABLESIZE]);
        ptr += coreCount * SCREENHEIGHT;
    }

    for (int core = 1; core < coreCount; ++core)
    {
        while (threadJob[core]);
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

    if (worstError == 0.0)
    {
        worstCol = colRand[++colCycle % SCREENWIDTH];
    }

    double sizeX = (SCREENWIDTH * size) / SCREENHEIGHT;
    double minX = targetX - sizeX/2;
    double pX = minX + (worstCol * sizeX) / SCREENWIDTH;
    mandelbrotcol(buffer+worstCol*SCREENHEIGHT, pX);
    colXValues[worstCol] = pX;
    colErrors[worstCol] = 0.0;

    return worstError != 0.0;
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

    if (worstError == 0.0)
    {
        worstRow = rowRand[++rowCycle % SCREENHEIGHT];
    }

    double minY = targetY - size/2;
    double pY = minY + (worstRow * size) / SCREENHEIGHT;
    mandelbrotrow(buffer+worstRow, pY);
    rowYValues[worstRow] = pY;
    rowErrors[worstRow] = 0.0;

    return worstError != 0.0;
}

void threadFunc(void *arg)
{
	int core = (int)arg;

	while (runThreads)
	{
        switch (threadJob[core])
        {
            case 1:
            {
                threadcol(threadBuffer+core, threadPX, core, coreCount);
                threadJob[core] = 0;
                break;
            }
            case 2:
            {
                threadrow(threadBuffer+core*SCREENHEIGHT, threadPY, core, coreCount);
                threadJob[core] = 0;
                break;
            }
        }
	}
}

void InitThreads()
{
    for (int percent = 100; percent > 1 && APT_SetAppCpuTimeLimit(percent) != 0; --percent) {}
    coreCount = new3DS ? 4 : 2;

	runThreads = true;
    for (int core = 1; core < coreCount; ++core)
    {
        threadJob[core] = 0;
        __dsb();
        threads[core] = threadCreate(threadFunc, (void*)core, THREADSTACKSIZE, 0x3F, core, false);
    }
}

void StopThreads()
{
	runThreads = false;
    for (int core = 1; core < coreCount; ++core)
    {
		threadJoin(threads[core], U64_MAX);
		threadFree(threads[core]);
    }
}

void onVBlank()
{
	if (!drawfull)
    {
        gfxSwapBuffers();
    }
    timeLeft = drawfull;
}

int main(void)
{
    osSetSpeedupEnable(true);
    APT_CheckNew3DS(&new3DS);
    iterations = new3DS ? ITERATIONSNEW3DS : ITERATIONSOLD3DS;

    gfxInit(GSP_RGB565_OES,GSP_RGB565_OES,false);
    gfxSetDoubleBuffering(GFX_TOP, true);
    gfxSetDoubleBuffering(GFX_BOTTOM, false);
    gfxSet3D(false);
    hidSetRepeatParameters(16,1);

    InitConsole();
    InitRand(colRand, SCREENWIDTH);
    InitRand(rowRand, SCREENHEIGHT);
    InitColourTable();
    Instructions(true);
    InitThreads();

    u16* backBuffer;
    u16* frontBuffer = (u16*)gfxGetFramebuffer(GFX_TOP, GFX_LEFT, NULL, NULL);
    gfxSwapBuffers();

    gspSetEventCallback(GSPGPU_EVENT_VBlank0, onVBlank, NULL, false);

    int strips = 0;
    startTime = getusec();
	while (aptMainLoop() && !quit)
	{
        backBuffer = frontBuffer;
        frontBuffer = (u16*)gfxGetFramebuffer(GFX_TOP, GFX_LEFT, NULL, NULL);

        if (drawfull)
        {
            timeLeft = true;
            gfxSwapBuffers();
            fullmandelbrot(backBuffer);
            gspWaitForVBlank();
            gfxSwapBuffers();
            memcpy(frontBuffer, backBuffer, sizeof(u16) * SCREENWIDTH * SCREENHEIGHT);
            drawfull = false;
            while (timeLeft);
            Instructions(false);
        }
        else
        {
            scale(frontBuffer, backBuffer);
            timeLeft = true;
            strips = 0;
            bool anyFixesNeeded = false;
            while (timeLeft)
            {
                bool fixNeeded = fixRowNext ? fixworstrow(backBuffer) : fixworstcol(backBuffer);
                anyFixesNeeded |= fixNeeded;
                if (timeLeft)
                {
                    fixRowNext = !fixRowNext;
                    ++strips;
                }
            }
            if (!anyFixesNeeded)
            {
                // Allows home button and lid closing to work
                svcSleepThread(1);
            }
        }

        u64 usec = getusec()-startTime;
        startTime = getusec();
        printf("\x1b[u%d  \x1b[u\x1b[1B%.16f  \x1b[u\x1b[2B%.16f  \x1b[u\x1b[3B%lld  \x1b[u\x1b[5B%d  \x1b[u\x1b[6B%lldms %lldfps  ", iterations, targetX, targetY, (u64)(4/size), strips, usec/1000, (2000000/usec+1)/2);

        Controls();
	}

    StopThreads();
	gfxExit();
	return 0;
}
