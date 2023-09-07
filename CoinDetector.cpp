// CoinDetector.cpp : This file contains the 'main' function. Program execution begins and ends there.
// William Jiang
// 4/24/21

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <list>
#include <iterator>
#include <chrono>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;
using namespace std::chrono;

const int k = 3;
int ROWSIZE = 600;
int COLSIZE = 600;
int** hysteresis;

class Point //a class that stores the x and y coordinates of a 2-D point in the Cartesian plane
{
private:
    int x;
    int y;

public:
    Point(int x, int y)
    {
        this->x = x;
        this->y = y;
    }
    ~Point()
    {

    }
    int getX()
    {
        return x;
    }
    int getY()
    {
        return y;
    }
};

class Pixel
{
private:
    int point[k]; // To store k dimensional point 

public:
    Pixel(int point[k])
    {
        for (int i = 0; i < k; ++i)
            this->point[i] = point[i];
    }
    Pixel()
    {
        for (int i = 0; i < k; ++i)
            this->point[i] = 0;
    }
    int* getCoordinates()
    {
        return point;
    }
    int getR()
    {
        return point[0];
    }
    int getG()
    {
        return point[1];
    }
    int getB()
    {
        return point[2];
    }
    void setColors(int r, int g, int b)
    {
        point[0] = r;
        point[1] = g;
        point[2] = b;
    }
    bool equal_to(Pixel p)
    {
        for (int i = 0; i < k; ++i)
            if (this->point[i] != p.getCoordinates()[i])
                return false;
        return true;
    }
    int compare(Pixel p, int dim)
    {
        if (this->point[dim] < p.getCoordinates()[dim])
            return -1;
        if (this->point[dim] == p.getCoordinates()[dim])
            return 0;
        return 1;
    }
    bool operator==(const Pixel& p) const
    {
        return point[0] == p.point[0] && point[1] == p.point[1] && point[2] == p.point[2];
    }
};

class Drawing //a class that can draw circles and lines on a 2d array
{
public:
    Drawing();
    ~Drawing();
    //draws a line between a(x1,y1) and b(x2,y2) satisfying x2 > x1,y2 > y1,|x2 - x1| > |y2 - y1|
    void linedownrightx(Point a, Point b);
    //draws a line between a(x1,y1) and b(x2,y2) satisfying x2 > x1 && y1 > y2, | x2 - x1 | > | y2 - y1 |
    void linedownleftx(Point a, Point b);
    //draws a line between a(x1,y1) and b(x2,y2) satisfying x2 > x1 && y2 > y1, |x2 - x1| < |y2 - y1|
    void linedownrighty(Point a, Point b);
    //draws a line between a(x1,y1) and b(x2,y2) satisfying x2 > x1 && y1 > y2, |x2 - x1| < |y2 - y1|
    void linedownlefty(Point a, Point b);
    //draws line between a1 and b1 using Bresenham's Algorithm
    void drawline(Point a1, Point b1);
    //draws circle with center c and radius r
    void drawcircle(Point c, int r);
    //draws the extension of the line between a and b that covers the whole screen 
    void extendedLine(Point a, Point b);
    //draws the points a,b,c, and d as circles with radius 2 and squares with vertices stored in arr1, then writes arr to ppm file
    void run(double** direction, string** r1, Pixel** pixel);
    int edgecount(Point c, int r, string** arr1);
private:
    int** arr;
    string** ppm;
};

Drawing::Drawing()
{

    arr = new int* [ROWSIZE];            //allocate ROWSIZE pointers/rows
    for (int i = 0; i < ROWSIZE; ++i)   //allocate COLSIZE columns
        arr[i] = new int[COLSIZE];

    for (int i = 0; i < ROWSIZE; i++)
        for (int j = 0; j < COLSIZE; j++)
            arr[i][j] = 0;

    ppm = new string * [ROWSIZE];
    for (int i = 0; i < ROWSIZE; ++i)
        ppm[i] = new string[COLSIZE];

    for (int i = 0; i < ROWSIZE; i++)
        for (int j = 0; j < COLSIZE; j++)
            ppm[i][j] = "0 0 0 ";

}
Drawing::~Drawing()
{
    for (int i = 0; i < ROWSIZE; i++)
        delete[] arr[i];
    delete[] arr;
}

void Drawing::linedownrightx(Point a, Point b)
{
    int x1 = a.getX();
    int x2 = b.getX();
    int y1 = a.getY();
    int y2 = b.getY();
    int x = x2 - x1;
    int y = y2 - y1;
    int j = y1;
    int e = y - x;
    for (int i = x1; i <= x2; i++)
    {
        arr[j][i] = arr[j][i] + 1;
        if (e >= 0)
        {
            j += 1;
            e -= x;
        }
        e += y;
    }
}
void Drawing::linedownleftx(Point a, Point b)
{
    int x1 = a.getX();
    int x2 = b.getX();
    int y1 = a.getY();
    int y2 = b.getY();
    int x = x2 - x1;
    int y = y1 - y2;
    int j = y1;
    int e = y - x;
    for (int i = x1; i <= x2; i++)
    {
        arr[j][i] = arr[j][i] + 1;
        if (e >= 0)
        {
            j -= 1;
            e -= x;
        }
        e += y;
    }
}
void Drawing::linedownrighty(Point a, Point b)
{
    int x1 = a.getX();
    int x2 = b.getX();
    int y1 = a.getY();
    int y2 = b.getY();
    int x = x2 - x1;
    int y = y2 - y1;
    int j = x1;
    int e = x - y;
    for (int i = y1; i <= y2; i++)
    {
        arr[i][j] = arr[i][j] + 1;
        if (e >= 0)
        {
            j += 1;
            e -= y;
        }
        e += x;
    }
}
void Drawing::linedownlefty(Point a, Point b)
{
    int x1 = a.getX();
    int x2 = b.getX();
    int y1 = a.getY();
    int y2 = b.getY();
    int x = x2 - x1;
    int y = y1 - y2;
    int j = x2;
    int e = x - y;
    for (int i = y2; i <= y1; i++)
    {
        arr[i][j] = arr[i][j] + 1;
        if (e >= 0)
        {
            j -= 1;
            e -= y;
        }
        e += x;
    }
}

void Drawing::drawline(Point a, Point b) //draws line between a1 and b1 using Bresenham's Algorithm
{
    int ax = a.getX();
    int bx = b.getX();
    int ay = a.getY();
    int by = b.getY();
    if (abs(bx - ax) >= abs(by - ay))
    {
        if (bx > ax && by > ay)
            linedownrightx(a, b);
        else if (bx < ax && by < ay)
            linedownrightx(b, a);
        else if (bx > ax && by < ay)
            linedownleftx(a, b);
        else if (bx < ax && by > ay)
            linedownleftx(b, a);
        else
        {
            for (int i = std::min(ax, bx); i <= std::max(ax, bx); i++)
                arr[ay][i] = arr[ay][i] + 1;
        }
    }
    else
    {
        if (bx > ax && by > ay)
            linedownrighty(a, b);
        else if (bx < ax && by < ay)
            linedownrighty(b, a);
        else if (bx > ax && by < ay)
            linedownlefty(a, b);
        else if (bx < ax && by > ay)
            linedownlefty(b, a);
        else
            for (int i = std::min(ay, by); i <= std::max(ay, by); i++)
                arr[i][ax] = arr[i][ax] + 1;
    }
}

void Drawing::drawcircle(Point c, int r) //draws red circle with center c and radius r by making each cell in the array "0 0 0 "
{
    int cx = c.getX();
    int cy = c.getY();
    int x, y, xmax, y2, y2_new, ty;
    xmax = (int)(r * 0.70710678); // maximum x at radius/sqrt(2)
    y = r;
    y2 = y * y;
    ty = (2 * y) - 1; y2_new = y2;
    for (x = 0; x <= xmax; x++)
    {
        if ((y2 - y2_new) >= ty)
        {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        if ((cy + x < ROWSIZE) && (cy + x >= 0))
        {
            if ((cx + y < COLSIZE) && (cx + y >= 0))
                ppm[cy + x][cx + y] = "255 0 0 ";
            if ((cx - y >= 0) && (cx - y < COLSIZE))
                ppm[cy + x][cx - y] = "255 0 0 ";
        }
        if ((cy - x >= 0) && (cy - x < ROWSIZE))
        {
            if ((cx + y < COLSIZE) && (cx + y >= 0))
                ppm[cy - x][cx + y] = "255 0 0 ";
            if ((cx - y >= 0) && (cx - y < COLSIZE))
                ppm[cy - x][cx - y] = "255 0 0 ";
        }
        if ((cy + y < ROWSIZE) && (cy + y >= 0))
        {
            if ((cx + x < COLSIZE) && (cx + x >= 0))
                ppm[cy + y][cx + x] = "255 0 0 ";
            if ((cx - x >= 0) && (cx - x < COLSIZE))
                ppm[cy + y][cx - x] = "255 0 0 ";
        }
        if ((cy - y >= 0) && (cy - y < ROWSIZE))
        {
            if ((cx + x < COLSIZE) && (cx + x >= 0))
                ppm[cy - y][cx + x] = "255 0 0 ";
            if ((cx - x >= 0) && (cx - x < COLSIZE))
                ppm[cy - y][cx - x] = "255 0 0 ";
        }
        y2_new -= (2 * x) - 3;
    }
}

int Drawing::edgecount(Point c, int r,string** arr1)
{
    int count = 0;
    int cx = c.getX();
    int cy = c.getY();
    int x, y, xmax, y2, y2_new, ty;
    xmax = (int)(r * 0.70710678); // maximum x at radius/sqrt(2)
    y = r;
    y2 = y * y;
    ty = (2 * y) - 1; y2_new = y2;
    for (x = 0; x <= xmax; x++)
    {
        if ((y2 - y2_new) >= ty)
        {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }
        if ((cy + x < ROWSIZE) && (cy + x >= 0))
        {
            if ((cx + y < COLSIZE) && (cx + y >= 0))
                if (arr1[cy + x][cx + y] == "255 255 255")
                    count++;
            if ((cx - y >= 0) && (cx - y < COLSIZE))
                if (arr1[cy + x][cx - y] == "255 255 255 ")
                    count++;
        }
        if ((cy - x >= 0) && (cy - x < ROWSIZE))
        {
            if ((cx + y < COLSIZE) && (cx + y >= 0))
                if (arr1[cy - x][cx + y] == "255 255 255 ")
                    count++;
            if ((cx - y >= 0) && (cx - y < COLSIZE))
                if (arr1[cy - x][cx - y] == "255 255 255 ")
                    count++;
        }
        if ((cy + y < ROWSIZE) && (cy + y >= 0))
        {
            if ((cx + x < COLSIZE) && (cx + x >= 0))
                if (arr1[cy + y][cx + x] == "255 255 255 ")
                    count++;
            if ((cx - x >= 0) && (cx - x < COLSIZE))
                if (arr1[cy + y][cx - x] == "255 255 255 ")
                    count++;
        }
        if ((cy - y >= 0) && (cy - y < ROWSIZE))
        {
            if ((cx + x < COLSIZE) && (cx + x >= 0))
                if (arr1[cy - y][cx + x] == "255 255 255 ")
                    count++;
            if ((cx - x >= 0) && (cx - x < COLSIZE))
                if (arr1[cy - y][cx - x] == "255 255 255 ")
                    count++;
        }
        y2_new -= (2 * x) - 3;
    }
    return count;
}

void Drawing::extendedLine(Point a, Point b) //draws the extension of the line between a,b that covers the whole screen
{
    double x1 = a.getX();
    double x2 = b.getX();
    double y1 = a.getY();
    double y2 = b.getY();

    int leftx = 0;
    int lefty = (-x1 * (y1 - y2) / (x1 - x2) + y1);
    int topx = (-y1 * (x1 - x2) / (y1 - y2) + x1);
    int topy = 0;
    int rightx = COLSIZE - 1;
    int righty = (ROWSIZE - x1) * (y1 - y2) / (x1 - x2) + y1;
    int bottomx = (COLSIZE - y1) * (x1 - x2) / (y1 - y2) + x1;
    int bottomy = ROWSIZE - 1;

    Point l(leftx, lefty);
    Point t(topx, topy);
    Point r(rightx, righty);
    Point bot(bottomx, bottomy);

    if (0 <= lefty && lefty < ROWSIZE)
    {
        if (0 <= topx && topx < COLSIZE)
        {
            drawline(l, t);
        }
        else if (0 <= righty && righty < ROWSIZE)
        {
            drawline(l, r);
        }
        else if (0 <= bottomx && bottomx < COLSIZE)
        {
            drawline(l, bot);
        }
    }
    else if (0 <= topx && topx < COLSIZE)
    {
        if (0 <= righty && righty < ROWSIZE)
        {
            drawline(t, r);
        }
        else if (0 <= bottomx && bottomx < COLSIZE)
        {
            drawline(t, bot);
        }
    }
    else if (0 <= bottomx && bottomx < COLSIZE && 0 <= righty && righty < ROWSIZE)
    {
        drawline(r, bot);
    }
}

void Drawing::run(double** direction, string** arr1, Pixel** pixel)
{
    Point a(0, 0);
    Point b(0, 0);
    int y;
    string temp;
    for (int i = 0; i < ROWSIZE; i++)
        for (int j = 0; j < COLSIZE; j++)
        {
            temp = to_string(pixel[i][j].getR()) + " ";
            temp += to_string(pixel[i][j].getG());
            temp += " ";
            temp += to_string(pixel[i][j].getB());
            temp += " ";
            ppm[i][j] = temp;
        }
    for (int i = 0; i < ROWSIZE; i++)
    {
        for (int j = 0; j < COLSIZE; j++)
        {
            if (arr1[i][j] == "1 1 1 ")
            {
                a = Point(i, j);
                if (i > 200)
                {
                    y = (int)(j - 200 * direction[i][j]);
                    if (y > 0 && y < COLSIZE)
                    {
                        b = Point(i - 200, y);
                        extendedLine(a, b);
                    }
                }
                
                else
                {
                    y = (int)(j + 200 * direction[i][j]);
                    if (y > 0 && y < COLSIZE)
                    {
                        b = Point(i + 200, y);
                        extendedLine(a, b);
                    }
                }
            }
        }
    }


    for (int i = 0; i < ROWSIZE; i++)
        for (int j = 0; j < COLSIZE; j++)
            if (arr[i][j] > 120)
            {
                a = Point(j, i);
                for (int r = 1; r <= 50; r++)
                    drawcircle(a, r);
            }
    int max = 0;
    for (int i = 0; i < ROWSIZE; i++)
        for (int j = 0; j < COLSIZE; j++)
            if (arr[i][j] > max)
                max = arr[i][j];
    std::ofstream myfile;
    myfile.open("imagev.ppm");
    myfile << "P3" << " " << ROWSIZE << " " << COLSIZE << " " << max << "\n";
    //writes array to ppm file
    for (int i = 0; i < COLSIZE; i++)
    {
        for (int j = 0; j < ROWSIZE; j++)
        {
            myfile << arr[j][i] << " " << arr[j][i] << " " << arr[j][i] << " ";
        }
        myfile << "\n";
    }
    myfile.close();
    myfile.open("imageCC.ppm");
    myfile << "P3" << " " << ROWSIZE << " " << COLSIZE << " " << 255 << "\n";
    for (int i = 0; i < COLSIZE; i++)
    {
        for (int j = 0; j < ROWSIZE; j++)
        {
            myfile << ppm[j][i];
        }
        myfile << "\n";
    }
    myfile.close();

    for (int i = 0; i < COLSIZE; i++)
    {
        for (int j = 0; j < ROWSIZE; j++)
        {
            ppm[i][j] = arr1[i][j];
        }
    }
    myfile.close();
    int min = 100;
    max = 350;
    int penny = 0;
    int nickel = 0;
    int dime = 0;
    int quarter = 0;
    int dollar = 0;
    for(int i = 0; i < ROWSIZE; ++i)
        for (int j = 0; j < COLSIZE; ++j)
        {
            if (arr[i][j] > 120)
                for (int k = min; k < max; k++)
                {
                    if (edgecount(Point(i, j), k, arr1) > 200)
                    {
                        if (k < 120)
                            dime++;
                        else if (k < 140)
                            penny++;
                        else if (k < 170)
                            nickel++;
                        else if (k < 240)
                            quarter++;
                        else
                            dollar++;
                        drawcircle(Point(i, j), k);
                    }
                }
        }
    myfile.open("coins.ppm");
    myfile << "P3" << " " << ROWSIZE << " " << COLSIZE << " " << 255 << "\n";
    for (int i = 0; i < COLSIZE; i++)
    {
        for (int j = 0; j < ROWSIZE; j++)
        {
            myfile << ppm[j][i];
        }
        myfile << "\n";
    }
    myfile.close();
    cout << dollar << " silver dollars, " << quarter << " quarters, " << dime << " dimes, " << nickel << " nickels, " << penny << " pennies" << endl;
    double total = dollar + 0.25 * quarter + 0.1 * dime + 0.05 * nickel + 0.01 * penny;
    cout << "Total sum: $" << total;

    myfile.open("results.txt");
    myfile << dollar << " silver dollars, " << quarter << " quarters, " << dime << " dimes, " << nickel << " nickels, " << penny << " pennies" << endl;
    myfile << "Total sum: $" << total;
}



//recursive
void edgeFill(int i, int j) //replaces any 1's and 2's adjacent to a given position (i,j) with a temporary character (3)
{
    if (i < 0 || i >= ROWSIZE || j < 0 || j >= COLSIZE)
        return;
    if (hysteresis[i][j] == 3 || hysteresis[i][j] == 0)
        return;
    else if (hysteresis[i][j] == 1 || hysteresis[i][j] == 2)
        hysteresis[i][j] = 3;
    else
        cout << hysteresis[i][j];
    edgeFill(i + 1, j + 1);
    edgeFill(i + 1, j);
    edgeFill(i + 1, j - 1);
    edgeFill(i - 1, j + 1);
    edgeFill(i - 1, j);
    edgeFill(i - 1, j - 1);
    edgeFill(i, j + 1);
    edgeFill(i, j - 1);
}

void part2()
{
    ifstream myfile;
    myfile.open("image.ppm");
    string s;
    getline(myfile, s);
    getline(myfile, s);
    string slength = s.substr(0, s.find(" "));
    string swidth = s.substr(s.find(" ") + 1);
    int length = stoi(slength);
    int width = stoi(swidth);
    ROWSIZE = length;
    COLSIZE = width;
    Pixel** points = new Pixel * [ROWSIZE];
    for (int i = 0; i < ROWSIZE; ++i)
        points[i] = new Pixel[COLSIZE];
    getline(myfile, s);
    int i, j, k;
    int colcount = -1;
    int rowcount = 0;
    while (myfile >> i >> j >> k)
    {
        colcount++;
        if (colcount % length == 0 && colcount != 0)
        {
            rowcount++;
            colcount = colcount % length;
        }

        points[colcount][rowcount].setColors(i, j, k);

    }
    myfile.close();

    int** angle = new int* [ROWSIZE]; //array storing arctans of gradients
    for (int i = 0; i < ROWSIZE; ++i)
        angle[i] = new int[COLSIZE];

    double** direction = new double* [ROWSIZE];
    for (int i = 0; i < ROWSIZE; ++i)
        direction[i] = new double[COLSIZE];

    int** sobel = new int* [ROWSIZE]; //array with gradients
    for (int i = 0; i < ROWSIZE; ++i)
        sobel[i] = new int[COLSIZE];

    int** nms = new int* [ROWSIZE]; //array storing values after non maximum suppression
    for (int i = 0; i < ROWSIZE; ++i)
        nms[i] = new int[COLSIZE];

    hysteresis = new int* [ROWSIZE]; //array storing values after hysteresis
    for (int i = 0; i < ROWSIZE; ++i)
        hysteresis[i] = new int[COLSIZE];

    for (int i = 1; i < ROWSIZE - 1; ++i)
    {
        for (int j = 1; j < COLSIZE - 1; j++)
        {
            //calculates sobel operator for x and y
            int x = -points[i - 1][j - 1].getR() + points[i + 1][j - 1].getR()
                - 2 * points[i - 1][j].getR() + 2 * points[i + 1][j].getR() - points[i - 1][j + 1].getR() + points[i + 1][j + 1].getR();
            int y = -points[i - 1][j - 1].getR() + points[i - 1][j + 1].getR()
                - 2 * points[i][j - 1].getR() + 2 * points[i][j + 1].getR() - points[i + 1][j - 1].getR() + points[i + 1][j + 1].getR();
            sobel[i][j] = (int)sqrt(x * x + y * y);

            direction[i][j] = (double)y / x;

            //calculates arctan for non maximum suppression
            double thisAngle = atan2(y, x) / 3.14159 * 180;

            //rounds angle to multiple of 45 degrees for non maximum suppression
            if (((thisAngle <= 22.5) && (thisAngle > -22.5)) || (thisAngle > 157.5) || (thisAngle <= -157.5))
                angle[i][j] = 0;
            if (((thisAngle > 22.5) && (thisAngle <= 67.5)) || ((thisAngle <= -112.5) && (thisAngle > -157.5)))
                angle[i][j] = 45;
            if (((thisAngle > 67.5) && (thisAngle <= 112.5)) || ((thisAngle <= -67.5) && (thisAngle > -112.5)))
                angle[i][j] = 90;
            if (((thisAngle > 112.5) && (thisAngle <= 157.5)) || ((thisAngle <= -22.5) && (thisAngle > -67.5)))
                angle[i][j] = 135;
        }
    }

    //implements non maximum suppression
    for (int i = 1; i < ROWSIZE - 1; ++i) {
        for (int j = 1; j < COLSIZE - 1; ++j) {
            switch (angle[i][j]) {
            case 0:
                if (sobel[i][j] >= sobel[i - 1][j] && sobel[i][j] >= sobel[i + 1][j])
                    nms[i][j] = 1;
                else
                    nms[i][j] = 0;
                break;
            case 45:
                if (sobel[i][j] >= sobel[i - 1][j - 1] && sobel[i][j] >= sobel[i + 1][j + 1])
                    nms[i][j] = 1;
                else
                    nms[i][j] = 0;
                break;
            case 90:
                if (sobel[i][j] >= sobel[i][j - 1] && sobel[i][j] >= sobel[i][j + 1])
                    nms[i][j] = 1;
                else
                    nms[i][j] = 0;
                break;
            case 135:
                if (sobel[i][j] >= sobel[i - 1][j + 1] && sobel[i][j] >= sobel[i + 1][j - 1])
                    nms[i][j] = 1;
                else
                    nms[i][j] = 0;
                break;
            default:
                cout << "something wrong";
            }
        }
    }

    for (int i = 0; i < ROWSIZE; ++i)
    {
        nms[i][0] = 0;
        nms[i][COLSIZE - 1] = 0;
    }
    for (int j = 0; j < COLSIZE; ++j)
    {
        nms[0][j] = 0;
        nms[ROWSIZE - 1][j] = 0;
    }

    for (int i = 0; i < ROWSIZE; ++i)
    {
        sobel[i][0] = 0;
        sobel[i][COLSIZE - 1] = 0;
    }
    for (int j = 0; j < COLSIZE; ++j)
    {
        sobel[0][j] = 0;
        sobel[ROWSIZE - 1][j] = 0;
    }
    for (int i = 0; i < ROWSIZE; ++i)
    {
        for (int j = 0; j < COLSIZE; ++j)
        {
            if (sobel[i][j] > 160) //upper threshold value is 160
            {
                sobel[i][j] = 1;
                hysteresis[i][j] = 2;
            }
            else
            {
                if (sobel[i][j] > 60) //lower threshold is 60
                    hysteresis[i][j] = 1;
                else
                    hysteresis[i][j] = 0;
                sobel[i][j] = 0;
            }
        }
    }

    for (int i = 0; i < ROWSIZE; ++i)
    {
        for (int j = 0; j < COLSIZE; ++j)
        {
            if (hysteresis[i][j] != 0 && hysteresis[i][j] != 1 && hysteresis[i][j] != 2) //threshold value is 180
            {
                cout << hysteresis[i][j];
            }
        }
    }

    for (int i = 0; i < ROWSIZE; ++i)
        for (int j = 0; j < COLSIZE; ++j)
        {
            if (hysteresis[i][j] == 2)
                edgeFill(i, j);
        }

    for (int i = 0; i < ROWSIZE; ++i)
        for (int j = 0; j < COLSIZE; ++j)
        {
            if (hysteresis[i][j] == 3)
                hysteresis[i][j] = 2;
            else
                hysteresis[i][j] = 0;
        }

    //creates imagem.ppm with edges
    string** arr = new string * [ROWSIZE];
    for (int i = 0; i < ROWSIZE; ++i)
        arr[i] = new string[COLSIZE];

    string s1;
    ofstream file1;
    file1.open("image1.ppm");
    file1 << "P3" << " " << length << " " << width << " " << 1 << "\n";
    for (int i = 0; i <= colcount; i++)
    {
        for (int j = 0; j <= rowcount; j++)
        {
            s1 = "";
            s1 += to_string(nms[i][j]);
            s1 += " ";
            s1 += to_string(nms[i][j]);
            s1 += " ";
            s1 += to_string(nms[i][j]);
            s1.append(" ");
            arr[i][j] = s1;
        }
    }

    for (int i = 0; i < COLSIZE; i++)
    {
        for (int j = 0; j < ROWSIZE; j++)
        {
            file1 << arr[j][i];
        }
        file1 << "\n";
    }

    file1.close();

    file1.open("image2.ppm");
    file1 << "P3" << " " << length << " " << width << " " << 1 << "\n";
    for (int i = 0; i <= colcount; i++)
    {
        for (int j = 0; j <= rowcount; j++)
        {
            s1 = "";
            s1 += to_string(hysteresis[i][j]);
            s1 += " ";
            s1 += to_string(hysteresis[i][j]);
            s1 += " ";
            s1 += to_string(hysteresis[i][j]);
            s1.append(" ");
            arr[i][j] = s1;
        }
    }

    for (int i = 0; i < COLSIZE; i++)
    {
        for (int j = 0; j < ROWSIZE; j++)
        {
            file1 << arr[j][i];
        }
        file1 << "\n";
    }

    file1.close();

    int final;
    file1.open("imagef.ppm");
    file1 << "P3" << " " << length << " " << width << " " << 1 << "\n";
    for (int i = 0; i <= colcount; i++)
    {
        for (int j = 0; j <= rowcount; j++)
        {
            if (hysteresis[i][j] == 2 && nms[i][j] == 1)
                final = 1;
            else
                final = 0;
            s1 = "";
            s1 += to_string(final);
            s1 += " ";
            s1 += to_string(final);
            s1 += " ";
            s1 += to_string(final);
            s1.append(" ");
            arr[i][j] = s1;
        }
    }

    for (int i = 0; i < COLSIZE; i++)
    {
        for (int j = 0; j < ROWSIZE; j++)
        {
            file1 << arr[j][i];
        }
        file1 << "\n";
    }

    file1.close();

    Drawing d;
    d.run(direction, arr, points);
}

int main()
{
    part2();
}
