//
// Created by nathnel on 4/15/18.
//
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
#ifndef INNOVATE_FPGA_LANES_WINDOW_H
#define INNOVATE_FPGA_LANES_WINDOW_H
using namespace cv;
using namespace std;

class Window {
public:
    int tolerance, m, y2, y1;
    int x;
    int mean_x;

    Window(int y1_, int y2_, int x_, int m_=100, int tolerance_ = 50) ;
    vector<Point> coordinates();
    vector<int> pixels_in(vector<Mat> nonzero, int x = 1000000);
};


#endif //INNOVATE_FPGA_LANES_WINDOW_H
