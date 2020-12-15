//
// Created by nathnel on 4/13/18.
//
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <cmath>
#include <vector>
#include "Window.h"
#include "Line.h"

using namespace cv;
using namespace std;
#ifndef INNOVATE_FPGA_LANES_LANEFINDER_H
#define INNOVATE_FPGA_LANES_LANEFINDER_H


class LaneFinder {
private:
    Mat arctan2(Mat &a, Mat &b);
public:
    int width;
    int height;
    int win_n;
    Line left, right;
    vector<Window> l_windows;
    vector<Window> r_windows;

    LaneFinder(Mat first_frame, int num_windows=9);
    Mat get_edges(Mat frame);
    int initialize_lines(Mat frame);
    Mat gradient_abs_value_mask(Mat& image, int sobel_kernel=3, char axis='x', vector<int> threshold = {0, 255});
    Mat gradient_magnitude_mask(Mat& image, int sobel_kernel=3, vector<int> threshold = {0, 255});
    Mat gradient_direction_mask(Mat& image, int sobel_kernel=3, vector<double> threshold = {0.2, 1.5});
    Mat color_threshold_mask(Mat& image, vector<int> threshold = {0, 255});
    Mat process(Mat& image, bool drawLane = true, bool drawStatistics = true);
    vector<int> scanFrameWithWindows(Mat& frame, vector<Window> windows);
    Mat flattenPerspective(Mat& image, Mat& outImage);

    vector<int> slice(vector<int> vector, int start, int end);

    int argmax(vector<int> inp);
};


#endif //INNOVATE_FPGA_LANES_LANEFINDER_H
