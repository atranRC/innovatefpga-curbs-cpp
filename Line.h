//
// Created by nathnel on 4/15/18.
//
#include <deque>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;
#ifndef INNOVATE_FPGA_LANES_LINE_H
#define INNOVATE_FPGA_LANES_LINE_H


class Line {
    int h;
    int w;
    int frame_impact = 0;
    deque<vector<int>> coefficients;

    Line(vector<int>, vector<int>, int, int);

    void process_points(vector<int> x_, vector<int> y_);

    int max(vector<int> lis);

    int min(vector<int> lis);

    vector<vector<double>> getPoints();

    int radiusOfCurvature();

    int cameraDistance();

    vector<double> averagedFit();

    void fit(vector<int> x_, vector<int> y_);

    template<typename T>
    vector<T> polyfit(const vector<T> &oX, const vector<T> &oY, int nDegree);

    template<typename T>
    vector<double> linspace(T start_in, T end_in, int num_in);

    valarray<double> copyValVec(vector<double> input);
    valarray<int> copyValVec(vector<int> input);
    vector<double> copyValVec(valarray<double> input);
    vector<int> copyValVec(valarray<int> input);

    template<typename T>
    vector<vector<T>> transpose(vector<vector<T>> &b);

    double max(vector<double> lis);
};


#endif //INNOVATE_FPGA_LANES_LINE_H
