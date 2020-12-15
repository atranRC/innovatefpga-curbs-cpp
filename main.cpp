#include <iostream>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/opencv.hpp>
#include "LaneFinder.h"
using namespace cv;
using namespace std;

int main() {

    Mat image = imread("/home/nathnel/PycharmProjects/innovatefpgacv/curbs/images/main.jpg");
    LaneFinder im(image);
    return 0;
}