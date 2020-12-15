//
// Created by nathnel on 4/13/18.
//

#include "LaneFinder.h"
#include "Line.h"

LaneFinder::LaneFinder(cv::Mat first_frame, int num_windows) {
    width = first_frame.cols;
    height = first_frame.rows;
    win_n = num_windows;
    initialize_lines(first_frame);
//    process(first_frame, true, true);
}

int LaneFinder::initialize_lines(cv::Mat frame) {
    Mat edges = get_edges(frame);
    Mat flat_edges;
    Mat unwrap_matrix = flattenPerspective(edges, flat_edges);

    Mat histogram_part = flat_edges(Range((int) height / 2, flat_edges.rows), Range(0, flat_edges.cols));
    vector<int> histogram;
    reduce(histogram_part, histogram, 0, CV_REDUCE_SUM);
    vector<Mat> nonzero;
    findNonZero(flat_edges, nonzero);
    vector<int> l_indices;
    vector<int> r_indices;
    int window_height = height / win_n;
    for ( int i = 0; i < win_n; i++ ) {
        Window l_window(height - (i + 1) * window_height, height - i * window_height,
                        !l_windows.empty() ? l_windows[l_windows.size() - 1].x : argmax(
                                slice(histogram, 0, width / 2)));
        Window r_window(height - (i + 1) * window_height, height - i * window_height,
                        !l_windows.empty() ? r_windows[r_windows.size() - 1].x :
                        argmax(slice(histogram, width / 2, width)) + width / 2);
        vector<int> l_nonzeros = l_window.pixels_in(nonzero);
        vector<int> r_nonzeros = l_window.pixels_in(nonzero);
        for (auto &item : l_nonzeros)
            l_indices.push_back(item);
        for (auto &item : r_nonzeros)
            r_indices.push_back(item);
        l_windows.push_back(l_window);
        r_windows.push_back(r_window);
    }
    Mat x = nonzero[0], y = nonzero[1];

    left = Line(y(l_indices), x[l_indices], height, width);
    right = Line(y[r_indices], x[r_indices], height, width);
     // TODO: Appending
    //    l_indices = np.append(l_indices, l_window.pixels_in(nonzero), axis=0)
    //    r_indices = np.append(r_indices, r_window.pixels_in(nonzero), axis=0)
    //    self.l_windows.append(l_window)
    //    self.r_windows.append(r_window)
//    self.left = Line(x=nonzero[1][l_indices], y=nonzero[0][l_indices], h=self.h, w=self.w)
//    self.right = Line(x=nonzero[1][r_indices], y=nonzero[0][r_indices], h=self.h, w=self.w)


}

Mat LaneFinder::get_edges(Mat frame) {
    Mat hls;
    cvtColor(frame, hls, CV_RGB2HLS);
    vector<Mat> hlsChannels(3);
    split(hls, hlsChannels);
    Mat s_channel = hlsChannels[2];
//    imshow("Image2", s_channel);
//    imshow("Image", mask);

    Mat gradient_x, gradient_y, magnitude, direction;
    gradient_x = gradient_abs_value_mask(s_channel, 3, 'x', {20, 100});
    gradient_y = gradient_abs_value_mask(s_channel, 3, 'y', {20, 100});
    magnitude = gradient_magnitude_mask(s_channel, 3, {20, 100});
    direction = gradient_direction_mask(s_channel, 3, {0.2, 1.5});
    Mat gradient_mask = Mat::zeros(s_channel.rows, s_channel.cols, CV_8UC1); // NOLINT
    gradient_mask = (gradient_x == 255) & (gradient_y == 255) | (magnitude == 255) & (direction == 255);

    Mat color_mask = color_threshold_mask(s_channel, {55, 100});

    Mat mask;
    mask = (gradient_mask == 255) | (color_mask == 255);
//    imshow("Hello", mask);
//    waitKey(0);
    return mask;
}

Mat LaneFinder::gradient_abs_value_mask(Mat &image, int sobel_kernel, char axis, vector<int> threshold) {
    Mat sobeled, sobel;
    if (axis == 'x') {
        Sobel(image, sobeled, CV_64F, 1, 0, sobel_kernel);
        sobel = abs(sobeled);
    } else if (axis == 'y') {
        Sobel(image, sobeled, CV_64F, 0, 1, sobel_kernel);
        sobel = abs(sobeled);
    }
    sobel = sobel * 255 / max(sobel, 0.);
    Mat mask = Mat::zeros(sobel.rows, sobel.cols, CV_8UC1); // NOLINT
//    inRange(sobel, Scalar(threshold[0], threshold[0], threshold[0]), Scalar(threshold[1], threshold[1], threshold[1]),
//            mask);
    mask = sobel <= threshold[1] & sobel >= threshold[0];
    return mask;
}

Mat LaneFinder::gradient_magnitude_mask(Mat &image, int sobel_kernel, vector<int> threshold) {
    Mat sobel_x, sobel_y, mag, mask;
    Sobel(image, sobel_x, CV_64F, 1, 0, sobel_kernel);
    Sobel(image, sobel_y, CV_64F, 0, 1, sobel_kernel);

    magnitude(sobel_x, sobel_y, mag);
//    imshow("mag", arctan2(sobel_x, sobel_y));
    mag = mag * 255 / max(mag, 0.);
//    imshow("dir", mag);
//    waitKey(0);
    mask = Mat::zeros(mag.rows, mag.cols, CV_8UC1); // NOLINT
//    inRange(mag, threshold[0], threshold[1], mask);
    mask = mag <= threshold[1] & mag >= threshold[0];
    return mask;
}


Mat LaneFinder::gradient_direction_mask(Mat &image, int sobel_kernel, vector<double> threshold) {
    Mat sobel_x, sobel_y, dir, mask;
    Sobel(image, sobel_x, CV_64F, 1, 0, sobel_kernel);
    Sobel(image, sobel_y, CV_64F, 0, 1, sobel_kernel);


    dir = arctan2(sobel_x, sobel_y);
    mask = Mat::zeros(dir.rows, dir.cols, CV_8UC1); // NOLINT
//    inRange(dir, threshold[0], threshold[1], mask);
//    mask /= 255;
    mask = dir <= threshold[1] & dir >= threshold[0];
    return mask;
}

Mat LaneFinder::arctan2(Mat &x_edge, Mat &y_edge) {
    Mat angles = Mat::zeros(x_edge.rows, x_edge.cols, CV_16S);
    Mat mags = Mat::zeros(x_edge.rows, x_edge.cols, CV_16S);
//    for ( int i = 0; i < angles.rows; i++ ) {
//        for ( int j = 0; j < angles.cols; j++ ) {
//            angles.at<float>(i, j) = fastAtan2(y_edge.at<int>(i, j), x_edge.at<int>(i, j));
//        }
//    }
    cartToPolar(x_edge, y_edge, mags, angles, false);
    return angles;
}

Mat LaneFinder::color_threshold_mask(Mat &image, vector<int> threshold) {
    Mat mask = Mat::zeros(image.rows, image.cols, CV_8UC1);
    mask = (image > threshold[0]) & (image <= threshold[1]);
    return mask;
}

Mat LaneFinder::process(Mat &image, bool drawLane, bool drawStatistics) {
    Mat edges = get_edges(image);
    Mat flat_edges;
    Mat unwrap_matrix = flattenPerspective(edges, flat_edges);
    vector<int> ls = scanFrameWithWindows(flat_edges, l_windows);
//    (l_x, l_y) = self.scan_frame_with_windows(flat_edges, self.l_windows)
//    self.left.process_points(l_x, l_y)
//            (r_x, r_y) = self.scan_frame_with_windows(flat_edges, self.r_windows)
//    self.right.process_points(r_x, r_y)
//    imshow("Flattened", flat_edges);
//    waitKey(0);
//    destroyAllWindows();
}

Mat LaneFinder::flattenPerspective(Mat &image, Mat &outImage) {
    int h = image.rows, w = image.cols;
    vector<Point2f> source = {Point2f(w / 2 - 76, h * .625), Point2f(w / 2 + 76, h * .625), Point2f(-100, h),
                              Point2f(w + 100, h)};
    vector<Point2f> destination = {Point2f(100, 0), Point2f(w - 100, 0), Point2f(100, h), Point2f(w - 100, h)};
    Mat transform_matrix = getPerspectiveTransform(source, destination);
    Mat unwarp_matrix = getPerspectiveTransform(destination, source);
    warpPerspective(image, outImage, transform_matrix, Size(w, h));
    return unwarp_matrix;
}

vector<int> LaneFinder::scanFrameWithWindows(Mat &frame, vector<Window> windows) {
    vector<int> indices;
    Mat nonzero;
    findNonZero(frame, nonzero);
//    imshow("Non", nonzero);
//    waitKey(0);
    int window_x;
    cout << " Found non zero " << window_x << endl;
    for ( auto &window : windows ) {
        vector<int> pixs = window.pixels_in(nonzero, window_x);
        for ( auto &pix : pixs )
            indices.push_back(pix);
        window_x = window.mean_x;
    }
    cout << window_x;
//    return {nonzero.at(indices, 1), nonzero.at(0, 0)};
//    return ;

//    indices = np.empty([0], dtype=np.int)
//    nonzero = frame.nonzero()
//    window_x = None
//    for window in windows:
//    indices = np.append(indices, window.pixels_in(nonzero, window_x), axis=0)
//    window_x = window.mean_x
//    return (nonzero[1][indices], nonzero[0][indices])
}

vector<int> LaneFinder::slice(vector<int> vec, int start, int end) {
    vector<int> copied = vec;
    vector<int>(copied.begin() + start, copied.begin() + end).swap(copied);
    return copied;
}

int LaneFinder::argmax(vector<int> inp) {
    if (inp.size() <= 0)
        throw range_error("Size <= 0");
    int id = 0;
    for ( int i = 0; i < inp.size(); i++ ) {
        if (inp[id] < inp[i]) id = i;
    }
    return id;
}
