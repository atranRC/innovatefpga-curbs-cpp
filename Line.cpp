//
// Created by nathnel on 4/15/18.
//

#include <utility>
#include <stdexcept>
#include <valarray>
#include "Line.h"
#include "Matrix.h"

Line::Line(vector<int> x, vector<int> y, int h_, int w_) {
    coefficients.resize(5);
    process_points(move(x), move(y));

}

void Line::process_points(vector<int> x_, vector<int> y_) {
    float enough_points = y_.size() > 0 and max(y_) - min(y_) > h * .625;
    if (enough_points || coefficients.size() == 0) {
        fit(x_, y_);
    }

}


void Line::fit(vector<int> x_, vector<int> y_) {
    vector<int> fitted = polyfit(y_, x_, 2);
    coefficients.push_back(fitted);
}

int Line::radiusOfCurvature() {
    double ym_per_pix = 20.0 / 720;
    double xm_per_pix = 1.7 / 700;
    vector<vector<double>> points = getPoints();
    vector<double> xs, ys;
    for ( auto i : points ) {
        for ( int j = 0; j < i.size(); j++ ) {
            if (j == 1) ys.push_back(i[j]);
            else if (j == 0) xs.push_back(i[j]);
        }
    }
    valarray<double> yms = copyValVec(ys);
    auto fit_cr = polyfit(copyValVec(yms * ym_per_pix), copyValVec(copyValVec(xs) * xm_per_pix), 2);
    return (int) (pow((1 + pow(2 * fit_cr[0] * 720 * ym_per_pix + fit_cr[1], 2)), 1.5) / abs(2 * fit_cr[0]));
}

vector<vector<double>> Line::getPoints() {
    vector<double> y = linspace(0, h - 1, h);
    vector<double> current_fit = averagedFit();
    valarray<double> y_arr = copyValVec(y);
    valarray<double> part_2 = (current_fit[1] * y_arr) + current_fit[2];
    valarray<double> part_1 = (current_fit[0] * y_arr * y_arr);
    vector<valarray<double>> stacked(2);
    stacked[0] = part_1 + part_2;
    stacked[1] = y_arr;
    vector<vector<double>> stacked_ret;
    stacked_ret[0] = copyValVec(stacked[0]);
    stacked_ret[1] = copyValVec(stacked[1]);
    return transpose(stacked_ret);
}

// https://stackoverflow.com/questions/27028226/python-linspace-in-c
template<typename T>
vector<double> Line::linspace(T start_in, T end_in, int num_in) {
    std::vector<double> linspaced;
    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for ( int i = 0; i < num - 1; ++i ) {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end);
    return linspaced;
}


vector<double> Line::averagedFit() {
    int max_size;
    for ( auto i : coefficients )
        if (max_size < i.size())
            max_size = i.size();


    vector<double> means(max_size);

    for ( int i = 0; i < coefficients.size(); i++ ) {
        for ( int j = 0; j < coefficients[i].size(); j++ ) {
            means[i] += coefficients[i][j];
        }
    }
    for ( int i = 0; i < means.size(); i++ ) {
        means[i] = (means[i] / coefficients.size());
    }
    return means;
}

valarray<double> Line::copyValVec(vector<double> input) {
    int size = static_cast<int>(input.size());
    valarray<double> copied(size);
    for ( int i = 0; i < size; i++ ) {
        copied[i] = input[i];
    }
    return copied;
}

valarray<int> Line::copyValVec(vector<int> input) {
    int size = static_cast<int>(input.size());
    valarray<int> copied(size);
    for ( int i = 0; i < size; i++ ) {
        copied[i] = input[i];
    }
    return copied;
}

vector<double> Line::copyValVec(valarray<double> input) {
    int size = static_cast<int>(input.size());
    vector<double> copied(size);
    for ( int i = 0; i < size; i++ ) {
        copied[i] = input[i];
    }
    return copied;
}

vector<int> Line::copyValVec(valarray<int> input) {
    int size = static_cast<int>(input.size());
    vector<int> copied(size);
    for ( int i = 0; i < size; i++ ) {
        copied[i] = input[i];
    }
    return copied;
}

template<typename T>
vector<vector<T>> Line::transpose(vector<vector<T>> &b) {
    if (b.size() == 0)
        throw range_error("Size == 0");

    vector<vector<T> > trans_vec(b[0].size(), vector<T>());

    for ( int i = 0; i < b.size(); i++ ) {
        for ( int j = 0; j < b[i].size(); j++ ) {
            trans_vec[j].push_back(b[i][j]);
        }
    }

    return trans_vec;
}

int Line::max(vector<int> lis) {
    int max_ = lis[0];
    for ( auto &i : lis ) {
        if (max_ < i) max_ = i;
    }
    return max_;
}

double Line::max(vector<double> lis) {
    double max_ = lis[0];
    for ( auto &i : lis ) {
        if (max_ < i) max_ = i;
    }
    return max_;
}

int Line::min(vector<int> lis) {
    int min_ = lis[0];
    for ( auto &i : lis ) {
        if (min_ > i) min_ = i;
    }
    return min_;
}

// http://www.vilipetek.com/2013/10/07/polynomial-fitting-in-c-using-boost/
template<typename T>
vector<T> Line::polyfit(const vector<T> &oX, const vector<T> &oY, int nDegree) {
    if (oX.size() != oY.size())
        throw std::invalid_argument("X and Y vector sizes do not match");

    nDegree++;

    size_t nCount = oX.size();
    matrix<T> oXMatrix(nCount, nDegree);
    matrix<T> oYMatrix(nCount, 1);

    for ( size_t i = 0; i < nCount; i++ ) {
        oYMatrix(i, 0) = oY[i];
    }

    for ( size_t nRow = 0; nRow < nCount; nRow++ ) {
        T nVal = 1.0f;
        for ( int nCol = 0; nCol < nDegree; nCol++ ) {
            oXMatrix(nRow, nCol) = nVal;
            nVal *= oX[nRow];
        }
    }

    matrix<T> oXtMatrix(oXMatrix.transpose());
    matrix<T> oXtXMatrix(oXtMatrix * oXMatrix);
    matrix<T> oXtYMatrix(oXtMatrix * oYMatrix);

    Givens<T> oGivens;
    oGivens.Decompose(oXtXMatrix);
    matrix<T> oCoeff = oGivens.Solve(oXtYMatrix);
    return oCoeff.data();
}

int Line::cameraDistance() {
    vector<vector<double>> points = getPoints();
    double xm_per_pix = 2.7 / 700;
    vector<double> ys;
    for ( int i = 0; i < points.size(); i++ ) {
        for ( int j = 0; j < points[i].size(); j++ ) {
            if (j == 0) ys.push_back(points[i][i]);
        }
    }
    double x = points[(int) max(ys)][0];
    return (int) abs((w / 2 - x) * xm_per_pix);
}

