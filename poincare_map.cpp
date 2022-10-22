//// g++ rossler_attractor.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;


//-----------------------------------------------------------

void plot2D(std::tuple<std::vector<float>, std::vector<float>> data)
{

    std::vector<float> xX = std::get<0>(data);
    std::vector<float> yY = std::get<1>(data);

    plt::plot(xX, yY);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::show();
}
//-----------------------------------------------------------

std::vector<float> PoincareMap_dot(float x, float y, float theta)
{

    // constans

    float alpha = 1.0;
    float beta = -1.0;
    float gamma = 0.8;
    float k = 0.3;
    float omega = 1.25;

    float x_dot = y;
    float y_dot = -beta * x - k * y - alpha * std::pow(x,3) + gamma * std::cos(theta);
    float theta_dot = omega;

    std::vector<float> dot = {x_dot, y_dot, theta_dot};

    return dot;
}

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> computePoincare()
{

    int size = 4000;
    std::vector<float> dot;
    float dt = 0.01;

    std::vector<float> bCurveX(size);
    std::vector<float> bCurveY(size);
    std::vector<float> bCurveZ(size);
    std::vector<float> dtime(size);

    float init_x = 0.0;
    float init_y = 0.0;
    float init_z = 0.0;

    bCurveX[0] = init_x;
    bCurveY[0] = init_y;
    bCurveZ[0] = init_z;
    dtime[0]  = 0.0;

    for (int ii = 0; ii < size - 1; ii++)
    {
        dot = PoincareMap_dot(bCurveX[ii], bCurveY[ii], bCurveZ[ii]);
        bCurveX[ii + 1] = bCurveX[ii] + dot[0] * dt;
        bCurveY[ii + 1] = bCurveY[ii] + dot[1] * dt;
        bCurveZ[ii + 1] = bCurveZ[ii] + dot[2] * dt;
        dtime [ii+1] = dt * ii;
    }

    return std::make_tuple(bCurveX, bCurveY, bCurveZ, dtime);
}

//-----------------------------------------------------------

int main()
{
    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> map_w_t = computePoincare();

    
    std::tuple<std::vector<float>, std::vector<float>> map_xy = {std::get<0>(map_w_t), std::get<1>(map_w_t)};

    
    plot2D(map_xy);
}