//// g++ rossler_attractor.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//-----------------------------------------------------------

void plot3Dexample()
{

    std::vector<float> xX;
    std::vector<float> yY;
    std::vector<float> zZ;
    float theta;
    float r;
    float z_inc = 4.0 / 99.0;
    float theta_inc = (8.0 * M_PI) / 99.0;

    for (float i = 0; i < 100; i += 1)
    {
        theta = -4.0 * M_PI + theta_inc * i;
        zZ.push_back(-2.0 + z_inc * i);
        r = zZ[i] * zZ[i] + 1;
        xX.push_back(r * std::sin(theta));
        yY.push_back(r * std::cos(theta));
    }

    plt::plot3(xX, yY, zZ);
    plt::show();
}

//-----------------------------------------------------------

void plot3D(std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> data)
{

    std::vector<float> xX = std::get<0>(data);
    std::vector<float> yY = std::get<1>(data);
    std::vector<float> zZ = std::get<2>(data);

    plt::plot3(xX, yY, zZ);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("z");
    plt::show();
}

//-----------------------------------------------------------

void plot2D(std::tuple<std::vector<float>, std::vector<float>> data)
{

    std::vector<float> xX = std::get<0>(data);
    std::vector<float> yY = std::get<1>(data);

    plt::plot(xX, yY);
    plt::show();
}
//-----------------------------------------------------------

std::vector<float> Rossler_dot(float x, float y, float z)
{

    // constans

    float a = 15.395;
    float b = 28;
    float R = -1.143;
    float c_2 = -0.714;

    float g = c_2 * x + 0.5 * (R - c_2) * (std::abs(x+1) - std::abs(x-1)); 

    float x_dot = a * (y - x - g);
    float y_dot = x - y + z;
    float z_dot = -b * y;

    std::vector<float> dot = {x_dot, y_dot, z_dot};

    return dot;
}

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> computeRossler()
{

    int size = 900000;
    std::vector<float> dot;
    float dt = 0.0001;

    std::vector<float> bCurveX(size);
    std::vector<float> bCurveY(size);
    std::vector<float> bCurveZ(size);
    std::vector<float> dtime(size);

    float init_x = 0.1;
    float init_y = 0.0;
    float init_z = 0.1;

    bCurveX[0] = init_x;
    bCurveY[0] = init_y;
    bCurveZ[0] = init_z;
    dtime[0]  = 0.0;

    for (int ii = 0; ii < size - 1; ii++)
    {
        dot = Rossler_dot(bCurveX[ii], bCurveY[ii], bCurveZ[ii]);
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
    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> rosslerCurve_w_t = computeRossler();

    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> rosslerCurve = {std::get<0>(rosslerCurve_w_t), std::get<1>(rosslerCurve_w_t), std::get<2>(rosslerCurve_w_t)};
    
    std::tuple<std::vector<float>, std::vector<float>> rossler_xt = {std::get<3>(rosslerCurve_w_t), std::get<0>(rosslerCurve_w_t)};

    plot3D(rosslerCurve);

   // plot2D(rossler_xt);
}