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

    std::vector<float> yYT;

    for (auto &ii : yY)
    {
        yYT.push_back(1/ii);
    }


    plt::plot(xX, yY);
    plt::show();
}
//-----------------------------------------------------------

std::vector<float> Hamiltonian_dot(float p1, float p2, float q1, float q2)
{

    // constans

    float w1 = std::sqrt(2);
    float w2 = 1;

    float dp1 = -w1 * q1;
    float dp2 = -w2 * q2;
    float dq1 =  w1 * p1;
    float dq2 =  w2 * p2;

    std::vector<float> dot = {dp1, dp2, dq1, dq2};

    return dot;
}

std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> computeHamiltonian()
{

    int size = 100000;
    std::vector<float> dot;
    float dt = 0.001;

    std::vector<float> bCurveP1(size);
    std::vector<float> bCurveP2(size);
    std::vector<float> bCurveQ1(size);
    std::vector<float> bCurveQ2(size);
    std::vector<float> dtime(size);

    float init_p1 = 0.5;
    float init_p2 = 1.0;
    float init_q1 = 0.3;
    float init_q2 = 0.0;

    bCurveP1[0] = init_p1;
    bCurveP2[0] = init_p2;
    bCurveQ1[0] = init_q1;
    bCurveQ2[0] = init_q2;
    dtime[0]  = 0.0;

    for (int ii = 0; ii < size - 1; ii++)
    {
        dot = Hamiltonian_dot(bCurveP1[ii], bCurveP2[ii], bCurveQ1[ii], bCurveQ2[ii]);
        bCurveP1[ii + 1] = bCurveP1[ii] + dot[0] * dt;
        bCurveP2[ii + 1] = bCurveP2[ii] + dot[1] * dt;
        bCurveQ1[ii + 1] = bCurveQ1[ii] + dot[2] * dt;
        bCurveQ2[ii + 1] = bCurveQ2[ii] + dot[3] * dt;
        dtime [ii+1] = dt * ii;
    }

    return std::make_tuple(bCurveP1, bCurveP2, bCurveQ1, bCurveQ2, dtime);
}

//-----------------------------------------------------------

int main()
{
    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> hamiltonianCurve_w_t = computeHamiltonian();

    std::tuple<std::vector<float>, std::vector<float>, std::vector<float>> hamiltonianCurve = {std::get<0>(hamiltonianCurve_w_t), std::get<2>(hamiltonianCurve_w_t), std::get<3>(hamiltonianCurve_w_t)};
    
    std::tuple<std::vector<float>, std::vector<float>> rossler_xt = {std::get<0>(hamiltonianCurve_w_t), std::get<2>(hamiltonianCurve_w_t)};

    plot3D(hamiltonianCurve);

    plot2D(rossler_xt);
}