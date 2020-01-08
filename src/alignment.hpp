/*******************************************************
 * 
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Dai Bo (bdairobot@gmail.com)
 * file name: alignment.hpp
 * function: used to align odometry and imu to get scale and other initial paramters
 *******************************************************/ 
#pragma once
#include "integration_base.h"

class Alignment
{
public:
    Alignment(){

    }

    void solveGyroscopeBias();
    bool LinearAlignment();
private:


}