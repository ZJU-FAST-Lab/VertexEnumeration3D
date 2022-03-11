#ifndef POLYVE_H
#define POLYVE_H

#include "polyve/config.hpp"
#include "geo_utils/quickhull.hpp"
#include "geo_utils/geo_utils.hpp"

#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include <ros/ros.h>
#include <ros/console.h>
#include <std_msgs/Empty.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

class Visualization
{
public:
    Visualization(Config &conf, ros::NodeHandle &nh_);

    Config config;
    ros::NodeHandle nh;

    ros::Publisher meshPub;
    ros::Publisher edgePub;
    ros::Publisher verticesPub;
    ros::Publisher interiorPub;

    void visualizeMesh(const Eigen::Matrix3Xd &mesh);
    void visualizeVertices(const Eigen::Matrix3Xd &vertices);
    void visualizeInterior(const Eigen::Vector3d &interior);
};

class PolyVe
{
public:
    PolyVe(Config &conf, ros::NodeHandle &nh_);
    ~PolyVe();

    Config config;
    ros::NodeHandle nh;

    ros::Subscriber triggerSub;
    void triggerCallBack(const std_msgs::Empty::ConstPtr &msg);

    void conductVE();
    Visualization visualization;
};
#endif