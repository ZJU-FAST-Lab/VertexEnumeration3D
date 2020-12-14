#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <Eigen/Eigen>

#include <ros/ros.h>

#include <string>
#include <vector>

struct Config
{
    std::string triggerTopic;
    std::string meshTopic;
    std::string edgeTopic;
    std::string vertexTopic;
    std::string interiorTopic;
    int randomTries;
    double randomScale;
    int redundantTryH;
    double testRate;

    inline void loadParameters(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TriggerTopic", triggerTopic);
        nh_priv.getParam("MeshTopic", meshTopic);
        nh_priv.getParam("EdgeTopic", edgeTopic);
        nh_priv.getParam("VertexTopic", vertexTopic);
        nh_priv.getParam("InteriorTopic", interiorTopic);
        nh_priv.getParam("RandomTries", randomTries);
        nh_priv.getParam("RandomScale", randomScale);
        nh_priv.getParam("RedundantTryH", redundantTryH);
        nh_priv.getParam("TestRate", testRate);
        return;
    }
};

#endif