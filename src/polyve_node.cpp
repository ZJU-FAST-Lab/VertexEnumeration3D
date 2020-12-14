#include "polyve/polyve.h"

int main(int argc, char **argv)
{
    ros::init(argc, argv, "polyve_node");
    ros::NodeHandle nh_;

    Config config;
    ros::NodeHandle nh_priv("~");

    config.loadParameters(nh_priv);
    PolyVe polyVe(config, nh_);

    if (config.testRate > 0.0)
    {
        ros::Rate lr(config.testRate);
        while (ros::ok())
        {
            polyVe.conductVE();
            ros::spinOnce();
            lr.sleep();
        }
    }
    else
    {
        ros::spin();
    }

    return 0;
}
