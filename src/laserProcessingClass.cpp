// Author of FLOAM: Wang Han 
// Email wh200720041@gmail.com
// Homepage https://wanghan.pro
#include "laserProcessingClass.h"

void LaserProcessingClass::init(lidar::Lidar lidar_param_in){
    
    lidar_param = lidar_param_in;

}

void LaserProcessingClass::featureExtraction(const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_in, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_edge, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_surf){

    std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*pc_in, indices);


    int N_SCANS = lidar_param.num_lines;
    std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr> laserCloudScans;
    for(int i=0;i<N_SCANS;i++){
        laserCloudScans.push_back(pcl::PointCloud<pcl::PointXYZI>::Ptr(new pcl::PointCloud<pcl::PointXYZI>()));
    }

    for (int i = 0; i < (int) pc_in->points.size(); i++)
    {
        int scanID=0;
        double distance = sqrt(pc_in->points[i].x * pc_in->points[i].x + pc_in->points[i].y * pc_in->points[i].y);
        if(distance<lidar_param.min_distance || distance>lidar_param.max_distance)
            continue;
        double angle = atan(pc_in->points[i].z / distance) * 180 / M_PI;
        
        if (N_SCANS == 16)
        {
            scanID = int((angle + 15) / 2 + 0.5);
            if (scanID > (N_SCANS - 1) || scanID < 0)
            {
                continue;
            }
        }
        else if (N_SCANS == 32)
        {
            scanID = int((angle + 92.0/3.0) * 3.0 / 4.0);
            if (scanID > (N_SCANS - 1) || scanID < 0)
            {
                continue;
            }
        }
        else if (N_SCANS == 64)
        {   
            if (angle >= -8.83)
                scanID = int((2 - angle) * 3.0 + 0.5);
            else
                scanID = N_SCANS / 2 + int((-8.83 - angle) * 2.0 + 0.5);

            if (angle > 2 || angle < -24.33 || scanID > 63 || scanID < 0)
            {
                continue;
            }
        }
        else
        {
            printf("wrong scan number\n");
        }
        laserCloudScans[scanID]->push_back(pc_in->points[i]); 

    }

    for(int i = 0; i < N_SCANS; i++){
        if(laserCloudScans[i]->points.size()<131){
            continue;
        }
        
        std::vector<Double2d> cloudCurvature; 
        int total_points = laserCloudScans[i]->points.size()-10;
        for(int j = 5; j < (int)laserCloudScans[i]->points.size() - 5; j++){
            double diffX = laserCloudScans[i]->points[j - 5].x + laserCloudScans[i]->points[j - 4].x + laserCloudScans[i]->points[j - 3].x + laserCloudScans[i]->points[j - 2].x + laserCloudScans[i]->points[j - 1].x - 10 * laserCloudScans[i]->points[j].x + laserCloudScans[i]->points[j + 1].x + laserCloudScans[i]->points[j + 2].x + laserCloudScans[i]->points[j + 3].x + laserCloudScans[i]->points[j + 4].x + laserCloudScans[i]->points[j + 5].x;
            double diffY = laserCloudScans[i]->points[j - 5].y + laserCloudScans[i]->points[j - 4].y + laserCloudScans[i]->points[j - 3].y + laserCloudScans[i]->points[j - 2].y + laserCloudScans[i]->points[j - 1].y - 10 * laserCloudScans[i]->points[j].y + laserCloudScans[i]->points[j + 1].y + laserCloudScans[i]->points[j + 2].y + laserCloudScans[i]->points[j + 3].y + laserCloudScans[i]->points[j + 4].y + laserCloudScans[i]->points[j + 5].y;
            double diffZ = laserCloudScans[i]->points[j - 5].z + laserCloudScans[i]->points[j - 4].z + laserCloudScans[i]->points[j - 3].z + laserCloudScans[i]->points[j - 2].z + laserCloudScans[i]->points[j - 1].z - 10 * laserCloudScans[i]->points[j].z + laserCloudScans[i]->points[j + 1].z + laserCloudScans[i]->points[j + 2].z + laserCloudScans[i]->points[j + 3].z + laserCloudScans[i]->points[j + 4].z + laserCloudScans[i]->points[j + 5].z;
            Double2d distance(j,diffX * diffX + diffY * diffY + diffZ * diffZ);
            cloudCurvature.push_back(distance);

        }
        for(int j=0;j<6;j++){
            int sector_length = (int)(total_points/6);
            int sector_start = sector_length *j;
            int sector_end = sector_length *(j+1)-1;
            if (j==5){
                sector_end = total_points - 1; 
            }
            std::vector<Double2d> subCloudCurvature(cloudCurvature.begin()+sector_start,cloudCurvature.begin()+sector_end); 
            
            featureExtractionFromSector(laserCloudScans[i],subCloudCurvature, pc_out_edge, pc_out_surf);
            
        }

    }

}

void LaserProcessingClass::featureExtraction2(const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_in, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_edge, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_surf) {
    float verticalAngle, horizonAngle, range;
    size_t rowIdn, columnIdn, index, cloudSize;
    pcl::PointXYZI thisPoint;
    int N_SCANS = lidar_param.num_lines;
    cloudSize = pc_in->points.size();
    int Horizon_SCAN = 360.0 / lidar_param.horizontal_angle_resolution;
    auto rangeMat = cv::Mat(64, Horizon_SCAN, CV_32F, cv::Scalar::all(FLT_MAX));
    cv::normalize(rangeMat,rangeMat);
    auto XYZMat = cv::Mat_<cv::Vec3f>(64, Horizon_SCAN);

    for (size_t i = 0; i < cloudSize; ++i) {
        thisPoint.x = pc_in->points[i].x;
        thisPoint.y = pc_in->points[i].y;
        thisPoint.z = pc_in->points[i].z;

        double distance = sqrt(pc_in->points[i].x * pc_in->points[i].x + pc_in->points[i].y * pc_in->points[i].y);
        if (distance < lidar_param.min_distance || distance > lidar_param.max_distance)
            continue;
        double angle = atan(pc_in->points[i].z / distance) * 180 / M_PI;

        if (N_SCANS == 16) {
            rowIdn = int((angle + 15.0) / 2.0 + 0.5);
            if (rowIdn > (N_SCANS - 1) || rowIdn < 0) {
                continue;
            }
        } else if (N_SCANS == 32) {
            rowIdn = int((angle + 92.0 / 3.0) * 3.0 / 4.0);
            if (rowIdn > (N_SCANS - 1) || rowIdn < 0) {
                continue;
            }
        } else if (N_SCANS == 64) {
            if (angle >= -8.83)
                rowIdn = int((2.0 - angle) * 3.0 + 0.5);
            else
                rowIdn = N_SCANS / 2 + int((-8.83 - angle) * 2.0 + 0.5);

            if (angle > 2 || angle < -24.33 || rowIdn > 63 || rowIdn < 0) {
                continue;
            }
        } else {
            printf("wrong scan number\n");
        }
//        rowIdn = N_SCANS - 1 - rowIdn; //从上往下计算，对应深度图。
        horizonAngle = atan2(thisPoint.y, thisPoint.x) * 180 / M_PI;
        columnIdn = -round((horizonAngle - 180.0) / lidar_param.horizontal_angle_resolution);
        if (columnIdn >= Horizon_SCAN)
            columnIdn -= Horizon_SCAN;
        if (columnIdn < 0 || columnIdn >= Horizon_SCAN)
            continue;


        range = sqrt(thisPoint.x * thisPoint.x + thisPoint.y * thisPoint.y + thisPoint.z * thisPoint.z);
        rangeMat.at<float>(rowIdn, columnIdn) = range;
        XYZMat.at<cv::Vec3f>(rowIdn, columnIdn)[0] = thisPoint.x * 1000;
        XYZMat.at<cv::Vec3f>(rowIdn, columnIdn)[1] = -thisPoint.z * 1000;
        XYZMat.at<cv::Vec3f>(rowIdn, columnIdn)[2] = thisPoint.y * 1000;
    }
    int length_threshold = 20;
    int distance_threshold = 3;
    int canny_th1 = 4;
    int canny_th2 = 4;
    int canny_aperture_size = 3;
    bool do_merge = true;

    cv::Mat tmpMat = rangeMat.clone();
    cv::Mat tmpMat2 = rangeMat.clone();
    cv::resize(tmpMat, tmpMat, cv::Size(Horizon_SCAN, 51*4), 0, 0, cv::INTER_LINEAR);
    cv::Ptr<cv::ximgproc::FastLineDetector> fld = cv::ximgproc::createFastLineDetector(length_threshold,
                                                                                       distance_threshold, canny_th1,
                                                                                       canny_th2, canny_aperture_size,
                                                                                       do_merge);
    std::vector<cv::Vec4f> lines_fld;
    tmpMat.convertTo(tmpMat, CV_8U, 255. / 80.);
    fld->detect(tmpMat, lines_fld);
//    std::cout << "lines_fld.size() =" << lines_fld.size() << std::endl;
    cv::imshow("1",tmpMat);
    fld->drawSegments(tmpMat, lines_fld);
//
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_all_line(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr line_info_cloud(new pcl::PointCloud<pcl::PointXYZI>);
    int line_cnt = 0;
    for (auto &l: lines_fld) {
        bool flag = true;
        cv::LineIterator lit(rangeMat, cv::Point(l[0], l[1]/4), cv::Point(l[2], l[3]/4));
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>);
        cloud->reserve(100);
        for (int j = 0; j < lit.count; ++j, ++lit) {
            int col = lit.pos().x;
            int row = lit.pos().y;
            pcl::PointXYZI pt;
            pt.x = XYZMat[row][col][0] / 1000.; //mm -> m
            pt.y = XYZMat[row][col][2] / 1000.;
            pt.z = -XYZMat[row][col][1] / 1000.;
            if (pt.z == 0. || isnan(pt.z))
            {
                flag= false;
                break;
            }
//                pt.b = color_im.ptr<uchar>(row)[col * 3];
//                pt.g = color_im.ptr<uchar>(row)[col * 3 + 1];
//                pt.r = color_im.ptr<uchar>(row)[col * 3 + 2];
//                pt.label = num_of_line;
            cloud->emplace_back(pt);
        }
        if (cloud->size() < 5||!flag)
            continue;

        //-----------------------------拟合直线-----------------------------
        pcl::SampleConsensusModelLine<pcl::PointXYZI>::Ptr model_line(
                new pcl::SampleConsensusModelLine<pcl::PointXYZI>(cloud));
        pcl::RandomSampleConsensus<pcl::PointXYZI> ransac(model_line);
        ransac.setDistanceThreshold(0.2);    //内点到模型的最大距离
        ransac.setMaxIterations(1000);        //最大迭代次数
        ransac.computeModel();                //直线拟合
        //--------------------------根据索引提取内点------------------------
        std::vector<int> inliers;
        ransac.getInliers(inliers);
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_line(new pcl::PointCloud<pcl::PointXYZI>);
        pcl::copyPointCloud<pcl::PointXYZI>(*cloud, inliers, *cloud_line);
        *cloud_all_line += *cloud_line;

        Eigen::VectorXf coef;
        ransac.getModelCoefficients(coef);
        pcl::PointXYZI line_direction_info;
//        line_direction_info.x = coef[3]; //nx
//        line_direction_info.y = coef[4]; //ny
//        line_direction_info.z = coef[5]; //nz
//        line_direction_info.label = num_of_line;
//        line_info_cloud->push_back(line_direction_info);

//        num_of_line++;
        line_cnt++;
    }
    pcl::PointXYZI line_num;
    line_num.x = static_cast<float>(line_cnt);
    pcl::PointCloud<pcl::PointXYZI>::Ptr line_num_cloud(new pcl::PointCloud<pcl::PointXYZI>);
    line_num_cloud->push_back(line_num);
    *pc_out_edge = *cloud_all_line;
    std::cout << "pc_out_edge.size() =" << pc_out_edge->size() << std::endl;

    struct OrganizedImage3D {
        const cv::Mat_<cv::Vec3f> &cloud;

        //note: ahc::PlaneFitter assumes mm as unit!!!
        OrganizedImage3D(const cv::Mat_<cv::Vec3f> &c) : cloud(c) {}

        inline int width() const { return cloud.cols; }

        inline int height() const { return cloud.rows; }

        inline bool get(const int row, const int col, double &x, double &y, double &z) const {
            const cv::Vec3f &p = cloud.at<cv::Vec3f>(row, col);
            x = p[0];
            y = p[1];
            z = p[2];
            return isnan(z) == 0; //return false if current depth is NaN
        }
    };
    ahc::PlaneFitter<OrganizedImage3D> pf;
    pf.minSupport = 100;
    pf.windowWidth = 4;
    pf.windowHeight = 4;
    pf.doRefine = true;

    cv::Mat seg(XYZMat.rows, XYZMat.cols, CV_8UC3);
    std::vector<std::vector<int>> vSeg;
    OrganizedImage3D Ixyz(XYZMat);
    pf.run(&Ixyz, &vSeg, &seg);
//    cv::resize(seg, seg, cv::Size(Horizon_SCAN, 200), 0, 0, cv::INTER_LINEAR);
    cv::imshow("seg",seg);
//        cv::resize(rangeMat, rangeMat, cv::Size(1800, 160), 0, 0, cv::INTER_LINEAR);
//    cv::resize(tmpMat2, tmpMat2, cv::Size(Horizon_SCAN, 200), 0, 0, cv::INTER_LINEAR);
//
//    cv::imshow("depth",tmpMat2/20.0);
    cv::waitKey(1);
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_all_plane(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr plane_info_cloud(new pcl::PointCloud<pcl::PointXYZI>);
    int plane_cnt = 0;

    for (auto idx_plane = 0; idx_plane < vSeg.size(); idx_plane++) {
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZI>);
        cloud->reserve(1000);
        for (auto idx_idx = 0; idx_idx < vSeg[idx_plane].size(); idx_idx++) {
            pcl::PointXYZI pt;
            int pt_idx = vSeg[idx_plane].at(idx_idx);
            int row = pt_idx / XYZMat.cols;
            int col = pt_idx % XYZMat.cols;
            int gap_tmp = 1;
            if(row%gap_tmp==0 && col%gap_tmp==0)
            {
            pt.x = XYZMat[row][col][0] / 1000.; //mm -> m
            pt.y = XYZMat[row][col][2] / 1000.;
            pt.z = -XYZMat[row][col][1] / 1000.;

//                pt.b = seg.ptr<uchar>(row)[col * 3];
//                pt.g = seg.ptr<uchar>(row)[col * 3 + 1];
//                pt.r = seg.ptr<uchar>(row)[col * 3 + 2];
//                pt.label = num_of_plane;
                cloud->emplace_back(pt);
            }
        }

        if (cloud->size() < 3) {
            continue;
        }

        //--------------------------RANSAC拟合平面--------------------------
        pcl::SampleConsensusModelPlane<pcl::PointXYZI>::Ptr model_plane(
                new pcl::SampleConsensusModelPlane<pcl::PointXYZI>(cloud));
        pcl::RandomSampleConsensus<pcl::PointXYZI> ransac(model_plane);
        ransac.setDistanceThreshold(0.2);    //设置距离阈值，与平面距离小于0.1的点作为内点
        ransac.computeModel();                //执行模型估计
        //-------------------------根据索引提取内点--------------------------
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_plane(new pcl::PointCloud<pcl::PointXYZI>);
        std::vector<int> inliers;                //存储内点索引的容器
        ransac.getInliers(inliers);            //提取内点索引
        pcl::copyPointCloud<pcl::PointXYZI>(*cloud, inliers, *cloud_plane);

        *cloud_all_plane += *cloud_plane;

        //----------------------------输出模型参数---------------------------
//        Eigen::VectorXf coefficient;
//        ransac.getModelCoefficients(coefficient);
//        if(coefficient[3]<0)
//        {
//            coefficient = -coefficient;
//        }
//        pcl::PointXYZI plane_info;
//        plane_info.x = coefficient[0];
//        plane_info.y = coefficient[1];
//        plane_info.z = coefficient[2];
//        plane_info.data[3] = coefficient[3];
//        plane_info.rgb = cloud_plane->size();
//        plane_info.label = num_of_plane;
//        plane_info_cloud->push_back(plane_info);
//
//        num_of_plane++;
//        plane_cnt++;
//    }
//    pcl::PointXYZI plane_num;
//    plane_num.x = static_cast<float>(plane_cnt);
//    pcl::PointCloud<pcl::PointXYZI>::Ptr plane_num_cloud(new pcl::PointCloud<pcl::PointXYZI>);
//    plane_num_cloud->push_back(plane_num);
        *pc_out_surf = *cloud_all_plane;//num of plane + plane info + plane points
    }
}
void LaserProcessingClass::featureExtractionFromSector(const pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_in, std::vector<Double2d>& cloudCurvature, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_edge, pcl::PointCloud<pcl::PointXYZI>::Ptr& pc_out_surf){
    std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*pc_in, indices);
    std::sort(cloudCurvature.begin(), cloudCurvature.end(), [](const Double2d & a, const Double2d & b)
    { 
        return a.value < b.value; 
    });


    int largestPickedNum = 0;
    std::vector<int> picked_points;
    int point_info_count =0;
    for (int i = cloudCurvature.size()-1; i >= 0; i--)
    {
        int ind = cloudCurvature[i].id; 
        if(std::find(picked_points.begin(), picked_points.end(), ind)==picked_points.end()){
            if(cloudCurvature[i].value <= 0.1){
                break;
            }
            
            largestPickedNum++;
            picked_points.push_back(ind);
            
            if (largestPickedNum <= 20){
                pc_out_edge->push_back(pc_in->points[ind]);
                point_info_count++;
            }else{
                break;
            }

            for(int k=1;k<=5;k++){
                double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k - 1].x;
                double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k - 1].y;
                double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k - 1].z;
                if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
                    break;
                }
                picked_points.push_back(ind+k);
            }
            for(int k=-1;k>=-5;k--){
                double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k + 1].x;
                double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k + 1].y;
                double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k + 1].z;
                if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
                    break;
                }
                picked_points.push_back(ind+k);
            }

        }
    }

    //find flat points
    // point_info_count =0;
    // int smallestPickedNum = 0;
    
    // for (int i = 0; i <= (int)cloudCurvature.size()-1; i++)
    // {
    //     int ind = cloudCurvature[i].id; 

    //     if( std::find(picked_points.begin(), picked_points.end(), ind)==picked_points.end()){
    //         if(cloudCurvature[i].value > 0.1){
    //             //ROS_WARN("extracted feature not qualified, please check lidar");
    //             break;
    //         }
    //         smallestPickedNum++;
    //         picked_points.push_back(ind);
            
    //         if(smallestPickedNum <= 4){
    //             //find all points
    //             pc_surf_flat->push_back(pc_in->points[ind]);
    //             pc_surf_lessFlat->push_back(pc_in->points[ind]);
    //             point_info_count++;
    //         }
    //         else{
    //             break;
    //         }

    //         for(int k=1;k<=5;k++){
    //             double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k - 1].x;
    //             double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k - 1].y;
    //             double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k - 1].z;
    //             if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
    //                 break;
    //             }
    //             picked_points.push_back(ind+k);
    //         }
    //         for(int k=-1;k>=-5;k--){
    //             double diffX = pc_in->points[ind + k].x - pc_in->points[ind + k + 1].x;
    //             double diffY = pc_in->points[ind + k].y - pc_in->points[ind + k + 1].y;
    //             double diffZ = pc_in->points[ind + k].z - pc_in->points[ind + k + 1].z;
    //             if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05){
    //                 break;
    //             }
    //             picked_points.push_back(ind+k);
    //         }

    //     }
    // }
    
    for (int i = 0; i <= (int)cloudCurvature.size()-1; i++)
    {
        int ind = cloudCurvature[i].id; 
        if( std::find(picked_points.begin(), picked_points.end(), ind)==picked_points.end())
        {
            pc_out_surf->push_back(pc_in->points[ind]);
        }
    }
    


}
LaserProcessingClass::LaserProcessingClass(){
    
}

Double2d::Double2d(int id_in, double value_in){
    id = id_in;
    value =value_in;
};

PointsInfo::PointsInfo(int layer_in, double time_in){
    layer = layer_in;
    time = time_in;
};
