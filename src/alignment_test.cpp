#include "defs.h"
#include "database.h"
#include "utils.h"
#include "points_utils.h"
#include "solver.h"
#include "distance_map_correspondence_finder.h"



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace pr;

const char* banner[]={
    "Alignement test.",
    "Align points from different laser scans using distance map correspondances.",
    "Spacebar for 1 alignement round, enter for complete alignement, x to display the actual best transformation.",
    0
};
int main (int argc, char** argv) {
    printBanner(banner);
    
    
    
    if ((argc == 2)||(argc > 3)){
        cerr<< "You have provided "<< argc - 1<<" input parameter! Provide nothing to compare the 2 default scans or the 2 ids of the scans you want to compare. The default scans will be compared in this run."<<endl;}
    
    ifstream inputFile("../loop-detector-2d.txt");
    string odomFile("../trajectoryLaser.g2o");
    string line;
    int threshold = 5;
    ScanDatabase database(threshold);
    vector<int> sequence_list;
    int scan_count = 0;
    
    
    
    //Get image matrices from each scan line
    while (getline(inputFile, line))
    {
        int seq_number = database.extractData(line);
        sequence_list.push_back(seq_number);
        database.pointsFromScan(seq_number);
        scan_count++;
    }
    
    database.extractOdometry(odomFile);
    
    
    
    
    
    
    map<int, vector<float>> m1 = database.getScanList();
    map<int,Vector2fVector> m2 = database.getScanDrawings();
    map<int, Eigen::Vector3f> m3 = database.getOdometryList();
    
    Eigen::Vector3f lastPose = m3.begin()->second; //initialize at first known pose
    
    for (int iD : sequence_list){
        if (m3.find(iD)!=m3.end()){
            lastPose = m3.find(iD)->second;
        }
        else {
            m3.insert(pair<int,Eigen::Vector3f>(iD, lastPose));
        }
    }
    
    //Script to print all the iDs with their initial guess
    /*
     for(std::map<int, Eigen::Vector3f>::iterator it = m3.begin(); it != m3.end(); ++it) {
     
     vector<int>::iterator itera;
     itera=find(sequence_list.begin(),sequence_list.end(), it->first);
     std::cout<< it->first<<" " << (itera - sequence_list.begin())<< " " << (it->second).transpose()<<std::endl;
     
     }*/
    
    
    
    //From 30690 to 33709...TOT 2413 scan
    
    int id_1 = 33033;
    int id_2 = 31318;
    
    if (argc == 3){
        
        int x;
        
        istringstream ss1(argv[1]);
        if (!(ss1 >> x)){
            cerr << "Invalid number " << argv[1] << endl; }
        else{
            id_1 = x; }
        
        istringstream ss2(argv[2]);
        if (!(ss2 >> x)){
            cerr << "Invalid number " << argv[1] << endl; }
        else {
            id_2 = x; }
        
        if ((id_1 < sequence_list[0] ) && (id_1 < sequence_list.size()) && (id_2 < sequence_list[0] ) && (id_2 < sequence_list.size()) ) {
            id_1 = sequence_list[id_1];
            id_2 = sequence_list[id_2];
        }
        
    }
    
    Vector2fVector reference_points = m2.find(id_1)->second;
    Vector2fVector scan_points = m2.find(id_2)->second;
    
    
    
    
    Vector2fVector display_reference = database.pointsToDisplay(reference_points);
    
    
    Solver solver;
    
    Eigen::Vector3f ref_pose = (m3.find(id_1))->second;
    Eigen::Vector3f actual_pose = (m3.find(id_2))->second;
    
    
    Eigen::Matrix3f guess = vec2mat(ref_pose).inverse() * vec2mat(actual_pose);
    
    
    solver.init(reference_points, scan_points,guess);
    
    cout<<"Initial Guess Vector: "<<mat2vec(guess).transpose()<<endl;
    
    Vector2fVector corrected_points;
    for (Eigen::Vector2f vec: scan_points){
        vec = guess.block<2,2>(0,0)*vec + guess.block<2,1>(0,2);
        corrected_points.push_back(vec);
    }
    
    
    
    
    DistanceMapCorrespondenceFinder correspondence_finder;
    int rows = 500;
    int cols = 500;
    float max_distance=15;
    correspondence_finder.init(display_reference,
                               rows,
                               cols,
                               max_distance);
    correspondence_finder.compute(corrected_points);
    
    
    
    
    
    RGBImage shown_image(rows,cols);
    
    bool iterate = true;
    char key=0;
    const char ESC_key=27;
    
    while (key!=ESC_key) {
        shown_image=cv::Vec3b(255,255,255);
        
        Vector2fVector display_corrected = database.pointsToDisplay(corrected_points);
        Vector2fVector display_scan = database.pointsToDisplay(scan_points);
        
        correspondence_finder.compute(display_corrected);
        
        
        drawPoints(shown_image,display_reference, cv::Scalar(255,0,0),1);
        
        drawPoints(shown_image,display_scan, cv::Scalar(255,255,0),1);
        
        drawPoints(shown_image,display_corrected, cv::Scalar(0,255,0),1);
        
        
        
        drawCorrespondences(shown_image,
                            display_reference,
                            display_corrected,
                            correspondence_finder.correspondences(), cv::Scalar(0,0,255));
        
        float chi = solver.computeError(correspondence_finder.correspondences(),scan_points,reference_points );
        cout<<"Error: " << chi << " Mean Error: "<<chi/correspondence_finder.correspondences().size()<< " Matched points ratio: "<< 100*correspondence_finder.correspondences().size()/scan_points.size()<<"%"<< endl;
        
        cv::imshow("alignment_test", shown_image);
        key=cv::waitKey(0);
        switch(key) {
                
            case 13 : {//enter
                correspondence_finder.compute(display_corrected);
                
                while (iterate){
                    Vector2fVector new_points;
                    Eigen::Matrix3f mat = solver.getX();
                    for (Eigen::Vector2f point : scan_points){
                        Eigen::Vector2f new_point;
                        new_point = mat.block<2,2>(0,0)*point + mat.block<2,1>(0,2);
                        new_points.push_back(new_point);
                        
                    }
                    corrected_points = new_points;
                    Vector2fVector display_corrected = database.pointsToDisplay(corrected_points);
                    correspondence_finder.compute(display_corrected);
                    
                    iterate = solver.oneRound(correspondence_finder.correspondences(),false);
                    
                    if (!iterate){
                        

                        
}
                    
                }
                
            }break;
                
                
            case ' ' : {
                
                correspondence_finder.compute(display_corrected);
                
                solver.oneRound(correspondence_finder.correspondences(),false);
                Vector2fVector new_points;
                Eigen::Matrix3f mat = solver.getX();
                for (Eigen::Vector2f point : scan_points){
                    Eigen::Vector2f new_point;
                    new_point = mat.block<2,2>(0,0)*point + mat.block<2,1>(0,2);
                    new_points.push_back(new_point);
                    
                }
                corrected_points = new_points;
                
                
                
            }break;
                
                
            case 'x' : {
                cout<< solver.getX()<< endl;
                cout << "Resulting transformation: " << mat2vec(solver.getX()).transpose()<<endl;
            }break;
        }
    }
    
}
