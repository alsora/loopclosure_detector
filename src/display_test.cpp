#include "defs.h"
#include "display.h"
#include "utils.h"
#include "points_utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace pr;

const char* banner[]={
  "display test",
  "represent a laser scan.",
  0
};


int main (int argc, char** argv) {
    printBanner(banner);
    
    //ifstream inputFile("../prova.txt");
    ifstream inputFile("../03-2DLoopDetector/loop-detector-2d.txt");
    ifstream g2oFile;
    string line;
    string line2;
    ScanDatabase database(5);
    vector<int> sequence_list;
    int scan_count = 0;
    
    while (getline(inputFile, line))
    {
        
       int seq_number = database.extractData(line);
        sequence_list.push_back(seq_number);
        database.pointsFromScan(seq_number);
        scan_count++;

        
    }
    //m1 non usata
    map<int, vector<float>> m1 = database.getScanList();
   map<int,Vector2fVector> m2 = database.getScanDrawings();
    
    
    RGBImage shown_image(400,700);
    
    char key=0;
    int t = 0;
    const char ESC_key=27;
    while (key!=ESC_key) {
        cout << sequence_list[t] << " " << t  << endl;
        shown_image=cv::Vec3b(255,255,255);

        Vector2fVector new_points = m2.find(sequence_list[t])->second;
        Vector2fVector display_points = database.pointsToDisplay(new_points);
        drawPoints(shown_image,display_points, cv::Scalar(255,0,0),1);
                cv::imshow("camera_test", shown_image);
        //for (Eigen::Vector2f vec: new_points){
          //  cout<< vec.x() << " " <<vec.y() << endl;
       // }
    key=cv::waitKey(0);
        switch(key) {
            case 'd': t = min (t+1, scan_count-1); break;
            case 'a' : t = max(t-1,0);break;
            default: ;
        }
    }

}
