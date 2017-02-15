#include "defs.h"
#include "database.h"
#include "utils.h"
#include "points_utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace pr;

const char* banner[]={
  "Display test",
  "Browse through the laser scans using A and D.",
  0
};


int main (int argc, char** argv) {
    printBanner(banner);
    
    ifstream inputFile("../loop-detector-2d.txt");
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

    map<int,Vector2fVector> m2 = database.getScanDrawings();
    
    
    RGBImage shown_image(400,400);
    
    char key=0;
    int t = 0;
    const char ESC_key=27;
    while (key!=ESC_key) {
        cout << sequence_list[t] << " " << t  << endl;
        shown_image=cv::Vec3b(255,255,255);

        Vector2fVector new_points = m2.find(sequence_list[t])->second;
        Vector2fVector display_points = database.pointsToDisplay(new_points);
        drawPoints(shown_image,display_points, cv::Scalar(255,0,0),1);
                cv::imshow("display_test", shown_image);

        key=cv::waitKey(0);
            switch(key) {
                case 'd': t = min (t+1, scan_count-1); break;
                case 'a' : t = max(t-1,0);break;
                default: ;
                }
        }

}
