#include "Stag.h"
#include "opencv2/opencv.hpp"
#include <opencv2/ximgproc.hpp>

using namespace cv;
using namespace cv::ximgproc;
using namespace std;

int main(int argc, char** argv) {

    char* img = "00000.png";
    if (argc > 1)
    {
        img = argv[1];
    }

  cv::Mat image = cv::imread(img, cv::IMREAD_GRAYSCALE);
  if (image.empty())
  {
      cout << "empty" << endl;
      return 1;
  }


  Stag stag(23, 7, true);

  stag.detectMarkers(image);
  stag.logResults("");
  return 0;
}
