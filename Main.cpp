#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
////////////////////////////////////////////////// ÈÇ 1ñ
int strInInt(const char* str, int size) {
    int res = 0;
    for (int i = 0; i < size; i++) {
        res *= 256;
        res += (int)str[size - 1 - i];
    }
    return res;
}

int getWidth(int n) {
    return n + ((4 - n % 4) % 4);
}

class Image {
public:
    int indentation;
    int width;
    int height;
    int pixSize;
    int size;
    char* data;
    char*** mas;

    Image(const std::string& fName) {
        std::ifstream is(fName, std::ios::binary);
        if (!(is.is_open())) {
            std::cout << "File not found." << std::endl;
            return;
        }
        is.seekg(0, std::ios::end);
        size = is.tellg();
        data = new char[size];
        is.seekg(0, std::ios::beg);
        is.read(data, size);
        is.close();

        indentation = strInInt(data + 10, 4);
        width = strInInt(data + 18, 4);
        height = strInInt(data + 22, 4);
        pixSize = strInInt(data + 28, 2) / 8;
        std::cout << "indentation = " << indentation << std::endl
            << "pixSize = " << pixSize << std::endl
            << "width = " << width << std::endl
            << "height = " << height << std::endl
            << "size = " << size << std::endl << std::endl;

        if (getWidth(width) * height * pixSize + indentation != size) {
            std::cout << "ERROR " << std::endl;
        }
        ///////////////////////////
        mas = new char** [this->height];
        for (int i = 0; i < height; i++)mas[i] = new char* [width];
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                mas[i][j] = new char[3];
            }
        }
        for (int i = 0; i < height; ++i) {
            for (int j = 0; j < width; ++j) {
                int curIdx = indentation + (getWidth(width) * i + j) * pixSize;
                mas[i][j][0] = data[curIdx];
                mas[i][j][1] = data[curIdx + 1];
                mas[i][j][2] = data[curIdx + 2];
            }
        }


    }

    Image(const Image& im) {
        indentation = im.indentation;
        width = im.width;
        height = im.height;
        pixSize = im.pixSize;
        size = im.size;
        data = new char[size];
        for (int i = 0; i < size; ++i) {
            data[i] = im.data[i];
        }
        //mas = im.mas;
        getmas(im);
        mas = im.mas;
    }

    Image& operator=(const Image& im) {
        indentation = im.indentation;
        width = im.width;
        height = im.height;
        pixSize = im.pixSize;
        size = im.size;
        for (int i = 0; i < size; ++i) {
            data[i] = im.data[i];
        }
        return *this;
    }

    ~Image() {
        delete[] data;
    }
    
    void getmas(const Image& bmp) {
        /*char*** mas0 = new char** [this->height];
        // for (int i = 0; i < height; i++)mas0[i] = new char*[width];
         for (int i = 0; i < height; ++i) {
             for (int j = 0; j < width; ++j) {
                 mas[i][j] = new char [3];
             }
         } */
        for (int i = 0; i < bmp.height; ++i) {
            for (int j = 0; j < bmp.width; ++j) {
                int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
                bmp.mas[i][j][0] = bmp.data[curIdx];
                bmp.mas[i][j][1] = bmp.data[curIdx + 1];
                bmp.mas[i][j][2] = bmp.data[curIdx + 2];
            }
        }

    }
};



void getstring(const Image& bmp) {
    /* for (int i = 0; i < height; ++i) {
         for (int j = 0; j < im.width; ++j) {
             int curIdx = im.indentation + (getWidth(im.width) * i + j) * im.pixSize;
             im.data[curIdx] = 0;
             im.data[curIdx + 1] = 0;
             im.data[curIdx + 2] = 0;
         }
     }*/
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
            bmp.data[curIdx] = bmp.mas[i][j][0];
            bmp.data[curIdx + 1] = bmp.mas[i][j][1];
            bmp.data[curIdx + 2] = bmp.mas[i][j][2];
        }
    }

}

int getWealth(const Image& bmp, int idx) {
    return (int)bmp.data[idx] >= 0 ? (int)bmp.data[idx] : (int)bmp.data[idx] + 256;
}

void task3(const Image& bmp, const std::string& fName, int comp) {
    char *res = new char[bmp.size];
    //char res[bmp.size];
    for (int i = 0; i < bmp.size; ++i) {
        res[i] = 0;
    }
    for (int i = 0; i < bmp.indentation; ++i) {
        res[i] = bmp.data[i];
    }
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize + comp;
            res[curIdx] = bmp.data[curIdx];
        }
    }
    std::ofstream os(fName);
    os.write(res, bmp.size);
    os.close();
    delete res;
}

float mWaiting(const Image& bmp, int comp, int x, int y) {
    float res = 0;
    int height = bmp.height - std::abs(y);
    int width = bmp.width - std::abs(x);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * (i + (y > 0 ? y : 0)) + j + (x > 0 ? x : 0)) * bmp.pixSize + comp;
            res += getWealth(bmp, curIdx);
        }
    }
    return res / height / width;
}

float auto_correl(const Image& bmp, int comp1, int comp2, int x, int y) {
    float mWainitgMain = 0;
    float centerSquare1 = 0;
    float centerSquare2 = 0;
    float mWaiting1 = mWaiting(bmp, comp1, -x, -y);
    float mWaiting2 = mWaiting(bmp, comp2, x, y);
    int height = bmp.height - std::abs(y);
    int width = bmp.width - std::abs(x);
    int m = 0; int n = 0;
    m = y;
    for (int i = 0; i < height-y; ++i) {
        n = x;
        for (int j = 0; j < width-x; ++j) {
            int curIdx1 = bmp.indentation + (getWidth(bmp.width) * (i + (y < 0 ? -y : 0)) + j + (x < 0 ? -x : 0)) * bmp.pixSize + comp1;
            int curIdx2 = bmp.indentation + (getWidth(bmp.width) * (m + (y > 0 ? y : 0)) + n + (x > 0 ? x : 0)) * bmp.pixSize + comp2;
            mWainitgMain += (getWealth(bmp, curIdx1) - mWaiting1) * (getWealth(bmp, curIdx2) - mWaiting2);
            centerSquare1 += (getWealth(bmp, curIdx1) - mWaiting1) * (getWealth(bmp, curIdx1) - mWaiting1);
            centerSquare2 += (getWealth(bmp, curIdx2) - mWaiting2) * (getWealth(bmp, curIdx2) - mWaiting2);
            n++;
        }
        m++;
    }

    mWainitgMain /= height * width;
    centerSquare1 = sqrt(centerSquare1 / (height * width - 1));
    centerSquare2 = sqrt(centerSquare2 / (height * width - 1));
    return mWainitgMain / centerSquare1 / centerSquare2;
}


float correl(const Image& bmp, int comp1, int comp2, int x, int y) {
    float mWainitgMain = 0;
    float centerSquare1 = 0;
    float centerSquare2 = 0;
    float mWaiting1 = mWaiting(bmp, comp1, -x, -y);
    float mWaiting2 = mWaiting(bmp, comp2, x, y);
    int height = bmp.height - std::abs(y);
    int width = bmp.width - std::abs(x);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int curIdx1 = bmp.indentation + (getWidth(bmp.width) * (i + (y < 0 ? -y : 0)) + j + (x < 0 ? -x : 0)) * bmp.pixSize + comp1;
            int curIdx2 = bmp.indentation + (getWidth(bmp.width) * (i + (y > 0 ? y : 0)) + j + (x > 0 ? x : 0)) * bmp.pixSize + comp2;
            mWainitgMain += (getWealth(bmp, curIdx1) - mWaiting1) * (getWealth(bmp, curIdx2) - mWaiting2);
            centerSquare1 += (getWealth(bmp, curIdx1) - mWaiting1) * (getWealth(bmp, curIdx1) - mWaiting1);
            centerSquare2 += (getWealth(bmp, curIdx2) - mWaiting2) * (getWealth(bmp, curIdx2) - mWaiting2);
        }
    }
    mWainitgMain /= height * width;
    centerSquare1 = sqrt(centerSquare1 / (height * width - 1));
    centerSquare2 = sqrt(centerSquare2 / (height * width - 1));
    return mWainitgMain / centerSquare1 / centerSquare2;
}

Image getYCbCr(const Image& bmp) {
    Image YCbCr(bmp);
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
            YCbCr.data[curIdx] = std::round(0.299 * getWealth(bmp, curIdx + 2) + 0.587 * getWealth(bmp, curIdx + 1) + 0.114 * getWealth(bmp, curIdx));
            YCbCr.data[curIdx + 1] = std::round(0.5643 * (getWealth(bmp, curIdx) - getWealth(YCbCr, curIdx)) + 128);
            YCbCr.data[curIdx + 2] = std::round(0.7132 * (getWealth(bmp, curIdx + 2) - getWealth(YCbCr, curIdx)) + 128);
        }
    }

    return YCbCr;
}

void saveAutoCorrel(const Image& bmp, const std::string& fPath, const std::string& comps) {
    for (int k = 0; k < 3; ++k) {
        std::string fName = fPath + "/" + comps[k * 2];
        if (comps[k * 2 + 1] != '_') {
            fName += comps[k * 2 + 1];
        }
        std::ofstream os(fName + ".txt", std::ios::out);
        for (int y = -10; y <=10 ; y=y+5) {
            for (int x = -50; x < 50; x++) {
                os << correl(bmp, k, k, x, y) << std::endl;
            }
            os << ' ';
        }
        os.close();
    }
}

void saveCorrel(const Image& bmp, const std::string& fPath, const std::string& comps) {
    for (int k = 0; k < 3; ++k) {
        std::string fName = fPath + "/" + comps[k * 2];
        if (comps[k * 2 + 1] != '_') {
            fName += comps[k * 2 + 1];
        }
        std::ofstream os(fName + ".csv", std::ios::out);
        for (int y = 0; y < 5; ++y) {
            for (int x = -50; x < 50; ++x) {
                os << correl(bmp, k, k, x, (y - 2) * 5) << ',';
            }
            os << correl(bmp, k, k, 50, (y - 2) * 5) << std::endl;
        }
        os.close();
    }
}

void task6(const Image& bmp, const std::string& fName, int comp) {
    char* res = new char[bmp.size];
    //char res[bmp.size];
    for (int i = 0; i < bmp.size; ++i) {
        res[i] = 0;
    }
    for (int i = 0; i < bmp.indentation; ++i) {
        res[i] = bmp.data[i];
    }
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
            res[curIdx] = bmp.data[curIdx + comp];
            res[curIdx + 1] = bmp.data[curIdx + comp];
            res[curIdx + 2] = bmp.data[curIdx + comp];
        }
    }
    std::ofstream os(fName);
    os.write(res, bmp.size);
    os.close();
    delete res;
}

int wealthRound(int wealth) {
    if (wealth < 0) {
        wealth = 0;
    }
    else if (wealth > 255) {
        wealth = 255;
    }
    return wealth;
}

Image getRGB(const Image& bmp) {
    Image _sourse(bmp);
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
            _sourse.data[curIdx] = wealthRound(std::round(getWealth(bmp, curIdx) + 1.772 * (getWealth(bmp, curIdx + 1) - 128)));
            _sourse.data[curIdx + 1] = wealthRound(std::round(getWealth(bmp, curIdx) - 0.714 * (getWealth(bmp, curIdx + 2) - 128) - 0.334 * (getWealth(bmp, curIdx + 1) - 128)));
            _sourse.data[curIdx + 2] = wealthRound(std::round(getWealth(bmp, curIdx) + 1.402 * (getWealth(bmp, curIdx + 2) - 128)));
        }
    }
    return _sourse;
}

void saveImage(const Image& bmp, const std::string& fName) {
    std::ofstream os(fName);
    os.write(bmp.data, bmp.size);
    os.close();
}

float PSNR(const Image& bmp1, const Image& bmp2, int comp) {
    float noise = 0;
    for (int i = 0; i < bmp1.height; ++i) {
        for (int j = 0; j < bmp1.width; ++j) {
            int curIdx = bmp1.indentation + (getWidth(bmp1.width) * i + j) * bmp1.pixSize + comp;
            noise = noise + std::pow(getWealth(bmp1, curIdx) - getWealth(bmp2, curIdx), 2);
        }
    }
    return 10 * std::log10((float)bmp1.height * bmp1.width * 255 * 255 / noise);
}

Image decimation1(const Image& bmp, int step) {
    Image YCbCr(bmp);
   
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
            //int baseIdx = bmp.indentation + (getWidth(bmp.width) * (i - (i % step)) + j - (j % step)) * bmp.pixSize;
            for (int t = 1; t <= step; t++) {
                int baseIdx = bmp.indentation + (getWidth(bmp.width) * (i - (i % t -1)) + j - (j % t)) * bmp.pixSize;
                YCbCr.data[curIdx] = bmp.data[baseIdx];
                YCbCr.data[curIdx + 1] = bmp.data[baseIdx + 1];
                YCbCr.data[curIdx + 2] = bmp.data[baseIdx + 2];
            }
        }
    }
  
    /*int** masY = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masY[i] = new int[bmp.width];
    int** masYc = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masYc[i] = new int[bmp.width];
    int** masYb = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masYb[i] = new int[bmp.width];


    int** masYc0 = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masYc0[i] = new int[bmp.width];
    int** masYb0 = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masYb0[i] = new int[bmp.width];

    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int baseIdx = bmp.indentation + (getWidth(bmp.width) * (i - (i % step)) + j - (j % step)) * bmp.pixSize;
            masY[i][j] = bmp.data[baseIdx];
            masYc[i][j] = bmp.data[baseIdx + 1];
            masYb[i][j] = bmp.data[baseIdx + 2];
        }
    }*/

    return YCbCr;
}

Image decimation2(const Image& bmp, int step) {
    if (bmp.height % step != 0 || bmp.width % step != 0) {
        std::cout << "ERROR" << std::endl;
    }
    Image YCbCr(bmp);
   
    int** masY = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masY[i] = new int[bmp.width];
    int** masYc = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masYc[i] = new int[bmp.width];
    int** masYb = new int* [bmp.height];
    for (int i = 0; i < bmp.height; i++)masYb[i] = new int[bmp.width];
    
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int baseIdx = bmp.indentation + (getWidth(bmp.width) * (i - (i % step)) + j - (j % step)) * bmp.pixSize;
            masY [i][j] = bmp.data[baseIdx];
            masYc[i][j] = bmp.data[baseIdx + 1];
            masYb[i][j] = bmp.data[baseIdx + 2];
        }
    }

    for (int i = 0; i < bmp.height; i=i+step) {
        for (int j = 0; j < bmp.width; j= j+step) {
            int sum1 = 0;
            int sum2 = 0;
            for (int k = i; k < step; step++) {
                for (int t = j; t < step; step++) {
                    sum1 = sum1 + masYc[k][t];
                    sum2 = sum2 + masYb[k][t];
                }
            }
            sum1 = sum1 / pow(step, 2);
            sum2 = sum2 / pow(step, 2);
            for (int k = i; k < step; step++) {
                for (int t = j; t < step; step++) {
                    masYc[k][t] = sum1;
                    masYb[k][t] = sum1;
                }
            }

        }
    }
 
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
            YCbCr.data[curIdx] = masY[i][j];
            YCbCr.data[curIdx+1] = masYc[i][j];
            YCbCr.data[curIdx+2] = masYb[i][j];
        }
    }
    

    return YCbCr;
}

float getHist(const Image& bmp, const std::string& fName, int comp) {
    int hist[256] = { 0 };
    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize + comp;
            hist[getWealth(bmp, curIdx)]++;
        }
    }
    std::ofstream os(fName, std::ios::out);
    for (int i = 0; i < 255; i++) {
        os << i << '*' << hist[i] << std::endl;
    }

    os.close();
    float H = 0;
    int n = bmp.height * bmp.width;
    for (int i = 0; i < 256; ++i) {
        if (hist[i] > 0) {
            H += (float)hist[i] / n * std::log2((float)hist[i] / n);
        }
    }
    return -H;
}

float getDHist(const Image& bmp, const std::string& fName, int comp, int rule) {
    if (rule > 0 && rule < 5) {
        int** D = new int *[bmp.height - 1];
        for (int i = 0; i < bmp.height - 1; i++) D[i] = new int[bmp.width - 1];

        if (rule > 0 && rule < 4) {
            int stepi = rule > 1 ? 1 : 0;
            int strpj = rule % 2 > 0 ? 1 : 0;
            for (int i = 1; i < bmp.height; ++i) {
                for (int j = 1; j < bmp.width; ++j) {
                    int idx1 = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize + comp;
                    int idx2 = bmp.indentation + (getWidth(bmp.width) * (i - stepi) + j - strpj) * bmp.pixSize + comp;
                    D[i - 1][j - 1] = getWealth(bmp, idx1) - getWealth(bmp, idx2);
                }
            }
        }
        else {
            for (int i = 1; i < bmp.height; ++i) {
                for (int j = 1; j < bmp.width; ++j) {
                    int idx1 = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize + comp;
                    int idx2 = bmp.indentation + (getWidth(bmp.width) * (i - 1) + j) * bmp.pixSize + comp;
                    int idx3 = bmp.indentation + (getWidth(bmp.width) * (i - 1) + j - 1) * bmp.pixSize + comp;
                    int idx4 = bmp.indentation + (getWidth(bmp.width) * (i)+j - 1) * bmp.pixSize + comp;
                    D[i - 1][j - 1] = getWealth(bmp, idx1) - (getWealth(bmp, idx2) + getWealth(bmp, idx3) + getWealth(bmp, idx4)) / 3;
                }
            }
        }
        int hist[511] = { 0 };
        for (int i = 1; i < bmp.height; ++i) {
            for (int j = 1; j < bmp.width; ++j) {
                hist[255 + D[i - 1][j - 1]]++;
            }
        }
        std::ofstream os(fName, std::ios::out);
        for (int i = 0; i < 511; i++) {
           // os << (i - 255) << '*' << hist[i] << std::endl;
            os << hist[i] << std::endl;
        }
        os.close();
        float H = 0;
        int n = (bmp.height - 1) * (bmp.width - 1);
        for (int i = 0; i < 511; ++i) {
            if (hist[i] > 0) {
                H += (float)hist[i] / n * std::log2((float)hist[i] / n);
            }
        }
        for (int count = 0; count < bmp.height - 1; count++) delete[] D[count];
        return -H;
    }
    std::cout << "ERROR" << std::endl;
    return 0;
}

Image ret90 (const Image & bmp) {
    Image new_bmp(bmp);

    for (int i = 0; i < bmp.height; ++i) {
        for (int j = 0; j < bmp.width; ++j) {
            int curIdx = bmp.indentation + (getWidth(bmp.width) * i + j) * bmp.pixSize;
            int new_curIdx = bmp.indentation + (getWidth(bmp.width) * (bmp.width-j) + i) * bmp.pixSize;
            new_bmp.data[curIdx] = bmp.data[new_curIdx];
            new_bmp.data[curIdx + 1] = bmp.data[new_curIdx + 1];
            new_bmp.data[curIdx + 2] = bmp.data[new_curIdx + 2];
        }
    }
    return new_bmp;
}

int main() {
    Image sourse("base/1.bmp");

    /*task3(sourse, "task3/blue.bmp", 0);
    task3(sourse, "task3/green.bmp", 1);
    task3(sourse, "task3/red.bmp", 2);

    std::cout << "r(B, G) = " << correl(sourse, 0, 1, 0, 0) << std::endl;
    std::cout << "r(G, R) = " << correl(sourse, 1, 2, 0, 0) << std::endl;
    std::cout << "r(R, B) = " << correl(sourse, 2, 0, 0, 0) << std::endl << std::endl;
    saveCorrel(sourse, "task4", "G_B_R_");

    saveAutoCorrel(sourse, "task4.1", "G_B_R_");*/
    
    Image YCbCr = getYCbCr(sourse);
   /* std::cout << "r(Y, Cb) = " << correl(YCbCr, 0, 1, 0, 0) << std::endl;
    std::cout << "r(Y, Cr) = " << correl(YCbCr, 0, 2, 0, 0) << std::endl;
    std::cout << "r(Cb, Cr) = " << correl(YCbCr, 1, 2, 0, 0) << std::endl << std::endl;
    saveCorrel(YCbCr, "task5", "Y_CbCr");

    task6(YCbCr, "task6/Y.bmp", 0);
    task6(YCbCr, "task6/Cb.bmp", 1);
    task6(YCbCr, "task6/Cr.bmp", 2);
    Image _sourse = getRGB(YCbCr);
    saveImage(_sourse, "task7/recover.bmp");

   std::cout << "PSNR(B) = " << PSNR(sourse, _sourse, 0) << std::endl;
    std::cout << "PSNR(G) = " << PSNR(sourse, _sourse, 1) << std::endl;
    std::cout << "PSNR(R) = " << PSNR(sourse, _sourse, 2) << std::endl << std::endl;
    */
   
    Image decim = decimation1(YCbCr, 2);
    Image redecim = getRGB(decim);
    //
    //redecim= getYCbCr(redecim);
    std::cout << "decim(1) step 2 PSNR(Cb) = " << PSNR(YCbCr, decim, 1) << std::endl;
    std::cout << "decim(1) step 2 PSNR(Cr) = " << PSNR(YCbCr, decim, 2) << std::endl;
    //std::cout << "decim(1) step 2 PSNR(B) = " << PSNR(sourse, redecim, 0) << std::endl;
   // std::cout << "decim(1) step 2 PSNR(G) = " << PSNR(sourse, redecim, 1) << std::endl;
   // std::cout << "decim(1) step 2 PSNR(R) = " << PSNR(sourse, redecim, 2) << std::endl << std::endl;
    saveImage(redecim, "task10/1.bmp");
    
    Image decim1 = decimation2(YCbCr, 2);
    Image redecim1 = getRGB(decim1);
    std::cout << "decim(2) step 2 PSNR(Cb) = " << PSNR(YCbCr, decim1, 1) << std::endl;
    std::cout << "decim(2) step 2 PSNR(Cr) = " << PSNR(YCbCr, decim1, 2) << std::endl;
   /* std::cout << "decim(2) step 2 PSNR(B) = " << PSNR(sourse, redecim, 0) << std::endl;
    std::cout << "decim(2) step 2 PSNR(G) = " << PSNR(sourse, redecim, 1) << std::endl;
    std::cout << "decim(2) step 2 PSNR(R) = " << PSNR(sourse, redecim, 2) << std::endl << std::endl;
    saveImage(redecim1, "task10/2.bmp");
 
    decim = decimation1(YCbCr, 4);
    redecim = getRGB(decim);
    std::cout << "decim(1) step 4 PSNR(Cb) = " << PSNR(YCbCr, decim, 1) << std::endl;
    std::cout << "decim(1) step 4 PSNR(Cr) = " << PSNR(YCbCr, decim, 2) << std::endl;
    /*std::cout << "decim(1) step 4 PSNR(B) = " << PSNR(sourse, redecim, 0) << std::endl;
    std::cout << "decim(1) step 4 PSNR(G) = " << PSNR(sourse, redecim, 1) << std::endl;
    std::cout << "decim(1) step 4 PSNR(R) = " << PSNR(sourse, redecim, 2) << std::endl << std::endl;*/
   /* saveImage(redecim, "task11/1.bmp");

    decim = decimation2(YCbCr, 4);
    redecim = getRGB(decim);
    std::cout << "decim(2) step 4 PSNR(Cb) = " << PSNR(YCbCr, decim, 1) << std::endl;
    std::cout << "decim(2) step 4 PSNR(Cr) = " << PSNR(YCbCr, decim, 2) << std::endl;
   /* std::cout << "decim(2) step 4 PSNR(B) = " << PSNR(sourse, redecim, 0) << std::endl;
    std::cout << "decim(2) step 4 PSNR(G) = " << PSNR(sourse, redecim, 1) << std::endl;
    std::cout << "decim(2) step 4 PSNR(R) = " << PSNR(sourse, redecim, 2) << std::endl << std::endl;*/
    /*saveImage(redecim, "task11/2.bmp");
    
   
    saveImage(ret90(sourse),"dop/dop.bmp");
    std::cout << "H(B) = " << getHist(sourse, "task12/B.csv", 0) << std::endl;
    std::cout << "H(G) = " << getHist(sourse, "task12/G.csv", 1) << std::endl;
    std::cout << "H(R) = " << getHist(sourse, "task12/R.csv", 2) << std::endl << std::endl;
    std::cout << "H(Y)  = " << getHist(YCbCr, "task12/Y.csv", 0) << std::endl;
    std::cout << "H(Cb) = " << getHist(YCbCr, "task12/Cb.csv", 1) << std::endl;
    std::cout << "H(Cr) = " << getHist(YCbCr, "task12/Cr.csv", 2) << std::endl << std::endl;
    */
    std::cout << "D1H(B) = " << getDHist(sourse, "task15/r1/B.csv", 0, 1) << std::endl;
    std::cout << "D1H(G) = " << getDHist(sourse, "task15/r1/G.csv", 1, 1) << std::endl;
    std::cout << "D1H(R) = " << getDHist(sourse, "task15/r1/R.csv", 2, 1) << std::endl << std::endl;
    std::cout << "D1H(Y)  = " << getDHist(YCbCr, "task15/r1/Y.csv", 0, 1) << std::endl;
    std::cout << "D1H(Cb) = " << getDHist(YCbCr, "task15/r1/Cb.csv", 1, 1) << std::endl;
    std::cout << "D1H(Cr) = " << getDHist(YCbCr, "task15/r1/Cr.csv", 2, 1) << std::endl << std::endl;

    std::cout << "D2H(B) = " << getDHist(sourse, "task15/r2/B.csv", 0, 2) << std::endl;
    std::cout << "D2H(G) = " << getDHist(sourse, "task15/r2/G.csv", 1, 2) << std::endl;
    std::cout << "D2H(R) = " << getDHist(sourse, "task15/r2/R.csv", 2, 2) << std::endl << std::endl;
    std::cout << "D2H(Y)  = " << getDHist(YCbCr, "task15/r2/Y.csv", 0, 2) << std::endl;
    std::cout << "D2H(Cb) = " << getDHist(YCbCr, "task15/r2/Cb.csv", 1, 2) << std::endl;
    std::cout << "D2H(Cr) = " << getDHist(YCbCr, "task15/r2/Cr.csv", 2, 2) << std::endl << std::endl;

    std::cout << "D3H(B) = " << getDHist(sourse, "task15/r3/B.csv", 0, 3) << std::endl;
    std::cout << "D3H(G) = " << getDHist(sourse, "task15/r3/G.csv", 1, 3) << std::endl;
    std::cout << "D3H(R) = " << getDHist(sourse, "task15/r3/R.csv", 2, 3) << std::endl << std::endl;
    std::cout << "D3H(Y)  = " << getDHist(YCbCr, "task15/r3/Y.csv", 0, 3) << std::endl;
    std::cout << "D3H(Cb) = " << getDHist(YCbCr, "task15/r3/Cb.csv", 1, 3) << std::endl;
    std::cout << "D3H(Cr) = " << getDHist(YCbCr, "task15/r3/Cr.csv", 2, 3) << std::endl << std::endl;

    std::cout << "D4H(B) = " << getDHist(sourse, "task15/r4/B.csv", 0, 4) << std::endl;
    std::cout << "D4H(G) = " << getDHist(sourse, "task15/r4/G.csv", 1, 4) << std::endl;
    std::cout << "D4H(R) = " << getDHist(sourse, "task15/r4/R.csv", 2, 4) << std::endl << std::endl;
    std::cout << "D4H(Y)  = " << getDHist(YCbCr, "task15/r4/Y.csv", 0, 4) << std::endl;
    std::cout << "D4H(Cb) = " << getDHist(YCbCr, "task15/r4/Cb.csv", 1, 4) << std::endl; 
    std::cout << "D4H(Cr) = " << getDHist(YCbCr, "task15/r4/Cr.csv", 2, 4) << std::endl;

    return 0;
}
