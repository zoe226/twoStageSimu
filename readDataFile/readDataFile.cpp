// readDataFile.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <vector>
//#include <string>
//using namespace std;
//typedef unsigned char byte;

int main()
{
    // C++ 实现
    std::string filename = "ADCData_Frame_000000_20240314_17_22_19_2876_ShareMemory.bin";

    const size_t bufferSize = 83886080;

    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "无法打开文件" << std::endl;
        return 1;
    }

    
    std::vector<unsigned char> buffer(bufferSize,0);
    file.read(reinterpret_cast<char*>(buffer.data()), bufferSize);
    size_t bytesRead = file.gcount();

    // 按照帧头结构解析出参数


    std::cout << "文件内容：" << std::endl;
    for (size_t i = 270; i < std::min<size_t>(270+10, bytesRead);i++ ) {
        std::cout << static_cast<int>(buffer[i]) << " "; //以整数形式输出unsigned char值
    }
    std::cout << std::endl;

    // 检查是否发生错误或者到达文件结尾
    if (file.bad()) {
        std::cerr << "读取文件时发生错误" << std::endl;
    }
    else if (file.eof()) {
        std::cout << "已到达文件末尾" << std::endl;
    }
    file.close();
    return 0;

    // C 实现
    /*FILE* fid;
    byte *RadarInputDataShareMem = new byte[83886080];

    memset(RadarInputDataShareMem, 0, 83886080 * sizeof(byte));
    string filename = "ADCData_Frame_000000_20240314_17_22_19_2876_ShareMemory.bin";
    fopen_s(&fid, filename.c_str(), "rb");
    if (fid == NULL) {
        std::cout << "failed to open file!\n";
        return 0;
    }
    fread(RadarInputDataShareMem, sizeof(byte), 83886080, fid);
    fclose(fid);

    delete RadarInputDataShareMem;

    std::cout << "Hello World!\n";
    return 0;*/
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
