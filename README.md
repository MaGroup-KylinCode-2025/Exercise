基于C++实现的分子结构表示类，用于量子化学计算。

## 项目结构
molecular-class/
├── cpp/ # C++源代码目录
│ ├── CMakeLists.txt # 项目构建配置
│ ├── main.cpp # 主测试程序
│ ├── molecular.h # Molecular类声明
│ └── molecular.cpp # Molecular类实现
├── h2o.txt # 水分子测试数据
└── README.md # 项目说明文档

text

## 编译运行
```bash
cd cpp/
mkdir build && cd build
cmake .. && make
./mol_class ../../h2o.txt
