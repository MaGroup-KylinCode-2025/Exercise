#include "molecular.h"
#include <iostream>

int main(int argc, char *argv[]) {
<<<<<<< HEAD
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " filename" << std::endl;
        std::cerr << "Example: ./mol_class h2o.txt" << std::endl;
        return 1;
    }
    
    const char *filename = argv[1];
    try {
        Molecular molecule(filename);
        std::cout << molecule << std::endl;
        
        // 测试setter方法
        molecule.set_basis("6-31G*");
        std::cout << "After setting basis set:" << std::endl;
        std::cout << molecule << std::endl;
        
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
=======
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " filename" << std::endl;
    return 1;
  }
  /*
  const char *filename = argv[1];
  Molecular molecule(filename);
  std::cout << molecule << std::endl;
  */

  return 0;
}
>>>>>>> upstream/AAA-Source-tyb-branch
