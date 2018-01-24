#include "pfHome.h"
#include "pfLmpDrv.h"

typedef std::chrono::high_resolution_clock Clock;

int main(int argc, char* argv[]) {
	auto t1 = Clock::now();

	pfHome* pfdrv = new pfHome(argc, argv);
	pfdrv->run(argc, argv);
	delete pfdrv;
	auto t2 = Clock::now();
	std::cout << "Delta t2-t1: " 
              << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
              << " seconds" << std::endl;
	return (0);	
}