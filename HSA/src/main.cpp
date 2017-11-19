#include <cstdio>
#include "HSA.h"
//#define DEBUG
//#define GLOBAL
//HSA::HSA hsa;




int main(int argc, char** argv) {
	HSA::HSA *hsa = new HSA::HSA();
	

#ifdef DEBUG
    if (hsa->BuildCost("C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\cost2.tsv") && hsa->BuildTreeTS("C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\fun.alm", "C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\fun.alm")) {
#ifdef GLOBAL
        hsa->outputGResult("C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\fun.almg");
#else
		std::vector<std::tuple<std::string, std::string>> subTreeRootsTS;
		subTreeRootsTS = hsa->outputLResult("C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\fun.alm",4);
		delete hsa;
		kOutputAndGeneratePValue(subTreeRootsTS, "C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\cost2.tsv", "C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\fun.alm","C:\\Users\\ATPs\\Documents\\GitHub\\VS\\HSA\\x64\\Release\\fun.alm");

#endif // GLOBAL

    }
#else
    if (argc > 4) {
		int testNum = 100;
        if (hsa->BuildCost(argv[3]) && hsa->BuildTreeTS(argv[1], argv[2])) {
            if (*(argv[4]) == Setting::GlobalC) {
				if (argc > 5) {
					sscanf(argv[5], "%ut", &testNum);
					hsa->TEST_CASE = testNum;
				}
                hsa->outputGResult(argv[1]);
            } else {
				if (argc > 6) {
					sscanf(argv[6], "%ut", &testNum);
					hsa->TEST_CASE = testNum;
				}
                int num = 1;
                sscanf(argv[5], "%d", &num);
				std::vector<std::tuple<std::string, std::string>> subTreeRootsTS;
				subTreeRootsTS = hsa->outputLResult(argv[1], num);//vector of subTreeRoots ids
				kOutputAndGeneratePValue(subTreeRootsTS, argv[3], argv[1], argv[2], testNum);
            }
        }
    } else {
		fprintf(stderr, "Usage: HSA <TreeS file path> <TreeT file path> <Cost file path> <g>\n");
		fprintf(stderr, "or: HSA <TreeS file path> <TreeT file path> <Cost file path> <l> <n>\n");
		fprintf(stderr, "or: HSA <TreeS file path> <TreeT file path> <Cost file path> <g> <testNum>\n");
		fprintf(stderr, "or: HSA <TreeS file path> <TreeT file path> <Cost file path> <l> <n> <testNum>\n");
		fprintf(stderr, "n: number of local alignments to output\ntestNum: number of test for calculate p-value (default:100, max:10000,min:2)\n");
		
        return -1;
    }
#endif // !DEBUG
    return 0;
}
