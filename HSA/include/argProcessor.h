#include <map>
#include <string>
using namespace std;

class ArgProcessor {
private:
	map<string, string> argVal;
public:
	ArgProcessor(int argc, char* argv[]);


};

ArgProcessor::ArgProcessor(int argc, char* argv[]) {
	for (int i = 1; i < argc; i++) {
		char* arg;
	}
}