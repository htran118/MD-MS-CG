#include <iostream>
#include <string>

class BigObject;
using namespace std;
int main(int argc, char *argv[]) {
    BigObject one;

    string str = "this is a stirng";
    cout<<str;
}

class BigObject {
public:
    BigObject(){
        cout<<"big object \n";
    }
    virtual ~BigObject(){}      
};
