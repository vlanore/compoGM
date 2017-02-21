#include "model.hpp"


// ######################################
//                INTERFACES
// ######################################
class Message {
public:
    virtual std::string getMessage() = 0;
};


// ######################################
//                INSTANCES
// ######################################
class I_Int : public Instance {
public:
    int val;

    I_Int(int val) : val(val) {}

    void hello() override {
        printf("I_Int(%d)", val);
    }
};


// ######################################
//                  TEST
// ######################################
int main() {
    Assembly mymodel;

    mymodel.instantiate<Instance>("I1");
    mymodel.instantiate<I_Int>("I2", 7);
    mymodel.instantiate<Array<I_Int> >("I3", 5, [](int i) { return I_Int(i); });
    mymodel.instantiate<I_Int>("I4", 8);

    mymodel.set<I_Int, int>("I2", &I_Int::val, 37);

    mymodel.print_all();

    return 0;
}
