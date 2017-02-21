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
class I_Youpi : public Instance, public Message {
public:
    std::string name;
    int size;

    Message* Out;

    I_Youpi(std::string name, int size) : name(name), size(size) {}

    void hello() override {
        printf("Child class: %s(%d).\n", name.c_str(), size);
        if (Out != nullptr)
            printf("<Neighbour: %s>\n", Out->getMessage().c_str());
    }

    std::string getMessage() override {
        return "Hello!";
    }
};


// ######################################
//                  TEST
// ######################################
int main() {
    Assembly mymodel;

    typedef UseProvide<I_Youpi, Message> UseMessage;

    mymodel.instantiate<Instance>("I1");
    mymodel.instantiate<I_Youpi>("I2", "alice", 7);
    mymodel.instantiate<Array<I_Youpi> >("IArray", 5, [](int i) { return I_Youpi("legion", i); });
    mymodel.instantiate<I_Youpi>("I3", "bob", 8);

    mymodel.set<I_Youpi, int>("I2", &I_Youpi::size, 37);

    mymodel.connect<UseMessage>("I2", "I3", &I_Youpi::Out);

    mymodel.print_all();

    return 0;
}
