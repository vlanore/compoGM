#include <cstdio>
#include <map>
#include <memory>
#include <string>

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
class SuperInstance : public Instance, public Message {
public:
    std::string name;
    int size;
    Message* Out;
    SuperInstance(std::string name, int size) : name(name), size(size) {}
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
    Model mymodel;

    typedef UseProvide<SuperInstance, Message> UseMessage;

    mymodel.instantiate<Instance>("I1");
    mymodel.instantiate<SuperInstance>("I2", "alice", 7);
    mymodel.instantiate<SuperInstance>("I3", "bob", 8);

    mymodel.set<SuperInstance, int>("I2", &SuperInstance::size, 37);

    mymodel.connect<UseMessage>("I2", "I3", &SuperInstance::Out);

    mymodel.print_all();

    return 0;
}
