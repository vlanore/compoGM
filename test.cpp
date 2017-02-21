#include <cstdio>
#include <map>
#include <memory>
#include <string>


// ######################################
//                INSTANCES
// ######################################
class Instance {
public:
    std::shared_ptr<Instance> Out;
    virtual void hello() { printf("Base class.\n"); }
};

class SuperInstance : public Instance {
public:
    std::string name;
    int size;
    SuperInstance(std::string name, int size) : name(name), size(size) {}
    void hello() override {
        printf("Child class: %s(%d).\n", name.c_str(), size);
        neighbour();
    }

    void neighbour() {
        if (Out != nullptr) {
            printf("Neighbour: ");
            Out->hello();
        }
    }
};


// ######################################
//                INTERFACES
// ######################################
class Interface {
public:

};


// ######################################
//                  MODEL
// ######################################
class Model {
public:
    std::map<std::string, std::shared_ptr<Instance> > instances;

    void print_all() {
        for (auto i : instances) {
            printf("%s: ", i.first.c_str());
            i.second->hello();
        }
    }

    template <class T, class... Args>
    void instantiate(std::string name, Args &&... args) {
        instances.emplace(name, std::make_shared<T>(std::forward<Args>(args)...));
    }

    void connect(std::string i1, std::string i2) {
        instances[i1]->Out = instances[i2];
    }
};


// ######################################
//                  TEST
// ######################################
int main() {
    Model mymodel;

    mymodel.instantiate<Instance>("I1");
    mymodel.instantiate<SuperInstance>("I2", "alice", 7);
    mymodel.instantiate<SuperInstance>("I3", "bob", 8);
    mymodel.connect("I2", "I3");

    mymodel.print_all();

    return 0;
}
