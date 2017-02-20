#include <string>
#include <vector>
#include <memory>
#include <cstdio>


// ######################################
//                INSTANCES
// ######################################
class Instance {
public:
    virtual void hello() { printf("Base class.\n"); }
};

class SuperInstance: public Instance {
public:
    std::string name;
    int size;
    SuperInstance(std::string name, int size): name(name), size(size) {}
    void hello() override { printf("Child class: %s(%d).\n", name.c_str(), size); }
};


// ######################################
//                INTERFACES
// ######################################
class Interface {
public:
    std::weak_ptr<Instance> parent;
};


// ######################################
//                  MODEL
// ######################################
class Model {
public:
    std::vector<std::shared_ptr<Instance> > instances;

    void print_all() {
        for (auto i: instances) {
            i->hello();
        }
    }

    template <class T, class... Args> void instantiate(Args&&... args) {
        instances.emplace_back(std::make_shared<T>(std::forward<Args>(args)...));
    }
};


// ######################################
//                  TEST
// ######################################
int main () {
    Model mymodel;

    mymodel.instantiate<Instance>();
    mymodel.instantiate<SuperInstance>("alice", 7);

    mymodel.print_all();

    return 0;
}
