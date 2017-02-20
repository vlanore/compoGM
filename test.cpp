#include <string>
#include <vector>
#include <memory>
#include <cstdio>

class Instance {
public:
    std::string name;
    Instance(std::string name): name(name) {}
    virtual void hello() { printf("Base class: %s\n", name.c_str()); }
};

class SuperInstance: public Instance {
public:
    SuperInstance(std::string name): Instance(name) {}
    void hello() override { printf("Child class: %s\n", name.c_str()); }
};

class Model {
public:
    std::vector<std::shared_ptr<Instance> > instances;

    void print_all() {
        for (auto i: instances) {
            i->hello();
        }
    }

    template <class T> void instantiate(std::string name) {
        instances.emplace_back(std::make_shared<T>(name));
    }
};

int main () {
    Model mymodel;

    mymodel.instantiate<Instance>("bob");
    mymodel.instantiate<SuperInstance>("alice");

    mymodel.print_all();

    return 0;
}
