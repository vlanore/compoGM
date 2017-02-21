#include <cstdio>
#include <map>
#include <memory>
#include <string>


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
class Instance {
public:
    virtual void hello() { printf("Base class.\n"); }
};

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

    template <class I, typename V>
    void set(std::string name, V I::* member, V value) {
        dynamic_cast<I*>(instances[name].get())->*member = value;
    }

    template <class O, class D>
    void connect(std::string i1, std::string i2, D* O::* member) {
        set<O,D*>(i1, member, dynamic_cast<D*>(instances[i2].get()));
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
    mymodel.set<SuperInstance, int>("I2", &SuperInstance::size, 37);
    mymodel.connect<SuperInstance, Message>("I2", "I3", &SuperInstance::Out);

    mymodel.print_all();

    return 0;
}
