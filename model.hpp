#include <cstdio>
#include <map>
#include <memory>
#include <string>


// ######################################
//                INSTANCES
// ######################################
class Instance {
public:
    virtual void hello() { printf("Base class.\n"); }
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
    void point_connect(std::string i1, std::string i2, D* O::* member) {
        set<O,D*>(i1, member, dynamic_cast<D*>(instances[i2].get()));
    }

    template <class Connector, class... Args>
    void connect(Args &&... args) {
        Connector::_connect(*this, std::forward<Args>(args)...);
    }

};


template <class User, class Provider>
class UseProvide {
public:
    static void _connect(Model& model, std::string i1, std::string i2, Provider* User::* member) {
        model.point_connect<User, Provider>(i1, i2, member);
    }
};
