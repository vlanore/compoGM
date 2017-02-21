#include <cstdio>
#include <map>
#include <vector>
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
class Assembly {
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


// ######################################
//            SPECIAL INSTANCES
// ######################################
template <class E>
class Array : public Instance {
    std::vector<E> vec;

public:
    Array(int size, E (*init)(int)) {
        for (int i=0; i<size; i++)
            vec.push_back(init(i));
    }

    void hello () {
        for (auto i: vec)
            i.hello();
    }
};


// ######################################
//               CONNECTORS
// ######################################
template <class User, class Provider>
class UseProvide {
public:
    static void _connect(Assembly& model, std::string i1, std::string i2, Provider* User::* member) {
        model.point_connect<User, Provider>(i1, i2, member);
    }
};
