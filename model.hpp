#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>


// ######################################
//                INSTANCES
// ######################################
class Instance {
public:
    virtual void _debug() { printf("Empty instance"); }
};


// ######################################
//                  MODEL
// ######################################
template < class T >
class Node {
    std::function<std::shared_ptr<T>()> _constructor;

    // template < class... Args >
    // std::shared_ptr<T> _factory(Args... args) { return std::make_shared<T>(std::forward<Args>(args)...); }

public:
    template < class... Args >
    Node(Args&&... args): _constructor( [=]() { return std::make_shared<T>(args...); } )
    {}

    std::shared_ptr<Instance> instantiate() { return _constructor(); }
};


// ######################################
//                 ASSEMBLY
// ######################################
class Assembly {
  public:
    std::map< std::string, std::shared_ptr< Instance > > instances;

    void print_all() {
        for (auto i : instances) {
            printf("%s: ", i.first.c_str());
            i.second->_debug();
            std::cout << "." << std::endl;
        }
    }

    template < class T, class... Args >
    void instantiate(std::string name, Args&&... args) {
        instances.emplace(name, std::make_shared< T >(std::forward< Args >(args)...));
    }

    template < class I, typename V >
    void set(std::string name, V I::*member, V value) {
        dynamic_cast< I* >(instances[name].get())->*member = value;
    }

    template < class O, class D >
    void point_connect(std::string i1, std::string i2, D* O::*member) {
        set< O, D* >(i1, member, dynamic_cast< D* >(instances[i2].get()));
    }

    template < class Connector, class... Args >
    void connect(Args&&... args) {
        Connector::_connect(*this, std::forward< Args >(args)...);
    }
};


// ######################################
//            SPECIAL INSTANCES
// ######################################
template < class E >
class Array : public Instance {
  public:
    std::vector< E > vec;

    Array(int size, E (*init)(int)) {
        for (int i = 0; i < size; i++) vec.push_back(init(i));
    }

    void _debug() {
        printf("Array : [");
        for (auto i : vec) {
            i._debug();
            std::cout << ", ";
        }
        std::cout << "]";
    }
};


// ######################################
//               CONNECTORS
// ######################################
template < class User, class Interface >
class UseProvide {
  public:
    static void _connect(Assembly& model, std::string i1, std::string i2,
                         Interface* User::*member) {
        model.point_connect< User, Interface >(i1, i2, member);
    }
};

template < class User, class Provider, class Interface >
class UseProvideArray {
  public:
    static void _connect(Assembly& model, std::string i1, std::string i2,
                         Interface* User::*member) {
        auto ptrUser = dynamic_cast< Array< User >* >(model.instances[i1].get());
        auto ptrProvider = dynamic_cast< Array< Provider >* >(model.instances[i2].get());

        for (unsigned int i = 0; i < ptrUser->vec.size(); i++) {
            ptrUser->vec[i].*member = &ptrProvider->vec[i];
        }
    }
};
