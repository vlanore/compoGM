#include <mpi.h>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

template <typename T>
struct _Type {};


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
class Assembly;

class _Connection {
    std::function<void(Assembly&)> _connector;

  public:
    template <class T, class... Args>
    _Connection(_Type<T>, Args&&... args)
        : _connector([=](Assembly& a) { T::_connect(a, args...); }) {}

    void _connect(Assembly& a) { _connector(a); }
};

class _Node {
    std::function<std::unique_ptr<Instance>()> _constructor;

  public:
    std::string name;

    template <class T, class... Args>
    _Node(_Type<T>, std::string name, Args&&... args)
        : _constructor(
              [=]() { return std::unique_ptr<Instance>(dynamic_cast<Instance*>(new T(args...))); }),
          name(name) {}

    std::unique_ptr<Instance> _instantiate() { return _constructor(); }
};


// ######################################
//                 ASSEMBLY
// ######################################
class Assembly {
  private:
    // Model
    std::vector<_Node> nodes;
    std::vector<_Connection> connections;

  public:
    Assembly() {
#ifdef USE_MPI
        MPI_Init(NULL, NULL);

        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        char processor_name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(processor_name, &name_len);

        printf(
            "Hello world from processor %s, rank %d"
            " out of %d processors\n",
            processor_name, world_rank, world_size);

        MPI_Finalize();
#endif
    }

    // Actual instances
    std::map<std::string, std::unique_ptr<Instance> > instances;

    void print_all() {
        for (auto& i : instances) {
            printf("%s: ", i.first.c_str());
            i.second->_debug();
            std::cout << "." << std::endl;
        }
    }

    template <class Instance_Type, class... Args>
    void node(std::string name, Args&&... args) {
        nodes.emplace_back(_Type<Instance_Type>(), name, std::forward<Args>(args)...);
    }

    template <class Connector, class... Args>
    void connection(Args&&... args) {
        connections.emplace_back(_Type<Connector>(), std::forward<Args>(args)...);
    }

    void instantiate() {
        for (auto& n : nodes) {
            instances.emplace(n.name, n._instantiate());
        }
        for (auto& n : connections) {
            n._connect(*this);
        }
    }

    template <class I, typename V>
    void set(std::string name, V I::*member, V value) {
        dynamic_cast<I*>(instances[name].get())->*member = value;
    }

    template <class O, class D>
    void point_connect(std::string i1, std::string i2, D* O::*member) {
        set<O, D*>(i1, member, dynamic_cast<D*>(instances[i2].get()));
    }
};


class MPI_Assembly : public Assembly {
  public:
    MPI_Assembly() {}
};


// ######################################
//            SPECIAL INSTANCES
// ######################################
template <class E>
class Array : public Instance {
  public:
    std::vector<E> vec;

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
template <class User, class Interface>
class UseProvide {
  public:
    static void _connect(Assembly& model, std::string i1, std::string i2,
                         Interface* User::*member) {
        model.point_connect<User, Interface>(i1, i2, member);
    }
};

template <class User, class Provider, class Interface>
class UseProvideArray {
  public:
    static void _connect(Assembly& model, std::string i1, std::string i2,
                         Interface* User::*member) {
        auto ptrUser = dynamic_cast<Array<User>*>(model.instances[i1].get());
        auto ptrProvider = dynamic_cast<Array<Provider>*>(model.instances[i2].get());

        for (unsigned int i = 0; i < ptrUser->vec.size(); i++) {
            ptrUser->vec[i].*member = &ptrProvider->vec[i];
        }
    }
};

template <class User, class Provider, class Interface>
class MultiUseArray {
public:
    static void _connect(Assembly& model, std::string i1, std::string i2, std::vector<Interface*> User::*member) {

        auto ptrUser = dynamic_cast<User*>(model.instances[i1].get());
        auto ptrProvider = dynamic_cast<Array<Provider>*>(model.instances[i2].get());

        for (unsigned int i = 0; i < ptrProvider->vec.size(); i++) {
            (ptrUser->*member).push_back(&ptrProvider->vec[i]);
        }

    }
};
