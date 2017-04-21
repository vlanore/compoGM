#include "model.hpp"

// ######################################
//                INTERFACES
// ######################################
class GetInt {
  public:
    virtual int getInt() = 0;
};


// ######################################
//                INSTANCES
// ######################################
class MyInt : public Instance, public GetInt {
  public:
    int val;
    bool flag;

    MyInt(int val, bool flag) : val(val), flag(flag) {
        // printf("MyInt constructor: %d[%d].\n", val, (int)flag);
    }

    std::string _sdebug() override {
        std::stringstream result;
        result << "MyInt: " << val << "[" << flag << "]";
        return result.str();
        }

    int getInt() override { return val; }
};

class MyIntProxy : public Instance, public GetInt {
  public:
    GetInt* use;

    std::string _sdebug() override {
        std::stringstream result;
        result << "MyIntProxy: " << ((use!=nullptr)?use->getInt():-1) ;
        return result.str();
        }

    int getInt() final {
        if (use != nullptr) {
            return use->getInt();
        } else {
            return -1;
        }
    }
};

class MyReducer : public Instance {
  public:
    std::vector<GetInt*> use;

    std::string _sdebug() override {
        int sum = 0;
        for (auto i : use) {
            sum += i->getInt();
        }
        std::stringstream ss;
        ss << "SUM = " << sum;
        return ss.str();
    }
};

// ######################################
//                  TEST
// ######################################
int main() {
    Assembly mymodel;

    mymodel.node<Instance>("I1");
    mymodel.node<MyIntProxy>("I2");
    mymodel.node<Array<MyInt>>("IArray", 5, [](int i) { return MyInt(2 * i, false); });
    mymodel.node<Array<MyIntProxy>>("IArray2", 5, [](int) { return MyIntProxy(); });
    mymodel.node<Array<MyIntProxy>>("IArray3", 5, [](int) { return MyIntProxy(); });
    mymodel.node<MyInt>("I4", 8, true);
    mymodel.node<MyReducer>("Reducer");

    using ProxyToGetInt = UseProvide<MyIntProxy, GetInt>;
    using ProxyArrayToGetInt = UseProvideArray<MyIntProxy, MyInt, GetInt>;

    mymodel.connection<ProxyToGetInt>("I2", "I4", &MyIntProxy::use);
    mymodel.connection<ProxyArrayToGetInt>("IArray2", "IArray", &MyIntProxy::use);
    mymodel.connection<UseProvideArray<MyIntProxy, MyIntProxy, GetInt>>("IArray3", "IArray2",
                                                                        &MyIntProxy::use);

    mymodel.connection<MultiUseArray<MyReducer, MyIntProxy, GetInt>>("Reducer", "IArray3", &MyReducer::use);

    mymodel.instantiate();

    mymodel.print_all();

    mymodel.to_dot();
}
