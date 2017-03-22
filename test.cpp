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
class I_Int : public Instance, public GetInt {
  public:
    int val;
    bool flag;

    I_Int(int val, bool flag) : val(val), flag(flag) {
        // printf("I_Int constructor: %d[%d].\n", val, (int)flag);
    }

    void _debug() override { printf("I_Int: %d[%d]", val, (int)flag); }

    int getInt() override { return val; }
};

class I_IntProxy : public Instance {
  public:
    GetInt* use;

    void _debug() override {
        if (use != nullptr)
            printf("I_IntProxy: %d", use->getInt());
        else
            printf("I_IntProxy: INVALID");
    }
};


// ######################################
//                  TEST
// ######################################
int main() {
    Assembly mymodel;

    mymodel.node<Instance>("I1");
    mymodel.node<I_IntProxy>("I2");
    mymodel.node<Array<I_Int> >("IArray", 5, [](int i) { return I_Int(2 * i, false); });
    mymodel.node<Array<I_IntProxy> >("IArray2", 5, [](int) { return I_IntProxy(); });
    mymodel.node<I_Int>("I4", 8, true);

    typedef UseProvide<I_IntProxy, GetInt> ProxyToGetInt;
    typedef UseProvideArray<I_IntProxy, I_Int, GetInt> ProxyArrayToGetInt;

    mymodel.connection<ProxyToGetInt>("I2", "I4", &I_IntProxy::use);
    mymodel.connection<ProxyArrayToGetInt>("IArray2", "IArray", &I_IntProxy::use);

    mymodel.instantiate();

    mymodel.set("I4", &I_Int::val, 37);

    mymodel.print_all();
}
