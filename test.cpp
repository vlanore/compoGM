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

    I_Int(int val) : val(val) {}

    void hello() override {
        printf("I_Int(%d)", val);
    }

    int getInt() override { return val; }
};

class I_IntProxy : public Instance {
public:
    GetInt* use;

    void hello() override {
        if (use != nullptr)
            printf("I_IntProxy(%d)", use->getInt());
        else
            printf("I_IntProxy(%p)", use);
    }
};


// ######################################
//                  TEST
// ######################################
int main() {
    Assembly mymodel;

    mymodel.instantiate <Instance> ("I1");
    mymodel.instantiate <I_IntProxy> ("I2");
    mymodel.instantiate <Array<I_Int> > ("IArray", 5, [](int i) { return I_Int(2*i); });
    mymodel.instantiate <Array<I_IntProxy> > ("IArray2", 5, [](int) { return I_IntProxy(); });
    mymodel.instantiate <I_Int> ("I4", 8);

    mymodel.connect<UseProvide<I_IntProxy, GetInt> >("I2", "I4", &I_IntProxy::use);
    mymodel.connect<UseProvideArray<I_IntProxy, I_Int, GetInt> >("IArray2", "IArray", &I_IntProxy::use);

    mymodel.set<I_Int, int>("I4", &I_Int::val, 37);

    mymodel.print_all();

    return 0;
}
