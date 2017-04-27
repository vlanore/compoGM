#include "model.hpp"

// interface to get a Real
class Real_I {
  public:
    virtual double getValue() const = 0;
};

// little object to encapsulate having a constant OR a pointer to Real_I
class RealProp {
    int mode{0};  // 0:unset, 1:constant, 2:pointer
    Real_I *ptr{nullptr};
    double value{0};

  public:
    RealProp() = default;
    explicit RealProp(double value) : mode(1), value(value) {}
    explicit RealProp(Real_I *ptr) : mode(2), ptr(ptr) {}
    double getValue() const {
        if (mode == 1) {
            return value;
        } else if (mode == 2) {
            return ptr->getValue();
        } else {
            std::cerr << "Property is no set!\n";
            exit(1);
        }
    }
};

class UnaryReal_C : public Component, public Real_I {
    RealProp param{};
    double value{0.0};

  public:
    std::string name{};

    UnaryReal_C() = delete;
    explicit UnaryReal_C(const std::string &name) : name(name) {
        port("paramConst", &UnaryReal_C::setLambda<double>);
        port("paramPtr", &UnaryReal_C::setLambda<Real_I *>);
    };

    std::string _debug() const override {
        std::stringstream ss;
        ss << name << "(" << param.getValue() << "):" << value;
        return ss.str();
    }
    double getValue() const override { return value; }
    void setValue(double valuein) { value = valuein; }

    template <class... Args>
    void setLambda(Args... args) {
        param = RealProp(std::forward<Args>(args)...);
    }
};

class Exponential_C : public UnaryReal_C {
  public:
    explicit Exponential_C(double value = 0.0) : UnaryReal_C("Exponential") { setValue(value); }
};

class Gamma_C : public UnaryReal_C {
  public:
    explicit Gamma_C() : UnaryReal_C("Gamma") {}
};

class Product_C : public Component, public Real_I {
    RealProp a{};
    RealProp b{};

  public:
    Product_C() {
        port("aPtr", &Product_C::setA<Real_I *>);
        port("bPtr", &Product_C::setB<Real_I *>);
        port("bConst", &Product_C::setB<double>);
    }

    template <class... Args>
    void setA(Args... args) {
        a = RealProp(std::forward<Args>(args)...);
    }

    template <class... Args>
    void setB(Args... args) {
        b = RealProp(std::forward<Args>(args)...);
    }

    void setA(Real_I *ptr) { a = RealProp(ptr); }
    double getValue() const override { return a.getValue() * b.getValue(); }
    std::string _debug() const override {
        std::stringstream ss;
        ss << "Product(" << a.getValue() << "," << b.getValue() << "):" << getValue();
        return ss.str();
    }
};

int main() {
    Assembly model;

    model.component<Exponential_C>("Sigma");
    model.property("Sigma", "paramConst", 1.0);

    model.component<Exponential_C>("Theta");
    model.property("Theta", "paramConst", 1.0);

    model.component<Array<Gamma_C, 5>>("Omega");
    model.connect<Map<Real_I>>("Theta", "Omega", "paramPtr");

    model.component<Array<Product_C, 5>>("rate");
    model.connect<ArrayOneToOne<Real_I>>("rate", "aPtr", "Omega");
    model.connect<ArrayOneToOne<Real_I>>("rate", "bPtr", "Omega");

    model.instantiate();
    model.print_all();
}
