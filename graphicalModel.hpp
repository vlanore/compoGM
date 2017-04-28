#include <random>
#include "model.hpp"

/*

===================================================================================================
  INTERFACES
=================================================================================================*/
class Go : public Component {
  public:
    Go() { port("go", &Go::go); }
    virtual void go() = 0;
};

class Real {
  public:
    virtual double getValue() const = 0;
    virtual void setValue(double value) = 0;
};

class RandomNode : public Real {
    double clampedVal{0.0};

  public:
    RandomNode() = default;
    virtual void sample() = 0;
    void clamp(double val) { clampedVal = val; }
    double clampedValue() const { return clampedVal; }
    virtual bool isConsistent() const { return clampedVal == getValue(); }
};

/*

===================================================================================================
  Helper classes
=================================================================================================*/
// little object to encapsulate having a constant OR a pointer to Real
class RealProp {
    int mode{0};  // 0:unset, 1:constant, 2:pointer
    Real *ptr{nullptr};
    double value{0};

  public:
    RealProp() = default;
    explicit RealProp(double value) : mode(1), value(value) {}
    explicit RealProp(Real *ptr) : mode(2), ptr(ptr) {}
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

/*

===================================================================================================
  Graphical model nodes
=================================================================================================*/
class UnaryReal : public Component, public RandomNode {
  protected:
    RealProp param{};
    double value{0.0};

  public:
    std::string name{};

    // constructors
    UnaryReal() = delete;
    explicit UnaryReal(const std::string &name) : name(name) {
        port("paramConst", &UnaryReal::setParam<double>);
        port("paramPtr", &UnaryReal::setParam<Real *>);
        port("sample", &UnaryReal::sample);
    };

    // methods required from parent classes
    std::string _debug() const override {
        std::stringstream ss;
        ss << name << "(" << param.getValue() << "):" << value << "[" << clampedValue() << "]";
        return ss.str();
    }
    double getValue() const override { return value; }
    void setValue(double valuein) override { value = valuein; }

    // class-specific methods
    template <class... Args>
    void setParam(Args... args) {
        param = RealProp(std::forward<Args>(args)...);
    }
};

class Exponential : public UnaryReal {
  public:
    explicit Exponential(double value = 0.0) : UnaryReal("Exponential") { setValue(value); }

    void sample() override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::exponential_distribution<> d(param.getValue());
        setValue(d(gen));
    }
};

class Gamma : public UnaryReal {
  public:
    explicit Gamma() : UnaryReal("Gamma") {}

    void sample() override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::gamma_distribution<> d(param.getValue(), param.getValue());
        setValue(d(gen));
    }
};

class Poisson : public UnaryReal {
  public:
    explicit Poisson(double value = 0.0) : UnaryReal("Poisson") { setValue(value); }

    void sample() override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::poisson_distribution<> d(param.getValue());
        setValue(d(gen));
    }
};

class Product : public Component, public Real {
    RealProp a{};
    RealProp b{};

  public:
    Product() {
        port("aPtr", &Product::setA<Real *>);
        port("bPtr", &Product::setB<Real *>);
        port("bConst", &Product::setB<double>);
    }

    template <class... Args>
    void setA(Args... args) {
        a = RealProp(std::forward<Args>(args)...);
    }

    template <class... Args>
    void setB(Args... args) {
        b = RealProp(std::forward<Args>(args)...);
    }

    void setA(Real *ptr) { a = RealProp(ptr); }
    double getValue() const override { return a.getValue() * b.getValue(); }
    void setValue(double) override {
        std::cerr << "-- Warning! Trying to set a deterministic node!\n";
    }
    std::string _debug() const override {
        std::stringstream ss;
        ss << "Product(" << a.getValue() << "," << b.getValue() << "):" << getValue();
        return ss.str();
    }
};

/*

===================================================================================================
  Moves and MCMC-related things
=================================================================================================*/
class SimpleMove : public Go {
    RandomNode *target{nullptr};

  public:
    SimpleMove() { port("target", &SimpleMove::setTarget); }

    std::string _debug() const override { return "SimpleMove"; }
    void go() override { target->sample(); }

    void setTarget(RandomNode *targetin) { target = targetin; }
};

class Scheduler : public Go {
    std::vector<SimpleMove *> moves;

  public:
    Scheduler() { port("register", &Scheduler::registerMove); }

    void go() override {
        std::cout << "\n-- Scheduler started!\n"
                  << "-- Sampling everything!\n";
        for (auto &move : moves) {
            move->go();
        }
        std::cout << "-- Done.\n\n";
    }

    std::string _debug() const override {
        std::stringstream ss;
        ss << "Scheduler[" << moves.size() << "]";
        return ss.str();
    }

    void registerMove(SimpleMove *ptr) { moves.push_back(ptr); }
};

class MultiSample : public Go {
    std::vector<RandomNode *> nodes;

  public:
    MultiSample() { port("register", &MultiSample::registerNode); }

    void go() override {
        for (auto &node : nodes) {
            node->sample();
        }
    }

    std::string _debug() const override { return "MultiSample"; }

    void registerNode(RandomNode *ptr) { nodes.push_back(ptr); }
};

class RejectionSampling : public Go {
    std::vector<RandomNode *> observedData;
    Go *sampler{nullptr};
    int nbIter{0};

  public:
    explicit RejectionSampling(int iter = 5) {
        nbIter = iter;
        port("sampler", &RejectionSampling::setSampler);
        port("data", &RejectionSampling::addData);
    }

    void setSampler(Go *ptr) { sampler = ptr; }
    void addData(RandomNode *ptr) { observedData.push_back(ptr); }

    std::string _debug() const override { return "RejectionSampling"; }
    void go() override {
        std::cout << "-- Starting rejection sampling!\n";
        for (auto i = 0; i < nbIter; i++) {
            std::cout << "-- Iteration " << i << ". Sampling!\n";
            sampler->go();
        }
    }
};

/*

===================================================================================================
  Plates and custom connectors
=================================================================================================*/
class GraphicalModel : public Assembly {};
