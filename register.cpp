#include "register.h"


std::unordered_map<std::string, MyCallbackRegister::AlgCallback> MyCallbackRegister::eventMap;
boost::program_options::variables_map MyCallbackRegister::vm;


bool MyCallbackRegister::registerCallback(const std::string &name, const std::string &params, const std::function<void()> &callback)
{
    MyCallbackRegister::AlgCallback ac(name, params, callback);
    // event
    eventMap.emplace(name, ac);
    return true;
}

void MyCallbackRegister::setvm(boost::program_options::variables_map& vm_)
{
    MyCallbackRegister::vm = vm_;
}

void MyCallbackRegister::showAllRegisteredCallback()
{
    for(auto&p : MyCallbackRegister::eventMap){
        std::cout << p.first << ": ";
        printVec(p.second.params.begin(), p.second.params.end(), "", std::cout);
    }
}


void MyCallbackRegister::addArg(const std::string &name, const boost::any& anyArg)
{
    vm.emplace(name, boost::program_options::variable_value(anyArg, false));
}


void MyCallbackRegister::run(const std::string &name)
{
    auto it = eventMap.find(name);
    if(it == eventMap.end()){
        std::cout << name << " not found!" << std::endl;
    } else{
        auto& ac = it->second;
        bool allFound = true;
        for(auto& param:ac.params){
            if(!MyCallbackRegister::vm.count(param)){
                std::cout << "Require parameter " << param << std::endl;
                allFound = false;
            }
        }
        if(allFound){
            ac.callback();
        }
    }
}