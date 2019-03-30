#pragma once

#include <string> 
#include <functional>
#include <unordered_map>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/any.hpp>
#include "util.h"

namespace MyCallbackRegister
{
struct AlgCallback{
	std::string name;
	std::vector<std::string> params;
	std::function<void()> callback;
	AlgCallback(const std::string &name, const std::string &paramString, std::function<void()> callback)
			:name(name), callback(callback){
		std::istringstream ss(paramString);
		std::string token;
		while(ss>>token){
			params.push_back(token);
		}
	}
};
    
extern std::unordered_map<std::string, AlgCallback> eventMap;
extern boost::program_options::variables_map vm;

bool registerCallback(const std::string &name, const std::string &params, const std::function<void()> &callback);

void showAllRegisteredCallback();
    
void setvm(boost::program_options::variables_map& vm_);

void addArg(const std::string &name, const boost::any& anyArg);

void run(const std::string &name);

template<typename T>
T algAs(const std::string& name){
    assert(vm.find(name)!=vm.end());
    return vm[name].as<T>();
}
};
