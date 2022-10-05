//
//  main.cpp
//  LevelSetMethod2
//
//  Created by kumada on 2012/11/15.
//  Copyright (c) 2012年 kumada. All rights reserved.
//

#if(UNIT_TEST)
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_CASE(TEST_main)
{
        std::cout << "main\n";
}

#else // UNIT_TEST

#include "CommandLineInterface.h"
#include <boost/program_options.hpp>
#include <iostream>
namespace program_options = boost::program_options;

int main(int argc, char * argv[])
{
    std::cout << argv[0] << std::endl;
        try {
                lsm::CommandLineInterface::execute_level_set_method(argc, argv);
        } catch (std::exception& error) {
                std::cout << ">> " << error.what() << std::endl;
        }
        return 0;
}
#endif // UNIT_TEST
